#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

/* ==========================================================
Physical constants and literature-based parameters
============================================================*/

constexpr double PI = 3.14159265358979323846;
constexpr double E2 = 1.439964;                 // e^2/(4*pi*epsilon_0) [MeV fm]
constexpr double HBAR_C = 197.3269804;          // hbar*c [MeV fm]
constexpr double MASS_ALPHA = 3727.38;          // [MeV/c^2], with c=1
constexpr double Z_ALPHA = 2.0;
constexpr double Z_GOLD = 79.0;
constexpr double A_GOLD = 197.0;
constexpr double K = E2 * Z_ALPHA * Z_GOLD;

constexpr double R_UNIFORM_COEFF = 1.2;         // [fm], comparison model
constexpr double FERMI_C = 6.38;               // [fm], Au-197 charge half-density radius
constexpr double FERMI_A = 0.535;              // [fm], Au-197 charge diffuseness

// McFadden-Satchler (1966), global alpha optical potential at 24.7 MeV.
// These parameters are not a validated 5-100 MeV energy-dependent fit.
constexpr double WS_V0 = 185.0;                // [MeV], real volume depth
constexpr double WS_W0 = 25.0;                 // [MeV], imaginary volume depth
constexpr double WS_R0 = 1.4;                 // [fm], multiply A^(1/3)
constexpr double WS_A = 0.52;                 // [fm], diffuseness

constexpr double X_START = 100.0;             // [fm], exact Coulomb asymptotic matching
constexpr double DT = 0.1875;                  // [fm/c], check convergence separately
constexpr double B_MAX = 20.0;               // [fm], Monte Carlo beam aperture
constexpr int DETECTOR_EVENTS_PER_MODEL = 15000;
constexpr double DET_E = 50.0;               // [MeV], example detector run
constexpr double THETA_MIN_DEG = 20.0;
constexpr double THETA_MAX_DEG = 160.0;
constexpr double PHI_HALF_WIDTH_DEG = 45.0;   // angular sector, NOT a flat detector
constexpr double DET_EFFICIENCY = 0.85;
constexpr double SIGMA_THETA_DEG = 1.0;
constexpr double SIGMA_E_RELATIVE = 0.01;

const double R_UNIFORM = R_UNIFORM_COEFF * std::cbrt(A_GOLD);
const double R_WS = WS_R0 * std::cbrt(A_GOLD);

/* ==========================================================
State and trajectory results
============================================================*/

struct State {
    double x, y, vx, vy;
};

struct Derivative {
    double dx, dy, dvx, dvy;
};

enum class NuclearModel {
    Point,
    UniformSphere,
    FermiCharge,
    FermiWithNuclear
};

const char* modelName(NuclearModel model) {
    switch (model) {
        case NuclearModel::Point: return "point";
        case NuclearModel::UniformSphere: return "uniform_sphere";
        case NuclearModel::FermiCharge: return "fermi_coulomb";
        case NuclearModel::FermiWithNuclear: return "fermi_nuclear";
    }
    return "unknown";
}

struct Trajectory {
    double thetaDeg = std::numeric_limits<double>::quiet_NaN();
    double minRadius = std::numeric_limits<double>::infinity();
    double energyError = std::numeric_limits<double>::quiet_NaN();
    double angularMomentumError = std::numeric_limits<double>::quiet_NaN();
    double survival = 1.0;
    bool escaped = false;
    int steps = 0;
};

/* ==========================================================
Numerical electrostatic potential of a Fermi charge density
============================================================*/

class FermiPotential {
public:
    static constexpr double DR = 0.002;        // [fm]
    static constexpr double R_MAX = 100.0;    // [fm], negligible charge tail

    FermiPotential() {
        const std::size_t n = static_cast<std::size_t>(std::llround(R_MAX / DR));
        fractionalCharge.assign(n + 1, 0.0);
        radialForce.assign(n + 1, 0.0);
        potential.assign(n + 1, 0.0);
        densityShape.assign(n + 1, 0.0);

        // Density shape is 1/(1+exp((r-c)/a)).
        // Density is normalized by the integral of 4*pi*r^2*rho(r).
        for (std::size_t i = 0; i <= n; ++i) {
            const double r = i * DR;
            densityShape[i] = 1.0 / (1.0 + std::exp((r - FERMI_C) / FERMI_A));
        }

        for (std::size_t i = 1; i <= n; ++i) {
            const double r1 = (i - 1) * DR;
            const double r2 = i * DR;
            fractionalCharge[i] = fractionalCharge[i - 1] +
                DR / 6.0 * (r1*r1*densityShape[i - 1] +
                    4.0 * std::pow(0.5*(r1+r2), 2) /
                    (1.0 + std::exp((0.5*(r1+r2)-FERMI_C)/FERMI_A)) +
                    r2*r2*densityShape[i]);
        }
        const double normalization = fractionalCharge[n];
        shapeIntegral = normalization;
        for (std::size_t i = 0; i <= n; ++i) {
            fractionalCharge[i] /= normalization;
        }

        radialForce[0] = 0.0;
        for (std::size_t i = 1; i <= n; ++i) {
            const double r = i * DR;
            radialForce[i] = K * fractionalCharge[i] / (r*r);
        }

        // Integrate U(r) inward, with U(R_MAX)=K/R_MAX.
        potential[n] = K / R_MAX;
        for (std::size_t i = n; i > 0; --i) {
            potential[i - 1] = potential[i] +
                0.5 * DR * (radialForce[i - 1] + radialForce[i]);
        }
    }

    double enclosedFraction(double r) const {
        if (r >= R_MAX) return 1.0;
        return interpolate(fractionalCharge, std::max(0.0, r));
    }

    double force(double r) const {
        if (r >= R_MAX) return K / (r*r);
        return interpolate(radialForce, std::max(0.0, r));
    }

    double energy(double r) const {
        if (r >= R_MAX) return K / r;
        // Exact integral of the piecewise-linear force: F = -dU/dr.
        r = std::max(0.0, r);
        const auto i = static_cast<std::size_t>(r / DR);
        const double h = r - i*DR;
        return potential[i] - radialForce[i]*h -
               0.5*(radialForce[i+1]-radialForce[i])*h*h/DR;
    }

    // rho_charge(r) in units of elementary charge per fm^3.
    double chargeDensityOverE(double r) const {
        if (r >= R_MAX) return 0.0;
        const double normalization = shapeIntegral;
        return Z_GOLD / (4.0 * PI * normalization) /
               (1.0 + std::exp((r - FERMI_C) / FERMI_A));
    }

private:
    std::vector<double> fractionalCharge, radialForce, potential, densityShape;
    double shapeIntegral = 0.0;

    static double interpolate(const std::vector<double>& values, double r) {
        const double u = r / DR;
        const auto i = static_cast<std::size_t>(u);
        if (i + 1 >= values.size()) return values.back();
        const double t = u - i;
        return values[i] * (1.0 - t) + values[i + 1] * t;
    }
};

/* ==========================================================
Nuclear Woods-Saxon form factor and potential
============================================================*/

double woodsSaxonShape(double r) {
    return 1.0 / (1.0 + std::exp((r - R_WS) / WS_A));
}

double nuclearPotentialReal(double r) {
    return -WS_V0 * woodsSaxonShape(r);
}

double nuclearPotentialImaginary(double r) {
    return -WS_W0 * woodsSaxonShape(r);
}

// F_N=-dU_N/dr, where U_N=-V0*f(r).
double nuclearForceRadial(double r) {
    const double f = woodsSaxonShape(r);
    return -(WS_V0 / WS_A) * f * (1.0 - f);
}

double potentialEnergy(double r, NuclearModel model, const FermiPotential& fermi) {
    if (model == NuclearModel::Point) return K / r;
    if (model == NuclearModel::UniformSphere) {
        if (r >= R_UNIFORM) return K / r;
        return K * (3.0 - r*r / (R_UNIFORM*R_UNIFORM)) / (2.0 * R_UNIFORM);
    }
    const double U_C = fermi.energy(r);
    return model == NuclearModel::FermiWithNuclear ? U_C + nuclearPotentialReal(r) : U_C;
}

double radialForce(double r, NuclearModel model, const FermiPotential& fermi) {
    if (model == NuclearModel::Point) return K / (r*r);
    if (model == NuclearModel::UniformSphere) {
        return r >= R_UNIFORM ? K / (r*r) : K * r / std::pow(R_UNIFORM, 3);
    }
    const double F_C = fermi.force(r);
    return model == NuclearModel::FermiWithNuclear ? F_C + nuclearForceRadial(r) : F_C;
}

/* ==========================================================
RK4 integrator
============================================================*/

Derivative equations(const State& s, NuclearModel model, const FermiPotential& fermi) {
    const double r = std::hypot(s.x, s.y);
    if (r < 1.0e-12) {
        // A nonzero-impact-parameter point trajectory must never reach r=0.
        if (model == NuclearModel::Point) throw std::runtime_error("Point nucleus reached r=0");
        return {s.vx, s.vy, 0.0, 0.0};
    }
    const double scale = radialForce(r, model, fermi) / (MASS_ALPHA * r);
    return {s.vx, s.vy, scale * s.x, scale * s.y};
}

State addScaled(const State& s, const Derivative& d, double factor) {
    return {s.x + factor*d.dx, s.y + factor*d.dy,
            s.vx + factor*d.dvx, s.vy + factor*d.dvy};
}

State rk4Step(const State& s, double dt, NuclearModel model, const FermiPotential& fermi) {
    const Derivative k1 = equations(s, model, fermi);
    const Derivative k2 = equations(addScaled(s, k1, 0.5*dt), model, fermi);
    const Derivative k3 = equations(addScaled(s, k2, 0.5*dt), model, fermi);
    const Derivative k4 = equations(addScaled(s, k3, dt), model, fermi);
    return {s.x  + dt*(k1.dx  + 2*k2.dx  + 2*k3.dx  + k4.dx)/6.0,
            s.y  + dt*(k1.dy  + 2*k2.dy  + 2*k3.dy  + k4.dy)/6.0,
            s.vx + dt*(k1.dvx + 2*k2.dvx + 2*k3.dvx + k4.dvx)/6.0,
            s.vy + dt*(k1.dvy + 2*k2.dvy + 2*k3.dvy + k4.dvy)/6.0};
}

// Outside the charge table the force is exactly Coulomb (WS tail negligible).
// The repulsive Runge-Lenz vector A = v x L + K r_hat gives both
// asymptotic directions, removing finite-radius angular bias.
double asymptoticDirection(const State& s, double energy, bool outgoing) {
    const double r = std::hypot(s.x, s.y);
    const double L = MASS_ALPHA*(s.x*s.vy-s.y*s.vx);
    const double ax = s.vy*L + K*s.x/r;
    const double ay = -s.vx*L + K*s.y/r;
    return std::atan2(ay, ax) -
           std::atan2(-L*std::sqrt(2.0*energy/MASS_ALPHA), outgoing ? K : -K);
}

Trajectory simulate(double energy, double b, NuclearModel model, const FermiPotential& fermi,
                    double xStart = X_START, double dt = DT, std::vector<State>* trace = nullptr) {
    if (!(energy > 0.0 && b >= 0.0 && dt > 0.0 && xStart >= FermiPotential::R_MAX))
        throw std::runtime_error("Require E>0, b>=0, dt>0 and start outside the charge table");
    const double rInitial = std::hypot(xStart, b);
    const double kinetic = energy - potentialEnergy(rInitial, model, fermi);
    if (kinetic <= 0.0) throw std::runtime_error("Initial kinetic energy must be positive");
    const double vInitial = std::sqrt(2.0*kinetic/MASS_ALPHA);
    const double vInfinity = std::sqrt(2.0*energy/MASS_ALPHA);
    const double vt = -b*vInfinity/rInitial;
    const double vr2 = vInitial*vInitial-vt*vt;
    if (vr2 <= 0.0) throw std::runtime_error("Start radius is inside the turning point");
    const double vr = -std::sqrt(vr2);
    State s{-xStart, b, (-xStart*vr-b*vt)/rInitial,
                           (b*vr-xStart*vt)/rInitial};
    const double incomingAngle = asymptoticDirection(s, energy, false);
    const double L0 = MASS_ALPHA * (s.x*s.vy - s.y*s.vx);
    const double E0 = kinetic + potentialEnergy(rInitial, model, fermi);

    if(trace) {trace->clear(); trace->push_back(s);}
    Trajectory result;
    result.minRadius = rInitial;
    double imaginaryIntegral = 0.0;

    // Bound trajectories and very slow orbiting are not forced into a false angle.
    const int maxSteps = static_cast<int>(std::min(2000000.0,
        std::ceil((10.0*xStart / vInitial + 15000.0) / dt)));
    for (int i = 0; i < maxSteps; ++i) {
        const double rBefore = std::hypot(s.x, s.y);
        State next = rk4Step(s, dt, model, fermi);
        const double rAfter = std::hypot(next.x, next.y);
        result.minRadius = std::min({result.minRadius, rBefore, rAfter});
        if (model == NuclearModel::FermiWithNuclear) {
            imaginaryIntegral += 0.5*dt*
                (std::abs(nuclearPotentialImaginary(rBefore)) +
                 std::abs(nuclearPotentialImaginary(rAfter)));
        }
        s = next;
        if(trace) trace->push_back(s);
        result.steps = i + 1;
        if (rAfter >= rInitial && (s.x*s.vx + s.y*s.vy) > 0.0) {
            result.escaped = true;
            break;
        }
    }

    if (!result.escaped) return result;
    const double Ef = 0.5*MASS_ALPHA*(s.vx*s.vx + s.vy*s.vy) +
                      potentialEnergy(std::hypot(s.x, s.y), model, fermi);
    const double Lf = MASS_ALPHA*(s.x*s.vy - s.y*s.vx);
    result.thetaDeg = std::abs(std::remainder(
        asymptoticDirection(s, Ef, true)-incomingAngle, 2.0*PI)) * 180.0 / PI;
    result.energyError = std::abs(Ef - E0) / std::abs(E0);
    result.angularMomentumError = std::abs(L0) > 1e-12 ?
        std::abs(Lf-L0)/std::abs(L0) : std::abs(Lf-L0);

    // SEMICLASSICAL attenuation proxy, not a quantum optical-model S-matrix.
    // W=-25 f(r) MeV; probability survival exp[-2 int |W| dt / (hbar*c)].
    if (model == NuclearModel::FermiWithNuclear) {
        result.survival = std::exp(-2.0*imaginaryIntegral/HBAR_C);
    }
    return result;
}

/* ==========================================================
CSV files and Monte Carlo detector response
============================================================*/

#ifndef NUCLEAR_SCATTERING_NO_MAIN
int main(int argc, char** argv) {
    try {
        const bool quick = argc > 1 && std::string(argv[1]) == "--quick";
        const bool convergence = argc > 1 && std::string(argv[1]) == "--convergence";
        std::filesystem::create_directories("data");
        std::filesystem::create_directories("results");
        FermiPotential fermi;
        if (convergence) {
            std::ofstream out;
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out.open("data/nuclear_timestep_convergence.csv");
            out << "case,model,E,b,dt,theta_deg,energy_relative_error,escaped\n";
            out << std::setprecision(15);
            struct CheckCase { const char* name; NuclearModel model; double E, b; };
            const std::vector<CheckCase> cases{
                {"Coulomb outer", NuclearModel::Point, 10.0, 5.0},
                {"Fermi penetration", NuclearModel::FermiCharge, 50.0, 1.0},
                {"Nuclear interior", NuclearModel::FermiWithNuclear, 50.0, 1.0},
                {"Nuclear grazing", NuclearModel::FermiWithNuclear, 50.0, 8.25},
                {"Nuclear near-critical", NuclearModel::FermiWithNuclear, 50.0, 8.03816104495}
            };
            const std::vector<double> dts{1.5, 0.75, 0.375, 0.1875, 0.09375};
            for (const auto& c : cases) {
                for (double dt : dts) {
                    const auto result = simulate(c.E, c.b, c.model, fermi, X_START, dt);
                    out << c.name << ',' << modelName(c.model) << ',' << c.E
                        << ',' << c.b << ',' << dt << ',' << result.thetaDeg
                        << ',' << result.energyError << ',' << (result.escaped?1:0) << '\n';
                }
                std::cout << "Convergence case: " << c.name << "\n";
            }
            out.close();
            std::cout << "Saved data/nuclear_timestep_convergence.csv\n";
            return 0;
        }
        std::cout << std::setprecision(9)
                  << "Au-197: c=" << FERMI_C << " fm; a=" << FERMI_A
                  << " fm; Woods-Saxon R=" << R_WS << " fm\n"
                  << "Fermi U(0)=" << fermi.energy(0.0) << " MeV\n";

        std::ofstream profile;
        profile.exceptions(std::ios::failbit | std::ios::badbit);
        profile.open("data/fermi_radial_profile.csv");
        profile << "r,rho_fermi_e_per_fm3,rho_uniform_e_per_fm3,Q_fraction,"
                   "U_point,U_uniform,U_fermi,U_fermi_nuclear,U_nuclear_real,W_imag,"
                   "F_fermi,F_nuclear\n";
        profile << std::setprecision(12);
        for (int i = 0; i <= 1500; ++i) {
            const double r = i*0.02;
            const double UPoint = r == 0.0 ? std::numeric_limits<double>::quiet_NaN() : K/r;
            const double uniformDensity = r < R_UNIFORM ? Z_GOLD/(4.0/3.0*PI*std::pow(R_UNIFORM,3)) : 0.0;
            profile << r << ',' << fermi.chargeDensityOverE(r) << ',' << uniformDensity << ','
                    << fermi.enclosedFraction(r) << ',' << UPoint << ','
                    << potentialEnergy(r,NuclearModel::UniformSphere,fermi) << ','
                    << fermi.energy(r) << ','
                    << potentialEnergy(r,NuclearModel::FermiWithNuclear,fermi) << ','
                    << nuclearPotentialReal(r) << ',' << nuclearPotentialImaginary(r) << ','
                    << fermi.force(r) << ',' << nuclearForceRadial(r) << '\n';
        }

        const std::vector<NuclearModel> models{
            NuclearModel::Point, NuclearModel::UniformSphere,
            NuclearModel::FermiCharge, NuclearModel::FermiWithNuclear};
        const std::vector<double> energies = quick ?
            std::vector<double>{10.0, 50.0, 75.0} :
            std::vector<double>{5.0, 10.0, 25.0, 35.0, 50.0, 75.0, 100.0};
        const double db = quick ? 1.0 : 0.5;
        const int bSteps = static_cast<int>(std::round(B_MAX/db));

        std::ofstream scan;
        scan.exceptions(std::ios::failbit | std::ios::badbit);
        scan.open("data/nuclear_model_scan.csv");
        scan << "model,E,b,theta_deg,r_min,energy_relative_error,"
                "angular_momentum_error,survival_probability,escaped,steps\n";
        scan << std::setprecision(12);
        for (double E : energies) {
            std::cout << "Scanning E=" << E << " MeV...\n" << std::flush;
            for (std::size_t j=0; j<models.size(); ++j) {
                for (int i=0; i<=bSteps; ++i) {
                    double b=i*db;
                    Trajectory tr = simulate(E,b,models[j],fermi);
                    scan << modelName(models[j]) << ',' << E << ',' << b << ','
                         << tr.thetaDeg << ',' << tr.minRadius << ',' << tr.energyError << ','
                         << tr.angularMomentumError << ',' << tr.survival << ','
                         << (tr.escaped?1:0) << ',' << tr.steps << '\n';
                }
            }
        }
        scan.close();

        std::mt19937_64 gen(20261028);
        std::uniform_real_distribution<double> unif(0.0,1.0);
        std::normal_distribution<double> gaussian(0.0,1.0);
        const int nEvents=quick?3000:DETECTOR_EVENTS_PER_MODEL;
        std::ofstream events;
        events.exceptions(std::ios::failbit | std::ios::badbit);
        events.open("data/detector_events.csv");
        events << "model,event_id,E_true,b,phi_true_deg,theta_true_deg,"
                  "survival_probability,absorbed,geometric_hit,registered,"
                  "theta_measured_deg,E_measured_MeV\n";
        events << std::setprecision(12);
        std::ofstream summary;
        summary.exceptions(std::ios::failbit | std::ios::badbit);
        summary.open("data/detector_summary.csv");
        summary << "model,generated,valid_trajectory,geometric_hits,"
                   "survivors_geometric,registered,registration_fraction\n";
        for (std::size_t j=0;j<models.size();++j) {
            int nValid=0,nGeometric=0,nSurvivors=0,nRegistered=0;
            for (int i=0;i<nEvents;++i) {
                const double b=B_MAX*std::sqrt(unif(gen));
                const double phi=360.0*unif(gen)-180.0;
                const Trajectory tr=simulate(DET_E,b,models[j],fermi);
                const double theta=tr.thetaDeg;
                const double survival=tr.survival;
                const bool valid=tr.escaped && std::isfinite(theta) && std::isfinite(survival);
                if (valid) ++nValid;
                const bool geometric=valid && theta>=THETA_MIN_DEG &&
                    theta<=THETA_MAX_DEG && std::abs(phi)<=PHI_HALF_WIDTH_DEG;
                if (geometric) ++nGeometric;
                // A sampled absorption outcome is independent of detector efficiency.
                const bool absorbed=valid && unif(gen)>survival;
                if (geometric && !absorbed) ++nSurvivors;
                const bool registered=geometric && !absorbed && unif(gen)<DET_EFFICIENCY;
                if (registered) ++nRegistered;
                // -1 is a sentinel for non-registered measurements in ROOT/CSV.
                double thetaMeasured=-1.0;
                double eMeasured=-1.0;
                if (registered) {
                    thetaMeasured=theta+SIGMA_THETA_DEG*gaussian(gen);
                    eMeasured=DET_E*(1.0+SIGMA_E_RELATIVE*gaussian(gen));
                }
                events << modelName(models[j]) << ',' << i << ',' << DET_E << ',' << b << ','
                       << phi << ',' << theta << ',' << survival << ','
                       << (absorbed?1:0) << ',' << (geometric?1:0) << ','
                       << (registered?1:0) << ',' << thetaMeasured << ',' << eMeasured << '\n';
            }
            summary << modelName(models[j]) << ',' << nEvents << ',' << nValid << ','
                    << nGeometric << ',' << nSurvivors << ',' << nRegistered << ','
                    << static_cast<double>(nRegistered)/nEvents << '\n';
            std::cout << modelName(models[j]) << ": " << nRegistered << "/" << nEvents
                      << " registered\n";
        }
        profile.close(); events.close(); summary.close();
        std::cout << "Files saved in data/. 'absorbed' is only a semiclassical proxy.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}

#endif
