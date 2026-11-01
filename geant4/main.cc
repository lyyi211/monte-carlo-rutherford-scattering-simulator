// Independent integration of the SAME fixed-nucleus, nonrelativistic models.
// Geant4 supplies Cash-Karp RKF45 and ROOT output, not QGSP_BIC transport.
// Units of the custom ODE: fm, fm/c, MeV. No G4 geometry tolerances are used.
#include "G4AnalysisManager.hh"
#include "G4CashKarpRKF45.hh"
#include "G4EquationOfMotion.hh"
#include "G4MagneticField.hh"
#include "Randomize.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi=3.14159265358979323846, K=1.439964*2*79;
constexpr double mass=3727.38, hbar=197.3269804, dr=0.002;
const double sphereR=1.2*std::cbrt(197.), nuclearR=1.4*std::cbrt(197.);
const std::array<std::string,4> names{"point","uniform_sphere","fermi_coulomb","fermi_nuclear"};

// Independent quadrature: U = K [Q_fraction/r + integral_r^infinity s*rho(s) ds / N].
// Main C++ instead integrates a tabulated force inward on a 0.002 fm grid.
struct Potential {
    std::vector<double> charge, outer;
    double norm;
    static double density(double r) { return 1/(1+std::exp((r-6.38)/0.535)); }
    Potential() : charge(50001,0.), outer(50001,0.) {
        for (int i=1;i<=50000;++i) {
            double a=(i-1)*dr,b=i*dr,m=(a+b)/2;
            charge[i]=charge[i-1]+dr/6*(a*a*density(a)+4*m*m*density(m)+b*b*density(b));
        }
        norm=charge.back();
        for(int i=50000;i>0;--i) {
            double a=(i-1)*dr,b=i*dr,m=(a+b)/2;
            outer[i-1]=outer[i]+dr/6*(a*density(a)+4*m*density(m)+b*density(b));
        }
    }
    double interp(const std::vector<double>& v,double r) const {
        double u=r/dr; auto i=static_cast<std::size_t>(u);
        if(i+1>=v.size()) return v.back();
        return v[i]+(u-i)*(v[i+1]-v[i]);
    }
    double shape(double r) const { return 1/(1+std::exp((r-nuclearR)/0.52)); }
    double force(double r,int model) const {
        if(model==0 || r>=100) return K/(r*r);
        if(model==1) return r>=sphereR ? K/(r*r) : K*r/std::pow(sphereR,3);
        // Correct analytic central limit; avoid interpolating Q proportional to r^3 as r.
        double q=r<dr ? density(0)*r*r*r/3 : interp(charge,r);
        double f=K*q/(norm*r*r);
        if(model==3) { double w=shape(r); f-=185/0.52*w*(1-w); }
        return f;
    }
    double energy(double r,int model) const {
        if(model==0 || r>=100) return K/r;
        if(model==1) return r>=sphereR ? K/r : K/(2*sphereR)*(3-r*r/(sphereR*sphereR));
        double u=r<dr ? K*(outer[0]-density(0)*r*r/6)/norm :
                            K*(interp(charge,r)/r+interp(outer,r))/norm;
        return model==3 ? u-185*shape(r) : u;
    }
};
// G4EquationOfMotion requires a field object. The custom RHS below implements
// Newton's equation directly; this zero field is only an interface adapter.
struct NullField : G4MagneticField {
    void GetFieldValue(const G4double[4],G4double* b) const override { b[0]=b[1]=b[2]=0; }
};
struct Equation : G4EquationOfMotion {
    const Potential& p; int model;
    Equation(NullField* f,const Potential& pot,int m):G4EquationOfMotion(f),p(pot),model(m) {}
    void SetChargeMomentumMass(G4ChargeState,G4double,G4double) override {}
    void EvaluateRhsGivenB(const double y[],const double[],double d[]) const override {
        std::fill(d,d+8,0.);
        double r=std::hypot(y[0],y[1]);
        d[0]=y[2]; d[1]=y[3];
        if(r>1e-12) { double a=p.force(r,model)/(mass*r); d[2]=a*y[0]; d[3]=a*y[1]; }
        d[4]=model==3 ? 2*25*p.shape(r)/hbar : 0.; // optical depth
        d[7]=1.; // time in fm/c; never pass this state to G4 transportation
    }
};
using State=std::array<double,8>;
double direction(const State& y,double E,bool outgoing) {
    double r=std::hypot(y[0],y[1]), L=mass*(y[0]*y[3]-y[1]*y[2]);
    return std::atan2(-y[2]*L+K*y[1]/r,y[3]*L+K*y[0]/r)-
           std::atan2(-L*std::sqrt(2*E/mass),outgoing?K:-K);
}
struct Result {double theta,survival,energyError,momentumError;};
Result integrate(double E,double b,int model,const Potential& p,std::vector<State>* trace=nullptr) {
    if(!(std::isfinite(E) && E>0 && std::isfinite(b) && b>=0))
        throw std::runtime_error("Invalid incident energy or impact parameter");
    const double radius=std::hypot(100.,b),vinf=std::sqrt(2*E/mass);
    double vt=-b*vinf/radius,vr=-std::sqrt(2*(E-K/radius)/mass-vt*vt);
    State y{-100.,b,(-100*vr-b*vt)/radius,(b*vr-100*vt)/radius,0,0,0,0};
    double incoming=direction(y,E,false), L0=mass*(y[0]*y[3]-y[1]*y[2]);
    NullField field; Equation eq(&field,p,model); G4CashKarpRKF45 rk(&eq,8);
    double step=1.; bool escaped=false;
    if(trace) trace->push_back(y);
    double nextSample=0.1875;
    for(int count=0;count<2000000;++count) {
        if(trace) step=std::min(step,nextSample-y[7]);
        State deriv{},next{},error{}; eq.RightHandSide(y.data(),deriv.data());
        rk.Stepper(y.data(),deriv.data(),step,next.data(),error.data());
        double scaled=0;
        for(int k : {0,1,2,3,4}) {
            double atol=k<2?1e-9:1e-11;
            scaled=std::max(scaled,std::abs(error[k])/(atol+1e-10*std::max(std::abs(y[k]),std::abs(next[k]))));
        }
        if(!std::isfinite(scaled)) throw std::runtime_error("Nonfinite integration error");
        if(scaled<=1) {
            y=next;
            if(trace && std::abs(y[7]-nextSample)<1e-9) {trace->push_back(y); nextSample+=0.1875;}
            if(std::hypot(y[0],y[1])>=radius && y[0]*y[2]+y[1]*y[3]>0) {escaped=true;break;}
        }
        step=std::min(20.,step*std::clamp(0.9*std::pow(std::max(scaled,1e-16),-0.2),0.1,4.));
        if(step<1e-10) throw std::runtime_error("Adaptive step underflow");
    }
    if(!escaped) throw std::runtime_error("Trajectory did not escape");
    double finalE=mass/2*(y[2]*y[2]+y[3]*y[3])+p.energy(std::hypot(y[0],y[1]),model);
    double L=mass*(y[0]*y[3]-y[1]*y[2]);
    if(!std::isfinite(finalE) || std::abs(finalE/E-1)>1e-5 ||
       std::abs(L-L0)/std::max(1.,std::abs(L0))>1e-5)
        throw std::runtime_error("Geant4 integration conservation check failed");
    return {std::abs(std::remainder(direction(y,finalE,true)-incoming,2*pi))*180/pi,
            std::exp(-y[4]),std::abs(finalE/E-1),std::abs(L-L0)/std::max(1.,std::abs(L0))};
}
std::vector<std::string> split(const std::string& line) {
    std::stringstream in(line); std::vector<std::string> v; std::string token;
    while(std::getline(in,token,',')) v.push_back(token);
    return v;
}
int modelIndex(const std::string& name) {
    for(int i=0;i<4;++i) if(names[i]==name) return i;
    throw std::runtime_error("Unknown model: "+name);
}
}

#ifndef GEANT4_MODEL_NO_MAIN
int main(int argc,char** argv) {
    try {
        // Run from project root, or supply project directory and output path.
        if(argc>3) throw std::runtime_error("Usage: alpha_gold_geant4 [project_directory] [output.root]");
        std::filesystem::path project=argc>1?argv[1]:".";
        std::string output=argc>2?argv[2]:(project/"geant4/build/geant4_model_check.root").string();
        std::ifstream scans(project/"data/nuclear_model_scan.csv"), events(project/"data/detector_events.csv");
        if(!scans || !events) throw std::runtime_error("Run the classical simulation first; input CSV missing");
        Potential potential; auto* man=G4AnalysisManager::Instance();
        man->CreateNtuple("scan","Same-model Geant4 RKF45 vs classical RK4");
        man->CreateNtupleIColumn("model_id");
        for(const char* s:{"E","b","theta_deg","theta_cpp_deg","survival_probability","survival_cpp",
                           "energy_relative_error","angular_momentum_error"}) man->CreateNtupleDColumn(s);
        man->FinishNtuple();
        man->CreateNtuple("events","Matched impact parameters and ideal angular detector");
        for(const char* s:{"model_id","event_id","geometric_hit","registered","absorbed"}) man->CreateNtupleIColumn(s);
        for(const char* s:{"E_true","b","phi_true_deg","theta_true_deg","theta_measured_deg","E_measured_MeV",
                           "survival_probability","theta_cpp_deg","survival_cpp","registered_cpp","theta_measured_cpp","energy_measured_cpp"}) man->CreateNtupleDColumn(s);
        man->FinishNtuple();
        if(!man->OpenFile(output)) throw std::runtime_error("Cannot open ROOT output");
        CLHEP::HepRandom::setTheSeed(20261028);
        std::string line; std::getline(scans,line); int nscan=0,nevents=0; double maxAngle=0,maxSurvival=0;
        while(std::getline(scans,line)) {
            auto v=split(line); if(v.size()!=10) throw std::runtime_error("Invalid scan CSV");
            int model=modelIndex(v[0]); double E=std::stod(v[1]),b=std::stod(v[2]);
            if(v[8]!="1") throw std::runtime_error("Unescaped classical scan trajectory");
            Result r=integrate(E,b,model,potential);
            man->FillNtupleIColumn(0,0,model);
            double values[]{E,b,r.theta,std::stod(v[3]),r.survival,std::stod(v[7]),r.energyError,r.momentumError};
            for(int j=0;j<8;++j) man->FillNtupleDColumn(0,j+1,values[j]);
            man->AddNtupleRow(0); ++nscan;
            maxAngle=std::max(maxAngle,std::abs(r.theta-std::stod(v[3])));
            maxSurvival=std::max(maxSurvival,std::abs(r.survival-std::stod(v[7])));
        }
        std::getline(events,line);
        while(std::getline(events,line)) {
            auto v=split(line); if(v.size()!=12) throw std::runtime_error("Invalid events CSV");
            int model=modelIndex(v[0]); double E=std::stod(v[2]),b=std::stod(v[3]),phi=std::stod(v[4]);
            Result r=integrate(E,b,model,potential);
            bool hit=r.theta>=20 && r.theta<=160 && std::abs(phi)<=45;
            bool absorbed=G4UniformRand()>r.survival;
            bool registered=hit && !absorbed && G4UniformRand()<0.85;
            double thetaMeasured=registered?r.theta+G4RandGauss::shoot(0.,1.):-1.;
            double energyMeasured=registered?E*(1+0.01*G4RandGauss::shoot()):-1.;
            int ints[]{model,std::stoi(v[1]),hit,registered,absorbed};
            for(int j=0;j<5;++j) man->FillNtupleIColumn(1,j,ints[j]);
            double values[]{E,b,phi,r.theta,thetaMeasured,energyMeasured,r.survival,std::stod(v[5]),std::stod(v[6]),std::stod(v[9]),std::stod(v[10]),std::stod(v[11])};
            for(int j=0;j<12;++j) man->FillNtupleDColumn(1,j+5,values[j]);
            man->AddNtupleRow(1); ++nevents;
            maxAngle=std::max(maxAngle,std::abs(r.theta-std::stod(v[5])));
            maxSurvival=std::max(maxSurvival,std::abs(r.survival-std::stod(v[6])));
            if(nevents%15000==0) std::cout<<nevents<<" matched events checked\n"<<std::flush;
        }
        if(!man->Write() || !man->CloseFile()) throw std::runtime_error("ROOT output write failed");
        std::cout<<"Scan rows "<<nscan<<", events "<<nevents<<", max angle difference "<<maxAngle
                 <<" deg, max survival difference "<<maxSurvival<<'\n';
        // A numerical agreement test, not experimental validation of optical physics.
        if(nscan==0 || nevents==0 || maxAngle>0.05 || maxSurvival>0.001)
            throw std::runtime_error("Same-model numerical comparison failed");
        return 0;
    } catch(const std::exception& e) {std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}

#endif
