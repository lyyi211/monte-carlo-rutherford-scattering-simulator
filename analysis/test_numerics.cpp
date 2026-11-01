#define NUCLEAR_SCATTERING_NO_MAIN
#include "../nuclear_scattering.cpp"

void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
int main() {
    FermiPotential f;
    double worst=0;
    for(double E:{5.,10.,25.,50.,100.}) for(double b:{0.,.5,1.,5.,10.,20.}) {
        auto t=simulate(E,b,NuclearModel::Point,f);
        double exact=2*std::atan2(K,2*E*b)*180/PI;
        worst=std::max(worst,std::abs(t.thetaDeg-exact));
        require(t.escaped && std::abs(t.thetaDeg-exact)<5e-6,"Rutherford angle failed");
        double rmin=(K+std::sqrt(K*K+4*E*E*b*b))/(2*E);
        require(std::abs(t.minRadius-rmin)<.005,"Rutherford turning point failed");
    }
    for(double r:{.003,.013,.127,1.123,6.381,8.149,19.317,99.999,101.}) {
        double h=1e-6;
        require(std::abs((f.energy(r-h)-f.energy(r+h))/(2*h)-f.force(r))<2e-7,
                "Fermi force is not the potential gradient");
        require(std::abs((nuclearPotentialReal(r-h)-nuclearPotentialReal(r+h))/(2*h)-nuclearForceRadial(r))<2e-7,
                "Woods-Saxon force sign/gradient failed");
    }
    double central=4*PI*K*f.chargeDensityOverE(0)/ (3*Z_GOLD);
    require(std::abs(f.force(.001)/.001/central-1)<1e-5,"Central Gauss law failed");
    for(auto model:{NuclearModel::Point,NuclearModel::UniformSphere,NuclearModel::FermiCharge,NuclearModel::FermiWithNuclear}) {
        auto a=simulate(50,5,model,f,1000,.1875),b=simulate(50,5,model,f,200,.1875);
        require(std::abs(a.thetaDeg-b.thetaDeg)<.001,"Starting radius bias remains");
        auto fine=simulate(50,8.25,model,f,1000,.09375),coarse=simulate(50,8.25,model,f);
        require(std::abs(fine.thetaDeg-coarse.thetaDeg)<.01,"Timestep angle convergence failed");
        require(std::abs(fine.survival-coarse.survival)<1e-4,"Absorption quadrature failed");
    }
    std::cout<<"PASS: Rutherford, turning radii, potential gradients, central Gauss law, starting radius, timestep, absorption; max Rutherford error "<<worst<<" deg\n";
}
