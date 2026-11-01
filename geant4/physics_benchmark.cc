#include "G4Alpha.hh"
#include "G4DynamicParticle.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
double absorptionCrossSection(double,double,double,double);
int main(int argc,char** argv) {
    try {
        std::ofstream out; out.exceptions(std::ios::failbit|std::ios::badbit);
        out.open(argc>1?argv[1]:"data/nuclear_physics_benchmark.csv");
        out<<"E_lab_MeV,E_per_nucleon_MeV,semiclassical_mb,grid_difference_mb,timestep_difference_mb,aperture_difference_mb,shen_mb,glauber_gribov_mb,shen_relative_difference,gg_relative_difference\n"<<std::setprecision(12);
        auto* alpha=G4Alpha::AlphaDefinition();
        G4IonsShenCrossSection shen; G4ComponentGGNuclNuclXsc gg;
        std::cout<<G4Version<<"\nReaction cross sections (mb); model differences are NOT test failures\n";
        for(double energy:{5.,10.,25.,35.,50.,75.,100.}) {
            double coarse=absorptionCrossSection(energy,.1,20.,.1875);
            double fine=absorptionCrossSection(energy,.05,20.,.1875);
            double timeFine=absorptionCrossSection(energy,.05,20.,.09375);
            double wide=absorptionCrossSection(energy,.05,30.,.1875);
            // 1 fm^2 = 10 mb. Native models use physical projectile/target masses.
            double classical=10*fine;
            G4DynamicParticle projectile(alpha,G4ThreeVector(0,0,1),energy*MeV);
            double s=shen.GetIsoCrossSection(&projectile,79,197)/millibarn;
            double g=gg.GetInelasticIsotopeCrossSection(alpha,energy*MeV,79,197)/millibarn;
            if(!std::isfinite(s) || !std::isfinite(g) || s<0 || g<0) throw std::runtime_error("Invalid native Geant4 cross section");
            double grid=10*std::abs(fine-coarse),time=10*std::abs(fine-timeFine),tail=10*std::abs(fine-wide);
            out<<energy<<','<<energy/4<<','<<classical<<','<<grid<<','<<time<<','<<tail<<','<<s<<','<<g<<','<<(s>0?(classical-s)/s:NAN)<<','<<(g>0?(classical-g)/g:NAN)<<'\n';
            std::cout<<energy<<" MeV: classical="<<classical<<", Shen="<<s<<", GG="<<g<<"; grid/dt/aperture="<<grid<<'/'<<time<<'/'<<tail<<" mb\n"<<std::flush;
            if(grid>std::max(5.,.02*classical) || time>std::max(.5,.002*classical) || tail>.1) throw std::runtime_error("Absorption integral not numerically converged");
        }
        out.close();
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
