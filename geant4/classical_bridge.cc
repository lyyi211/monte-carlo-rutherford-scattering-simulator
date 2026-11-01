#define NUCLEAR_SCATTERING_NO_MAIN
#include "../nuclear_scattering.cpp"
#include <array>
std::vector<std::array<double,4>> classicalTrace(double energy,double b,int model) {
    FermiPotential potential; std::vector<State> trace;
    auto result=simulate(energy,b,static_cast<NuclearModel>(model),potential,X_START,DT,&trace);
    if(!result.escaped) throw std::runtime_error("Classical trace did not escape");
    std::vector<std::array<double,4>> out;
    for(const auto& s:trace) out.push_back({s.x,s.y,s.vx,s.vy});
    return out;
}
double absorptionCrossSection(double energy,double db,double bmax,double dt) {
    FermiPotential potential;
    const int n=static_cast<int>(std::ceil(bmax/db)); db=bmax/n; double sum=0;
    for(int i=0;i<n;++i) {
        double b=(i+.5)*db;
        auto t=simulate(energy,b,NuclearModel::FermiWithNuclear,potential,X_START,dt);
        if(!t.escaped || !std::isfinite(t.survival)) throw std::runtime_error("Invalid absorption trajectory");
        sum+=2*PI*b*(1-t.survival)*db;
    }
    return sum;
}
