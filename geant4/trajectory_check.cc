#define GEANT4_MODEL_NO_MAIN
#include "main.cc"
#include <iomanip>
std::vector<std::array<double,4>> classicalTrace(double,double,int);
int main(int argc,char** argv) {
    try {
        std::ofstream out; out.exceptions(std::ios::failbit|std::ios::badbit);
        out.open(argc>1?argv[1]:"data/trajectory_comparison.csv");
        out<<"case_id,model_id,E,b,t,x_cpp,y_cpp,x_g4,y_g4,position_error,velocity_error\n"<<std::setprecision(15);
        Potential potential; double worst=0; int id=0;
        for(int m=0;m<4;++m) for(double b:{0.,1.,5.,8.03816104495,8.25,15.}) {
            auto cpp=classicalTrace(50.,b,m);
            std::vector<State> g4; integrate(50.,b,m,potential,&g4);
            auto n=std::min(cpp.size(),g4.size());
            if(n<100 || std::abs(static_cast<double>(cpp.size())-g4.size())>2) throw std::runtime_error("Trace coverage mismatch");
            for(std::size_t i=0;i<n;++i) {
                double dx=std::hypot(cpp[i][0]-g4[i][0],cpp[i][1]-g4[i][1]);
                double dv=std::hypot(cpp[i][2]-g4[i][2],cpp[i][3]-g4[i][3]);
                if(!std::isfinite(dx) || !std::isfinite(dv)) throw std::runtime_error("Nonfinite trace");
                worst=std::max(worst,dx);
                out<<id<<','<<m<<",50,"<<b<<','<<i*.1875<<','<<cpp[i][0]<<','<<cpp[i][1]<<','<<g4[i][0]<<','<<g4[i][1]<<','<<dx<<','<<dv<<'\n';
            }
            ++id;
        }
        out.close();
        std::cout<<"Compared "<<id<<" complete trajectories; max position error "<<worst<<" fm\n";
        if(worst>.02) throw std::runtime_error("Trajectory position agreement failed");
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
