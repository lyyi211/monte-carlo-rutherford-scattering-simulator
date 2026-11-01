// A real Geant4 geometry/Qt viewer for independently integrated single-nucleus paths.
// Display mapping only: 1 physical fm is represented by 1 Geant4 mm.
// No material transport is run, and no solid boundary enters the equations of motion.
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4Geantino.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VUserVisAction.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include <array>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <filesystem>
#include <cmath>
#include <QApplication>
#include <QWidget>
#include <QTimer>
#include <QPixmap>

class NuclearScene : public G4VUserDetectorConstruction {
public:
    G4VPhysicalVolume* Construct() override {
        auto* vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        auto* world=new G4LogicalVolume(new G4Box("DisplayWorld",110*mm,110*mm,110*mm),vacuum,"DisplayWorld");
        world->SetVisAttributes(G4VisAttributes::GetInvisible());
        auto* placed=new G4PVPlacement(nullptr,{},world,"DisplayWorld",nullptr,false,0,true);
        // Fermi c is a half-density radius, not a hard scattering boundary.
        auto* core=new G4LogicalVolume(new G4Sphere("FermiHalfDensity",0,6.38*mm,0,360*deg,0,180*deg),vacuum,"Au197HalfDensityMarker");
        auto* coreStyle=new G4VisAttributes(G4Colour(1.,.75,.1,.55));
        coreStyle->SetForceSolid(true); core->SetVisAttributes(coreStyle);
        new G4PVPlacement(nullptr,{},core,"Au197HalfDensityMarker",world,false,0,true);
        auto* shell=new G4LogicalVolume(new G4Sphere("DiffuseReference",8.146*mm,8.166*mm,0,360*deg,0,180*deg),vacuum,"WoodsSaxonRadiusMarker");
        auto* shellStyle=new G4VisAttributes(G4Colour(.3,.7,1.,.25));
        shellStyle->SetForceWireframe(true);shell->SetVisAttributes(shellStyle);
        new G4PVPlacement(nullptr,{},shell,"WoodsSaxonRadiusMarker",world,false,0,true);
        return placed;
    }
};
class DisplayPhysics : public G4VUserPhysicsList {
    void ConstructParticle() override {G4Geantino::GeantinoDefinition();}
    void ConstructProcess() override {AddTransportation();}
    void SetCuts() override {SetCutsWithDefault();}
};
class NuclearPaths : public G4VUserVisAction {
    std::map<int,G4Polyline> paths;
public:
    NuclearPaths(const std::string& filename,int model) {
        std::ifstream in(filename);
        if(!in) throw std::runtime_error("Missing trajectory comparison CSV");
        std::string line; std::getline(in,line);
        while(std::getline(in,line)) {
            std::stringstream row(line);std::string word;std::vector<double> v;
            while(std::getline(row,word,',')) v.push_back(std::stod(word));
            if(v.size()!=11) throw std::runtime_error("Invalid trajectory CSV");
            if(static_cast<int>(v[1])==model) paths[static_cast<int>(v[0])].push_back({v[7]*mm,v[8]*mm,0});
        }
        if(paths.size()!=6) throw std::runtime_error("Expected six trajectories for selected model");
    }
    void Draw() override {
        auto* vis=G4VVisManager::GetConcreteInstance();
        const std::array<G4Colour,6> colors{G4Colour(1,1,1),G4Colour(1,.35,.25),G4Colour(.2,1,.5),G4Colour(.3,.65,1),G4Colour(1,.8,.2),G4Colour(1,.4,1)};
        int i=0;
        for(auto& entry:paths) {
            G4VisAttributes style(colors[i++]);style.SetLineWidth(2);
            entry.second.SetVisAttributes(style);vis->Draw(entry.second);
        }
        G4Text title("Alpha + Au-197: independent Geant4 RKF45 paths",G4Point3D(-.94,.90,0));
        title.SetScreenSize(18);title.SetVisAttributes(G4VisAttributes(G4Colour(1,1,1)));vis->Draw2D(title);
        G4Text caption("50 MeV | c = 6.38 fm | display: 1 mm = 1 fm | no material transport",G4Point3D(-.94,-.90,0));
        caption.SetScreenSize(14);caption.SetVisAttributes(G4VisAttributes(G4Colour(1,1,1)));vis->Draw2D(caption);
    }
};
int main(int argc,char** argv) {
    try {
        std::string input=argc>1?argv[1]:"data/trajectory_comparison.csv";
        int model=argc>2?std::stoi(argv[2]):3;
        if(model<0 || model>3) throw std::runtime_error("Model must be 0..3");
        auto run=std::make_unique<G4RunManager>();
        run->SetUserInitialization(new NuclearScene);
        run->SetUserInitialization(new DisplayPhysics);
        run->Initialize();
        int uiArgc=1; auto ui=std::make_unique<G4UIExecutive>(uiArgc,argv,"Qt");
        auto vis=std::make_unique<G4VisExecutive>();vis->Initialize();
        NuclearPaths paths(input,model);
        vis->RegisterRunDurationUserVisAction("NuclearPaths",&paths,G4VisExtent(-105*mm,105*mm,-105*mm,105*mm,-10*mm,10*mm));
        auto* commands=G4UImanager::GetUIpointer();
        for(const auto& command:{
            "/vis/open OGLSQt 1000x750+0+0",
            "/vis/viewer/set/background 0.035 0.05 0.09",
            "/vis/drawVolume",
            "/vis/scene/add/userAction NuclearPaths",
            "/vis/viewer/set/viewpointVector 0 0 1",
            "/vis/viewer/zoom 1.3",
            "/vis/viewer/flush"}) {
            if(commands->ApplyCommand(command)!=0) throw std::runtime_error(std::string("Viewer command failed: ")+command);
        }
        if(argc>3 && std::string(argv[3])=="--capture") {
            std::filesystem::create_directories("results/screenshots/frames");
            // Capture actual rendered Qt widgets, including the Geant4 UI.
            // The scene is rotated through Geant4's own camera commands.
            auto* timer=new QTimer(qApp);
            auto frame=std::make_shared<int>(0);
            QObject::connect(timer,&QTimer::timeout,[timer,frame,commands]() {
                QWidget* main=nullptr;
                for(auto* w:QApplication::topLevelWidgets())
                    if(w->isVisible() && (!main || w->width()*w->height()>main->width()*main->height())) main=w;
                if(!main) return;
                main->resize(1150,820);
                main->setWindowTitle("Alpha-Gold | Geant4 nuclear geometry");
                const int index=(*frame)++;
                if(index>=39) {
                    timer->stop();
                    std::ofstream("results/screenshots/capture_complete.txt")<<"39 real Qt frames captured\n";
                    return;
                }
                double theta=index<2?0.:55.;
                double phi=index<3?30.:(index-3)*10.;
                double zoom=index==0?1.3:(index==1?6.:3.5);
                double t=theta*3.141592653589793/180,p=phi*3.141592653589793/180;
                std::ostringstream direction;
                direction<<"/vis/viewer/set/viewpointVector "<<std::sin(t)*std::cos(p)<<' '<<std::sin(t)*std::sin(p)<<' '<<std::cos(t);
                commands->ApplyCommand(direction.str());
                commands->ApplyCommand("/vis/viewer/zoomTo "+std::to_string(zoom));
                commands->ApplyCommand("/vis/viewer/flush");
                std::string filename=index==0?"results/screenshots/geant4_wide.png":
                    index==1?"results/screenshots/geant4_close.png":
                    index==2?"results/screenshots/geant4_tilted.png":
                    "results/screenshots/frames/frame_"+std::to_string(index-3)+".png";
                QTimer::singleShot(350,[main,filename]() {
                    if(!main->grab().save(QString::fromStdString(filename)))
                        std::cerr<<"Screenshot save failed: "<<filename<<'\n';
                });
            });
            timer->start(1000);
        }
        ui->SessionStart();
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
