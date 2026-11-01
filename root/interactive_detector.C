// Run with: root -l root/interactive_detector.C (no -b or -q).
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TStyle.h>
void interactive_detector() {
    gStyle->SetOptStat(0);
    auto* file=TFile::Open("data/geant4_model_check.root");
    if(!file || file->IsZombie()) throw std::runtime_error("Missing data/geant4_model_check.root");
    auto* events=file->Get<TTree>("events");
    if(!events) throw std::runtime_error("Missing events tree");
    auto* canvas=new TCanvas("detector_gui","Alpha-Gold | ROOT interactive detector comparison",1200,850);
    canvas->Divide(2,2);
    const char* labels[]{"Point Coulomb","Uniform sphere","Fermi Coulomb","Fermi + Woods-Saxon"};
    for(int model=0;model<4;++model) {
        canvas->cd(model+1);
        auto* cpp=new TH1D(Form("cpp_%d",model),Form("%s, 50 MeV;Measured angle (deg);Registered events / 5 deg",labels[model]),36,0,180);
        auto* g4=new TH1D(Form("g4_%d",model),"",36,0,180);
        events->Draw(Form("theta_measured_cpp>>cpp_%d",model),Form("model_id==%d && registered_cpp",model),"goff");
        events->Draw(Form("theta_measured_deg>>g4_%d",model),Form("model_id==%d && registered",model),"goff");
        cpp->SetLineColor(kRed+1);g4->SetLineColor(kBlue+1);
        cpp->SetLineWidth(2);g4->SetLineWidth(2);
        cpp->SetMaximum(1.25*std::max(cpp->GetMaximum(),g4->GetMaximum()));
        cpp->Draw("E");g4->Draw("E same");
        auto* legend=new TLegend(.51,.72,.89,.88);
        legend->AddEntry(cpp,"C++ RK4","l");legend->AddEntry(g4,"Geant4 RKF45","l");legend->Draw();
    }
    canvas->Update();
}
