#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

// Run from project root. Input contains the SAME four single-nucleus models.
// This is a numerical cross-check, not experimental validation of the potential.
void analyze_geant4(const char* filename="geant4/build/geant4_model_check.root") {
    gSystem->mkdir("results",true);
    gStyle->SetOptStat(0);
    ROOT::RDataFrame rawScan("scan",filename), rawEvents("events",filename);
    auto scan=rawScan.Define("angle_difference","theta_deg-theta_cpp_deg")
                     .Define("survival_difference","survival_probability-survival_cpp")
                     .Define("abs_angle_difference","abs(angle_difference)")
                     .Define("abs_survival_difference","abs(survival_difference)");
    auto events=rawEvents.Define("angle_difference","theta_true_deg-theta_cpp_deg")
                         .Define("abs_angle_difference","abs(angle_difference)")
                         .Define("abs_survival_difference","abs(survival_probability-survival_cpp)")
                         .Define("angular_residual","theta_measured_deg-theta_true_deg")
                         .Define("energy_residual","(E_measured_MeV/E_true-1.)/0.01")
                         .Define("expected_registration","geometric_hit*survival_probability*0.85");
    if(*scan.Count()==0 || *events.Count()==0) throw std::runtime_error("Empty input");
    const char* labels[]{"Point Coulomb","Uniform sphere","Fermi Coulomb","Fermi + Woods-Saxon"};
    TCanvas angles("angles","Matched angular distributions",1200,900); angles.Divide(2,2);
    TCanvas checks("checks","Numerical differences",1200,900); checks.Divide(2,2);
    TCanvas direct("direct","C++ versus Geant4 detector response",1200,1600); direct.Divide(2,4);
    TCanvas response("response","Matched detector",1200,900); response.Divide(2,2);
    // Keep lazy histogram results alive until the canvases are saved.
    std::vector<ROOT::RDF::RResultPtr<TH1D>> keep;
    for(int m=0;m<4;++m) {
        auto sample=events.Filter("model_id == "+std::to_string(m));
        auto generated=sample.Histo1D({("g"+std::to_string(m)).c_str(),labels[m],90,0,180},"theta_true_deg");
        auto reference=sample.Histo1D({("r"+std::to_string(m)).c_str(),"C++ RK4",90,0,180},"theta_cpp_deg");
        auto registered=sample.Filter("registered == 1");
        auto detected=registered.Histo1D({("d"+std::to_string(m)).c_str(),labels[m],90,0,180},"theta_measured_deg");
        auto expected=sample.Histo1D({("w"+std::to_string(m)).c_str(),"Expected before smearing",90,0,180},
                                     "theta_true_deg","expected_registration");
        angles.cd(m+1); generated->SetLineColor(kBlue); reference->SetLineColor(kRed);
        generated->GetXaxis()->SetTitle("Asymptotic angle (deg)"); generated->GetYaxis()->SetTitle("Events");
        generated->DrawCopy("hist"); reference->SetLineStyle(2); reference->DrawCopy("hist same");
        auto* legend=new TLegend(.5,.72,.88,.88); legend->AddEntry(generated.GetPtr(),"Geant4 RKF45","l");
        legend->AddEntry(reference.GetPtr(),"C++ RK4","l"); legend->Draw();
        response.cd(m+1); expected->SetLineColor(kBlue); detected->SetLineColor(kRed);
        expected->SetTitle(labels[m]); expected->GetXaxis()->SetTitle("Angle (deg)");
        expected->DrawCopy("hist"); detected->DrawCopy("E same");
        auto* leg=new TLegend(.39,.72,.89,.88); leg->AddEntry(expected.GetPtr(),"Expected, before smearing","l");
        leg->AddEntry(detected.GetPtr(),"Registered, smeared","l"); leg->Draw();
        const auto n=*sample.Count(), hits=*sample.Filter("geometric_hit == 1").Count(), count=*registered.Count();
        const double mean=*sample.Sum<double>("expected_registration");
        auto variance=sample.Define("variance","expected_registration*(1-expected_registration)").Sum<double>("variance");
        double pull=(count-mean)/std::sqrt(*variance);
        std::cout<<labels[m]<<": generated="<<n<<", geometric="<<hits<<", registered="<<count
                 <<", expected="<<mean<<", registration pull="<<pull<<'\n';
        auto bad=sample.Filter("registered && (!geometric_hit || absorbed)").Count();
        auto missing=sample.Filter("!registered && (theta_measured_deg != -1 || E_measured_MeV != -1)").Count();
        if(*bad || *missing || std::abs(pull)>6) throw std::runtime_error("Detector bookkeeping/statistics failed");
        if(count>100) {
            double angularSigma=*registered.StdDev<double>("angular_residual");
            double energySigma=*registered.StdDev<double>("energy_residual");
            if(std::abs(angularSigma-1)>0.2 || std::abs(energySigma-1)>0.2)
                throw std::runtime_error("Resolution check failed");
        }
        auto cppRegistered=sample.Filter("registered_cpp == 1");
        auto cppAngle=cppRegistered.Histo1D({("ca"+std::to_string(m)).c_str(),labels[m],36,0,180},"theta_measured_cpp");
        auto g4Angle=registered.Histo1D({("ga"+std::to_string(m)).c_str(),labels[m],36,0,180},"theta_measured_deg");
        auto cppEnergy=cppRegistered.Histo1D({("ce"+std::to_string(m)).c_str(),labels[m],30,47,53},"energy_measured_cpp");
        auto g4Energy=registered.Histo1D({("ge"+std::to_string(m)).c_str(),labels[m],30,47,53},"E_measured_MeV");
        double countDifference=static_cast<double>(count)-*cppRegistered.Count();
        double countPull=countDifference/std::sqrt(2*(*variance));
        std::cout<<labels[m]<<": direct detector count pull="<<countPull<<'\n';
        if(std::abs(countPull)>6) throw std::runtime_error("Independent detector counts disagree");
        for(int axis=0;axis<2;++axis) {
            auto a=axis==0?cppAngle:cppEnergy;
            auto b=axis==0?g4Angle:g4Energy;
            direct.cd(2*m+axis+1);
            a->SetLineColor(kRed); b->SetLineColor(kBlue);
            a->GetXaxis()->SetTitle(axis==0?"Measured angle (deg)":"Measured energy (MeV)");
            a->GetYaxis()->SetTitle("Registered events / bin");
            a->SetMaximum(1.3*std::max(a->GetMaximum(),b->GetMaximum()));
            a->DrawCopy("E"); b->DrawCopy("E same");
            auto* comparisonLegend=new TLegend(.52,.73,.88,.88);
            comparisonLegend->AddEntry(a.GetPtr(),"C++ RK4","l");
            comparisonLegend->AddEntry(b.GetPtr(),"Geant4 RKF45","l"); comparisonLegend->Draw();
        }
        keep.insert(keep.end(),{cppAngle,g4Angle,cppEnergy,g4Energy});
        keep.insert(keep.end(),{generated,reference,detected,expected});
    }
    auto da=scan.Histo1D({"da","Scan: angle difference;Geant4 - C++ (deg);Trajectories",100,-0.05,0.05},"angle_difference");
    auto ds=scan.Histo1D({"ds","Scan: survival difference;Geant4 - C++;Trajectories",100,-.001,.001},"survival_difference");
    auto de=scan.Histo1D({"de","Geant4 energy conservation;Relative energy error;Trajectories",100,0,1e-5},"energy_relative_error");
    auto dm=events.Histo1D({"dm","MC: angle difference;Geant4 - C++ (deg);Events",100,-.05,.05},"angle_difference");
    checks.cd(1); da->DrawCopy(); checks.cd(2); ds->DrawCopy(); checks.cd(3); de->DrawCopy(); checks.cd(4); dm->DrawCopy();
    direct.SaveAs("results/geant4_detector_direct.pdf");
    angles.SaveAs("results/geant4_model_angles.pdf");
    checks.SaveAs("results/geant4_model_errors.pdf");
    response.SaveAs("results/geant4_model_detector.pdf");
    double maxAngle=std::max(*scan.Max<double>("abs_angle_difference"),*events.Max<double>("abs_angle_difference"));
    double maxSurvival=std::max(*scan.Max<double>("abs_survival_difference"),*events.Max<double>("abs_survival_difference"));
    double maxEnergy=*scan.Max<double>("energy_relative_error");
    std::cout<<"Max angle difference "<<maxAngle<<" deg; max survival difference "<<maxSurvival
             <<"; max Geant4 relative energy error "<<maxEnergy<<'\n';
    if(maxAngle>.05 || maxSurvival>.001 || maxEnergy>1e-5)
        throw std::runtime_error("Numerical model agreement failed");
    std::cout<<"PASS: matched models, detector bookkeeping, efficiency and resolution\n";
}
