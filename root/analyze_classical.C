#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <TSystem.h>

// Run from the project root:
// root -l -q root/analyze_classical.C
void analyze_classical() {
    gSystem->mkdir("results",true);
    auto data=ROOT::RDF::FromCSV("data/detector_events.csv");
    auto sample=data.Filter("model == \"fermi_coulomb\"");
    auto accepted=sample.Filter("geometric_hit == 1");
    auto detected=sample.Filter("registered == 1");

    auto hGenerated=sample.Histo1D({"hGenerated","50 MeV classical detector;Angle (deg);Events",90,0,180},
                                  "theta_true_deg");
    auto hAccepted=accepted.Histo1D({"hAccepted","Geometric acceptance",90,0,180},"theta_true_deg");
    auto hMeasured=detected.Histo1D({"hMeasured","Registered, smeared",90,0,180},"theta_measured_deg");

    TCanvas canvas("c1","Classical detector response",1000,650);
    hGenerated->SetLineColor(kBlack);
    hAccepted->SetLineColor(kBlue);
    hMeasured->SetLineColor(kRed);
    hGenerated->SetStats(false);
    hGenerated->Draw("hist");
    hAccepted->Draw("hist same");
    hMeasured->Draw("hist same");
    canvas.BuildLegend();
    canvas.SaveAs("results/root_classical_detector.png");
    std::cout << "Registered events: " << *detected.Count() << std::endl;
}
