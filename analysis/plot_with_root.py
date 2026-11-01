"""Regenerate existing figures with PyROOT when matplotlib/pandas are absent."""
import csv
import math
from array import array
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
PROJECT = Path(__file__).resolve().parent.parent
RESULTS = PROJECT / "results"
RESULTS.mkdir(exist_ok=True)
keep = []
colors = [1, 4, 3, 2, 6]
models = ["point", "uniform_sphere", "fermi_coulomb", "fermi_nuclear"]
labels = ["Point", "Uniform sphere", "Fermi Coulomb", "Fermi + nuclear"]


def read(name):
    with (PROJECT / "data" / name).open(newline="") as f:
        return list(csv.DictReader(f))


profile = read("fermi_radial_profile.csv")
scan = read("nuclear_model_scan.csv")
events = read("detector_events.csv")
summary = read("detector_summary.csv")


def graph(rows, x, y, label, color):
    points = [(float(p[x]), float(p[y])) for p in rows]
    points = [(a, b) for a, b in points if math.isfinite(a) and math.isfinite(b)]
    g = ROOT.TGraph(len(points), array("d", (p[0] for p in points)), array("d", (p[1] for p in points)))
    g.SetTitle(label)
    g.SetLineColor(color)
    g.SetLineWidth(2)
    keep.append(g)
    return g


def hist(rows, column, name, color, bins=36, low=0, high=180):
    h = ROOT.TH1D(f"h{len(keep)}", name, bins, low, high)
    for p in rows:
        h.Fill(float(p[column]))
    h.SetLineColor(color)
    h.SetLineWidth(2)
    keep.append(h)
    return h


def plot(filename, panels, nx=2):
    ny = (len(panels)+nx-1)//nx
    c = ROOT.TCanvas(filename, filename, nx*650, ny*480)
    c.Divide(nx, ny)
    for i, (title, bounds, objects, log) in enumerate(panels, 1):
        c.cd(i)
        ROOT.gPad.SetGrid()
        ROOT.gPad.DrawFrame(*bounds, title)
        if log:
            ROOT.gPad.SetLogx(); ROOT.gPad.SetLogy()
        legend = ROOT.TLegend(.48,.68,.89,.89)
        for obj in objects:
            obj.Draw("HIST SAME" if obj.InheritsFrom("TH1") else "L SAME")
            legend.AddEntry(obj, obj.GetTitle(), "l")
        legend.Draw()
        keep.append(legend)
    temporary = RESULTS / ("updated_"+filename)
    c.SaveAs(str(temporary))
    if not temporary.exists():
        raise RuntimeError(f"Cannot save {temporary}")
    temporary.replace(RESULTS / filename)
    keep.append(c)


plot("nuclear_charge_density.png", [
    ("Au-197 charge density;r (fm);Density (e/fm^{3})", (0,0,14,.075),
     [graph(profile,"r",col,label,color) for col,label,color in
      [("rho_fermi_e_per_fm3","Fermi",3),("rho_uniform_e_per_fm3","Uniform sphere",4)]], False),
    ("Enclosed charge;r (fm);Q/(79e)", (0,0,14,1.05), [graph(profile,"r","Q_fraction","Fermi",3)], False)])
plot("nuclear_potential_energy.png", [
    ("Real potential energy;r (fm);U (MeV)", (0,-160,25,120),
     [graph(profile,"r",col,label,color) for col,label,color in
      zip(["U_point","U_uniform","U_fermi","U_fermi_nuclear"],labels,colors)], False)], nx=1)
panels=[]
for E in [10,35,50,75]:
    graphs=[graph([p for p in scan if float(p["E"])==E and p["model"]==m],
                   "b","theta_deg",label,color) for m,label,color in zip(models,labels,colors)]
    panels.append((f"E = {E} MeV;b at infinity (fm);Angle (deg)",(0,0,20,185),graphs,False))
plot("nuclear_scattering_angles.png",panels)
plot("nuclear_absorption.png",[
    ("Radial forces;r (fm);F (MeV/fm)",(0,-95,18,7),
     [graph(profile,"r","F_fermi","Coulomb",3),graph(profile,"r","F_nuclear","Nuclear",2)],False),
    ("Semiclassical survival proxy;b (fm);Survival",(0,0,20,1.05),
     [graph([p for p in scan if p["model"]=="fermi_nuclear" and float(p["E"])==E],
             "b","survival_probability",f"{E} MeV",color) for E,color in zip([25,35,50,75,100],colors)],False)])
panels=[]
for m,label in zip(models,labels):
    rows=[p for p in events if p["model"]==m]
    hs=[hist(rows,"theta_true_deg","Generated",1),
        hist([p for p in rows if p["geometric_hit"]=="1"],"theta_true_deg","Geometric",4),
        hist([p for p in rows if p["registered"]=="1"],"theta_measured_deg","Registered, smeared",2)]
    panels.append((label+";Angle (deg);Events",(0,0,180,hs[0].GetMaximum()*1.15),hs,False))
plot("detector_angular_response.png",panels)
res=[{"residual":float(p["theta_measured_deg"])-float(p["theta_true_deg"])} for p in events
     if p["model"]=="fermi_coulomb" and p["registered"]=="1"]
h=hist(res,"residual","Monte Carlo",4,40,-4,4)
h.Scale(1/h.Integral(),"width")
g=graph([{"x":i/100,"y":math.exp(-(i/100)**2/2)/math.sqrt(2*math.pi)} for i in range(-400,401)],
        "x","y","Gaussian sigma = 1 deg",2)
fractions=ROOT.TH1D("fractions","Registration fraction",4,0,4)
for i,(row,label) in enumerate(zip(summary,labels),1):
    fractions.SetBinContent(i,100*float(row["registration_fraction"]))
    fractions.GetXaxis().SetBinLabel(i,label)
keep.append(fractions)
plot("detector_resolution_efficiency.png",[
    ("Angular resolution;Measured - true angle (deg);Density",(-4,0,4,.6),[h,g],False),
    ("Registration fraction;Model index (point, uniform, Fermi, nuclear);Registered (%)",(0,0,4,12),[fractions],False)])
conv=read("nuclear_timestep_convergence.csv")
angle_graphs=[]; energy_graphs=[]
for case,color in zip(dict.fromkeys(p["case"] for p in conv),colors):
    rows=sorted([dict(p) for p in conv if p["case"]==case],key=lambda p:float(p["dt"]))
    for row in rows:
        row["difference"]=max(1e-12,abs(float(row["theta_deg"])-float(rows[0]["theta_deg"])))
        row["energy_relative_error"]=max(1e-14,float(row["energy_relative_error"]))
    angle_graphs.append(graph(rows,"dt","difference",case,color))
    energy_graphs.append(graph(rows,"dt","energy_relative_error",case,color))
plot("nuclear_timestep_convergence.png",[
    ("Step convergence;dt (fm/c);Angle difference (deg)",(.08,1e-12,2,10),angle_graphs,True),
    ("Energy conservation;dt (fm/c);Relative error",(.08,1e-14,2,.001),energy_graphs,True)])
print("Saved seven current-data figures using ROOT")
