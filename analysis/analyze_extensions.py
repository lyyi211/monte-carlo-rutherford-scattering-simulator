"""ROOT plots and conditional detector-bin checks for the added comparisons."""
from pathlib import Path
import csv
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
project = Path(__file__).resolve().parent.parent
results = project / "results"
results.mkdir(exist_ok=True)
keep = []

def graph(x, y, color, style=1):
    from array import array
    obj = ROOT.TGraph(len(x), array("d", x), array("d", y))
    obj.SetLineColor(color)
    obj.SetLineStyle(style)
    obj.SetLineWidth(2)
    keep.append(obj)
    return obj

traces = {}
with (project / "data/trajectory_comparison.csv").open() as stream:
    for row in csv.DictReader(stream):
        traces.setdefault(int(row["case_id"]), []).append({k: float(v) for k, v in row.items()})
assert len(traces) == 24, "Incomplete trajectory matrix"
canvas = ROOT.TCanvas("paths", "Full trajectory comparison", 1200, 1000)
canvas.Divide(2, 2)
labels = ["Point Coulomb", "Uniform sphere", "Fermi Coulomb", "Fermi + Woods-Saxon"]
for model in range(4):
    canvas.cd(model+1)
    frame = ROOT.gPad.DrawFrame(-105, -105, 105, 105, labels[model]+";x (fm);y (fm)")
    legend = ROOT.TLegend(.12, .72, .55, .89)
    for case in range(6):
        rows = traces[6*model+case]
        a = graph([r["x_cpp"] for r in rows], [r["y_cpp"] for r in rows], case+1)
        b = graph([r["x_g4"] for r in rows], [r["y_g4"] for r in rows], case+1, 2)
        a.Draw("L"); b.Draw("L")
        if case == 0:
            legend.AddEntry(a, "C++ RK4 (solid)", "l")
            legend.AddEntry(b, "Geant4 RKF45 (dashed)", "l")
    legend.Draw(); keep.append(legend)
canvas.SaveAs(str(results / "trajectory_comparison.pdf"))
max_position = max(r["position_error"] for rows in traces.values() for r in rows)
max_velocity = max(r["velocity_error"] for rows in traces.values() for r in rows)
assert max_position < .02

with (project / "data/nuclear_physics_benchmark.csv").open() as stream:
    rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
assert [r["E_lab_MeV"] for r in rows] == [5, 10, 25, 35, 50, 75, 100]
canvas2 = ROOT.TCanvas("physics", "Independent nuclear models", 1000, 700)
maximum = max(r[k] for r in rows for k in ("semiclassical_mb", "shen_mb", "glauber_gribov_mb"))
ROOT.gPad.DrawFrame(0, 0, 105, 1.2*maximum, "Model comparison, not experimental validation;Alpha lab energy (MeV);Reaction cross section (mb)")
legend = ROOT.TLegend(.13, .67, .66, .89)
for column, label, color in [
    ("semiclassical_mb", "Classical Fermi + WS absorption", 1),
    ("shen_mb", "Geant4 Shen (low-energy extrapolation)", 2),
    ("glauber_gribov_mb", "Geant4 Glauber-Gribov (low-energy model)", 4),
]:
    g = graph([r["E_lab_MeV"] for r in rows], [r[column] for r in rows], color)
    g.Draw("LP"); legend.AddEntry(g, label, "l")
legend.Draw()
canvas2.SaveAs(str(results / "nuclear_physics_benchmark.pdf"))

# Conditional on the shared incident sample, the two detector draws are independent.
# Bin probabilities integrate the Gaussian response, including under/overflow tails.
file = ROOT.TFile.Open(str(project / "geant4/build/geant4_model_check.root"))
tree = file.Get("events")
stats = {}
def cdf(x):
    return .5*(1+math.erf(x/math.sqrt(2)))
for event in tree:
    row_values = {name: float(getattr(event, name)) for name in ['model_id', 'theta_measured_cpp', 'theta_measured_deg', 'energy_measured_cpp', 'E_measured_MeV', 'theta_cpp_deg', 'theta_true_deg', 'E_true', 'phi_true_deg', 'survival_cpp', 'survival_probability', 'registered_cpp', 'registered']}
    possible_hit = abs(row_values['phi_true_deg']) <= 45 and any(20 <= row_values[name] <= 160 for name in ('theta_cpp_deg', 'theta_true_deg'))
    if not possible_hit:
        assert not row_values['registered_cpp'] and not row_values['registered'], 'Registration outside acceptance'
        continue
    model = int(row_values["model_id"])
    for axis, low, high, n in [("angle", 0., 180., 36), ("energy", 47., 53., 30)]:
        key = model, axis
        if key not in stats:
            stats[key] = [[0., 0., 0., 0.] for _ in range(n+2)]
        width = (high-low)/n
        values = (row_values["theta_measured_cpp"], row_values["theta_measured_deg"]) if axis == "angle" else (row_values["energy_measured_cpp"], row_values["E_measured_MeV"])
        means = (row_values["theta_cpp_deg"], row_values["theta_true_deg"]) if axis == "angle" else (row_values["E_true"], row_values["E_true"])
        sigma = 1. if axis == "angle" else .01*row_values["E_true"]
        for j in range(n+2):
            lo = -math.inf if j == 0 else low+(j-1)*width
            hi = math.inf if j == n+1 else low+j*width
            probs = []
            for source in range(2):
                theta = (row_values["theta_cpp_deg"], row_values["theta_true_deg"])[source]
                hit = 20 <= theta <= 160 and abs(row_values["phi_true_deg"]) <= 45
                survival = (row_values["survival_cpp"], row_values["survival_probability"])[source]
                probs.append(hit*.85*survival*(cdf((hi-means[source])/sigma)-cdf((lo-means[source])/sigma)))
            cell = stats[key][j]
            cell[0] += bool(row_values["registered_cpp"]) and lo <= values[0] < hi
            cell[1] += bool(row_values["registered"]) and lo <= values[1] < hi
            cell[2] += probs[1]-probs[0]
            cell[3] += sum(p*(1-p) for p in probs)
with (project / "data/detector_bin_comparison.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["model_id", "axis", "bin_with_flow", "cpp_count", "g4_count", "expected_difference", "variance", "pull"])
    max_pull = 0
    for (model, axis), bins in stats.items():
        for j, (a, b, expected, variance) in enumerate(bins):
            pull = (b-a-expected)/math.sqrt(variance) if variance > 1e-12 else 0.
            if variance >= 10:
                max_pull = max(max_pull, abs(pull))
            writer.writerow([model, axis, j, a, b, expected, variance, pull])
assert max_pull < 6, "Detector bin comparison failed (bins with variance >= 10)"
print(f"PASS: 24 trajectories, max dx={max_position:.6g} fm, max dv={max_velocity:.6g} c")
print(f"PASS: direct angle/energy detector distributions, max conditional bin pull={max_pull:.4g}")
print("Independent nuclear model differences are reported, not forced to pass an agreement threshold.")
