"""Physical and bookkeeping regressions; Python standard library only."""
from pathlib import Path
import csv
import math
import statistics

ROOT = Path(__file__).resolve().parent.parent


def read(name):
    with (ROOT / "data" / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def check(condition, message):
    if not condition:
        raise AssertionError(message)


profile = read("fermi_radial_profile.csv")
scan = read("nuclear_model_scan.csv")
events = read("detector_events.csv")
summary = read("detector_summary.csv")
models = {"point", "uniform_sphere", "fermi_coulomb", "fermi_nuclear"}
check({s["model"] for s in summary} == models, "Missing model")
r = [float(p["r"]) for p in profile]
rho = [float(p["rho_fermi_e_per_fm3"]) for p in profile]
charge = 4 * math.pi * sum((r[i+1]-r[i]) * (rho[i]*r[i]**2 + rho[i+1]*r[i+1]**2)/2
                           for i in range(len(r)-1))
check(abs(charge-79) < 1e-3, f"Incorrect charge: {charge}")
q = [float(p["Q_fraction"]) for p in profile]
check(all(0 <= x <= 1.00000001 for x in q), "Charge bounds")
check(all(a <= b+1e-12 for a, b in zip(q, q[1:])), "Charge monotonicity")
# Independent central potential integral, not a hardcoded implementation value.
u0 = 4*math.pi*227.514312/79 * sum((r[i+1]-r[i]) * (rho[i]*r[i]+rho[i+1]*r[i+1])/2
                                  for i in range(len(r)-1))
check(abs(float(profile[0]["U_fermi"])-u0) < .001, "Central electrostatic potential")
max_rutherford = 0
for row in scan:
    E, b, theta = (float(row[k]) for k in ("E", "b", "theta_deg"))
    check(row["escaped"] == "1" and math.isfinite(theta), "Invalid trajectory")
    check(0 <= theta <= 180.000001, "Invalid polar angle")
    check(float(row["energy_relative_error"]) < 1e-3, "Energy conservation")
    check(float(row["angular_momentum_error"]) < 1e-4, "Angular momentum conservation")
    check(0 <= float(row["survival_probability"]) <= 1, "Survival bounds")
    if row["model"] != "fermi_nuclear":
        check(float(row["survival_probability"]) == 1, "Spurious absorption")
    if row["model"] == "point":
        error = abs(theta-math.degrees(2*math.atan2(227.514312, 2*E*b)))
        max_rutherford = max(max_rutherford, error)
        check(error < 5e-6, "Rutherford asymptotic angle")
low = {(p["model"], p["b"]): float(p["theta_deg"]) for p in scan if float(p["E"]) == 10}
for (_, b), theta in low.items():
    check(abs(theta-low["point", b]) < .001, "Low-energy Coulomb limit")
head = {p["model"]: float(p["theta_deg"]) for p in scan if float(p["E"]) == 50 and float(p["b"]) == 0}
check(head["point"] > 179 and head["fermi_coulomb"] > 179 and head["uniform_sphere"] < 1,
      "Head-on barrier physics")
for row in summary:
    sample = [p for p in events if p["model"] == row["model"]]
    check(len(sample) == int(row["generated"]), "Generated count")
    check(len({p["event_id"] for p in sample}) == len(sample), "Duplicate event IDs")
    geometric = sum(int(p["geometric_hit"]) for p in sample)
    survivors = sum(int(p["geometric_hit"]) and not int(p["absorbed"]) for p in sample)
    detected = [p for p in sample if int(p["registered"])]
    check(geometric == int(row["geometric_hits"]), "Geometric count")
    check(survivors == int(row["survivors_geometric"]), "Survivor count")
    check(len(detected) == int(row["registered"]), "Registration count")
    check(abs(float(row["registration_fraction"])-len(detected)/len(sample)) < 1e-6, "Registration fraction")
    check(int(row["valid_trajectory"]) == len(sample), "Invalid MC trajectories")
    for p in sample:
        theta, phi, survival = (float(p[k]) for k in ("theta_true_deg", "phi_true_deg", "survival_probability"))
        check(math.isfinite(theta) and 0 <= survival <= 1, "Invalid MC state")
        check(int(p["geometric_hit"]) == (20 <= theta <= 160 and abs(phi) <= 45), "Acceptance")
        if int(p["registered"]):
            check(int(p["geometric_hit"]) and not int(p["absorbed"]), "Registered absorbed/missed event")
        else:
            check(float(p["theta_measured_deg"]) == -1 and float(p["E_measured_MeV"]) == -1, "Missing sentinel")
    if survivors:
        check(abs(len(detected)-.85*survivors) < 6*math.sqrt(.85*.15*survivors), "85% efficiency")
    for key, scale in (("theta_measured_deg", 1.), ("E_measured_MeV", .5)):
        residual = [(float(p[key])-float(p["theta_true_deg"] if key.startswith("theta") else p["E_true"]))/scale
                    for p in detected]
        check(abs(statistics.mean(residual)) < .2 and abs(statistics.stdev(residual)-1) < .2, "Gaussian resolution")
    mean_b2 = statistics.mean(float(p["b"])**2/400 for p in sample)
    check(abs(mean_b2-.5) < 6/math.sqrt(12*len(sample)), "Area-uniform impact parameters")
convergence = read("nuclear_timestep_convergence.csv")
for case in {p["case"] for p in convergence}:
    rows = sorted((p for p in convergence if p["case"] == case), key=lambda p: float(p["dt"]))
    ref = float(rows[0]["theta_deg"])
    nominal = next(p for p in rows if float(p["dt"]) == .1875)
    check(abs(float(nominal["theta_deg"])-ref) < .01, "Timestep convergence")
print(f"PASS: charge={charge:.9f} e; max Rutherford angle error={max_rutherford:.3g} deg")
print(f"Max RK4 energy error={max(float(p['energy_relative_error']) for p in scan):.3g}")
print("PASS: forces/barriers, low-energy limit, conservation, absorption, detector counts/resolution, convergence")
