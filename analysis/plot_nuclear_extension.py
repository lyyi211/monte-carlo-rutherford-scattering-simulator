from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    import runpy
    import sys
    runpy.run_path(str(Path(__file__).with_name("plot_with_root.py")), run_name="__main__")
    sys.exit(0)

PROJECT = Path(__file__).resolve().parent.parent
DATA = PROJECT / "data"
RESULTS = PROJECT / "results"
RESULTS.mkdir(exist_ok=True)

profile = pd.read_csv(DATA / "fermi_radial_profile.csv")
scan = pd.read_csv(DATA / "nuclear_model_scan.csv")
events = pd.read_csv(DATA / "detector_events.csv")
summary = pd.read_csv(DATA / "detector_summary.csv")

model_labels = {
    "point": "Point nucleus",
    "uniform_sphere": "Uniform sphere",
    "fermi_coulomb": "Fermi charge (Coulomb)",
    "fermi_nuclear": "Fermi + nuclear potential",
}


def save(fig, name):
    fig.tight_layout()
    fig.savefig(RESULTS / name, dpi=220, bbox_inches="tight")
    plt.close(fig)


# PLOT 1: Electrostatic charge density and cumulative charge.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
r = profile.r
axes[0].plot(r, profile.rho_uniform_e_per_fm3, label="Uniform sphere")
axes[0].plot(r, profile.rho_fermi_e_per_fm3, label="Fermi charge distribution")
axes[0].set(xlabel="Radius r (fm)", ylabel=r"Charge density $\rho/e$ ($\mathrm{fm}^{-3}$)",
            title="Gold-197 nuclear charge density", xlim=(0, 14))
axes[0].legend()
axes[1].plot(r, profile.Q_fraction, label="Enclosed fraction Q(r)/Ze")
axes[1].axvline(6.38, ls="--", lw=1, label="Fermi c = 6.38 fm")
axes[1].set(xlabel="Radius r (fm)", ylabel="Enclosed charge fraction",
            title="Charge accumulation", xlim=(0, 14), ylim=(0, 1.04))
axes[1].legend()
for ax in axes:
    ax.grid(alpha=0.25)
save(fig, "nuclear_charge_density.png")


# PLOT 2: Coulomb and real nuclear potential energies.
fig, ax = plt.subplots(figsize=(10, 6))
for name, col in [("Point", "U_point"), ("Uniform sphere", "U_uniform"),
                  ("Fermi Coulomb", "U_fermi"),
                  ("Fermi + real nuclear", "U_fermi_nuclear")]:
    ax.plot(r, profile[col], label=name)
ax.axhline(0, color="gray", lw=0.9)
ax.axvline(6.38, color="gray", ls=":", lw=1, label="Fermi c")
ax.set(xlabel="Radius r (fm)", ylabel="Potential energy U(r) (MeV)",
       title="Electrostatic and real nuclear potential energy", xlim=(0, 25), ylim=(-160, 120))
ax.grid(alpha=0.25)
ax.legend()
save(fig, "nuclear_potential_energy.png")


# PLOT 3: Scattering angles at different energies.
selected_energies = [10, 35, 50, 75]
fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
for ax, energy in zip(axes.flat, selected_energies):
    for model, label in model_labels.items():
        subset = scan[(scan.E == energy) & (scan.model == model) & (scan.escaped == 1)]
        ax.plot(subset.b, subset.theta_deg, label=label, lw=1.7)
    ax.set(title=f"E = {energy:g} MeV", xlim=(0, 20), ylim=(0, 185),
           xlabel="Impact parameter b (fm)", ylabel="Scattering angle (degrees)")
    ax.grid(alpha=0.25)
axes[0, 0].legend(fontsize=8)
fig.suptitle("Classical alpha-gold scattering: four interaction models", fontsize=14)
save(fig, "nuclear_scattering_angles.png")


# PLOT 4: Nuclear-force contribution and semi-classical survival.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].plot(r, profile.F_fermi, label="Coulomb force (Fermi)")
axes[0].plot(r, profile.F_nuclear, label="Attractive nuclear force")
axes[0].axhline(0, color="gray", lw=0.9)
axes[0].set(xlabel="Radius r (fm)", ylabel="Radial force (MeV/fm)",
            title="Opposing radial-force components", xlim=(0, 18))
axes[0].legend()
for energy in [25, 35, 50, 75, 100]:
    s = scan[(scan.model == "fermi_nuclear") & (scan.E == energy)]
    axes[1].plot(s.b, s.survival_probability, label=f"{energy:g} MeV")
axes[1].set(xlabel="Impact parameter b (fm)", ylabel="Survival proxy",
            title="Semiclassical absorption (not a quantum S-matrix)",
            ylim=(-0.04, 1.04), xlim=(0, 20))
axes[1].legend(title="Alpha energy")
for ax in axes:
    ax.grid(alpha=0.25)
save(fig, "nuclear_absorption.png")


# PLOT 5: Generated, geometric, and registered angular distributions.
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
bins = np.arange(0, 181, 5)
for ax, (model, label) in zip(axes.flat, model_labels.items()):
    s = events[events.model == model]
    ax.hist(s.theta_true_deg, bins=bins, histtype="step", label="Generated", lw=1.6)
    ax.hist(s.loc[s.geometric_hit == 1, "theta_true_deg"], bins=bins,
            histtype="step", label="Within angular acceptance", lw=1.6)
    ax.hist(s.loc[s.registered == 1, "theta_measured_deg"].dropna(), bins=bins,
            histtype="step", label="Registered (smeared)", lw=1.6)
    ax.set(xlabel="Scattering angle (degrees)", ylabel="Event count", title=label)
    ax.grid(alpha=0.25)
axes[0, 0].legend(fontsize=8)
fig.suptitle("50 MeV: detector acceptance, efficiency and angular resolution", fontsize=14)
save(fig, "detector_angular_response.png")


# PLOT 6: Angular resolution and registration fraction.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
s = events[(events.model == "fermi_coulomb") & (events.registered == 1)]
residual = s.theta_measured_deg - s.theta_true_deg
axes[0].hist(residual, bins=np.linspace(-4, 4, 41), density=True, alpha=0.65,
             label="Monte Carlo residuals")
xs = np.linspace(-4, 4, 200)
axes[0].plot(xs, np.exp(-xs**2 / 2) / np.sqrt(2*np.pi), label=r"Gaussian $\sigma=1^\circ$")
axes[0].set(xlabel=r"$\theta_{measured}-\theta_{true}$ (degrees)",
            ylabel="Probability density", title="Angular resolution")
axes[0].legend()
axes[1].bar([model_labels[m] for m in summary.model],
            summary.registration_fraction * 100)
axes[1].set(ylabel="Registered / generated (%)", title="Overall registration fraction")
axes[1].tick_params(axis="x", rotation=22)
for ax in axes:
    ax.grid(axis="y", alpha=0.25)
save(fig, "detector_resolution_efficiency.png")


# Numerical comparison table for the documentation.
print("\nDetector summary (E=50 MeV):")
print(summary.to_string(index=False))
print("\nSaved six figures to", RESULTS)

# PLOT 7: Fixed-initial-radius RK4 step-size sensitivity.
conv_path = DATA / "nuclear_timestep_convergence.csv"
if conv_path.exists():
    conv = pd.read_csv(conv_path)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for name, subset in conv.groupby("case"):
        subset = subset.sort_values("dt")
        reference = subset.iloc[0].theta_deg
        delta = np.abs(subset.theta_deg.to_numpy() - reference)
        axes[0].loglog(subset.dt, np.maximum(delta, 1e-12), marker="o", label=name)
        axes[1].loglog(subset.dt, np.maximum(subset.energy_relative_error, 1e-14), marker="o", label=name)
    axes[0].set(xlabel="RK4 step dt (fm/c)", ylabel="Absolute angle difference (degrees)",
                title="Difference relative to smallest available dt")
    axes[1].set(xlabel="RK4 step dt (fm/c)", ylabel="Relative energy-conservation error",
                title="Numerical energy consistency")
    for ax in axes:
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    save(fig, "nuclear_timestep_convergence.png")
    print("Saved timestep convergence figure")
