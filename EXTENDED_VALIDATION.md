# Geant4 and ROOT comparisons

[Physics and results](NUCLEAR_README.md)

## Trajectory integration

The independent calculation uses Geant4's `G4CashKarpRKF45` stepper, a custom
equation of motion and adaptive error control. Its eight-component state uses
time in fm/c. The force equations are integrated separately from Geant4's
material-transport engine.

Both implementations use the same fixed nucleus, four interaction models and
detector settings. The Geant4 calculation obtains the Fermi potential from
Simpson quadrature of enclosed charge and exterior shells on a 0.002 fm grid;
the C++ RK4 calculation obtains it by integrating the tabulated force.

The comparison covers 1,148 scan rows and 60,000 incident events. It also follows
24 complete paths at 50 MeV: all four models at
$b=0,1,5,8.03816104495,8.25,15$ fm. Positions and velocities are compared at
common 0.1875 fm/c intervals. The adaptive solver shortens its steps to reach
those times. Maximum position differences are below the 0.02 fm comparison
tolerance; escape-time coverage differs by at most two samples.

| Quantity | Maximum difference |
|---|---:|
| Position along the 24 selected paths | 0.000731565 fm |
| Velocity along those paths | $4.32458\times10^{-6}c$ |
| Final scattering angle | 0.000564158 degrees |
| Survival probability | $1.17108\times10^{-6}$ |

## Detector statistics

ROOT compares measured angles and energies with identical bins and absolute
counts. Both calculations start from the same incident sample and use independent
random draws for the detector response. For one measured bin, the event
probability is

```math
p_{\rm bin}=A\,P_{\rm survive}\,\varepsilon
\left[\Phi\left(\frac{x_{\rm upper}-x_{\rm true}}{\sigma}\right)
-\Phi\left(\frac{x_{\rm lower}-x_{\rm true}}{\sigma}\right)\right],
\qquad \varepsilon=0.85.
```

Here $A$ is the angular-acceptance indicator, $\Phi$ is the standard-normal
cumulative distribution function, and $x$ is angle or energy. Conditional on the
incident sample, the expected count difference is
$\sum_i(p_{{\rm G4},i}-p_{{\rm C++},i})$, with variance
$\sum_i[p_{{\rm G4},i}(1-p_{{\rm G4},i})+
p_{{\rm C++},i}(1-p_{{\rm C++},i})]$.

The shared incident events make this conditional comparison more appropriate
than treating the two histograms as independent Poisson samples.
Pulls are evaluated for bins with variance at least 10; sparser bins retain
their counts and expectations without a Gaussian significance assignment.
The largest absolute pull is 2.363, below the comparison threshold of 6.
[Bin counts and pulls](data/detector_bin_comparison.csv).

## Reaction cross sections

Geant4's Shen and Glauber-Gribov models provide independent reaction-cross-section
predictions for alpha + Au-197. These parameterizations are evaluated directly,
separately from the classical force and survival calculation. The corresponding
classical quantity is the loss of elastic-channel flux:

```math
\sigma_R(E)=2\pi\int_0^{b_{\max}} b\,[1-P_{\rm survive}(E,b)]\,db.
```

The input is the total laboratory alpha energy. Cross sections are reported
in mb, with $1\ {\rm fm}^2=10\ {\rm mb}$. The fixed-target, nonrelativistic
trajectory model and the native Geant4 models use different kinematic
approximations, so the comparison probes physical model dependence as well as
the scale of absorption.

Midpoint grids with $\Delta b=0.1$ and $0.05$ fm, a halved time step and
$b_{\max}=30$ instead of 20 fm measure the numerical sensitivity. The optical
parameters are fixed throughout the comparison.

| Alpha lab energy (MeV) | Classical absorption (mb) | Shen (mb) | Glauber-Gribov (mb) |
|---|---:|---:|---:|
| 5 | 0 | 0 | 0 |
| 10 | 1.05e-8 | 0 | 0 |
| 25 | 745.18 | 578.53 | 769.25 |
| 35 | 1528.47 | 1427.89 | 1557.20 |
| 50 | 2075.03 | 2026.94 | 2070.88 |
| 75 | 2435.31 | 2434.85 | 2375.67 |
| 100 | 2586.22 | 2593.19 | 2464.81 |

At 50 MeV, the classical cross section is 2.37% above Shen and 0.20% above
Glauber-Gribov. Its 14.93 mb grid sensitivity is about 0.72%.
The Shen difference reaches 28.8% at 25 MeV; the Glauber-Gribov difference
is 4.93% at 100 MeV. These predictions compare physical models without fitting
their parameters to one another.

The 5-100 MeV alpha range corresponds to 1.25-25 MeV per nucleon. Glauber theory
has a high-energy basis. In Geant4 11.4.3, Shen applies a Coulomb-barrier
condition without an active 10 MeV/nucleon cutoff; its energy correction extends
below the 30 MeV/nucleon lower edge of the underlying Kox data. This is a
comparison of reaction loss, with tunneling and elastic interference outside
the classical calculation.

[Reaction data](data/nuclear_physics_benchmark.csv) |
[Detector-bin data](data/detector_bin_comparison.csv) |
[Trajectory data](data/trajectory_comparison.csv)
