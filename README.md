# Diabetic Kidney Latent Reserve Model

## Dynamic Risk Surface Model Integrating Albuminuria, eGFR, Cystatin C, and Frailty to Infer Renal Reserve Potential in Diabetes

**Complete theoretical manuscript, reproducible simulation, sensitivity analysis, and validation roadmap**

[![Research Status](https://img.shields.io/badge/status-theoretical%20proof--of--concept-blue)]()
[![Python](https://img.shields.io/badge/Python-reproducible%20simulation-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

Albuminuria and estimated glomerular filtration rate (eGFR) are established axes for evaluating diabetic kidney disease, while cystatin C and frailty may provide additional information on filtration uncertainty and physiological vulnerability.

This repository contains a theoretical mathematical-biology and computational-medicine framework developed to examine whether these interacting signals can be represented as a coupled dynamical system with a latent renal reserve potential.

The model combines:

- six coupled dimensionless state variables;
- bounded nonlinear ordinary differential equations;
- a burden-adjusted renal reserve potential;
- a nonlinear dynamic risk surface;
- a logistic transition hazard;
- three provisional latent reserve regimes;
- an optimal-stopping interpretation;
- equilibrium and Jacobian stability analysis;
- numerical continuation;
- two-parameter risk-surface analysis;
- intervention scenarios;
- Latin hypercube uncertainty analysis;
- partial rank correlation coefficients (PRCC);
- local elasticity analysis;
- solver-tolerance robustness analysis; and
- finite-horizon basin analysis.

The objective is not to provide a clinical score, but to establish a transparent and mathematically testable framework for studying renal vulnerability, recoverability, nonlinear interaction, and intervention timing in diabetic kidney disease.

---

## Research Status and Integrity Statement

> **Working theoretical model — not patient-validated.**

This study is a **literature-parameterized theoretical proof-of-concept**.

Albuminuria, eGFR, cystatin C, and frailty are clinically established observable constructs. However, the renal reserve potential represented in this model is a **latent dimensionless model quantity** and is not equivalent to a directly measured stimulated renal functional reserve test.

No participants were recruited, no biological samples were collected, and no patient-level outcomes were analysed.

Numerical thresholds, transition hazards, reserve regimes, intervention effects, and decision boundaries reported in this repository are outputs of the stated mathematical equations and provisional parameterization.

They are **not**:

- patient data;
- clinical-effect estimates;
- diagnostic cut-offs;
- validated prognostic thresholds; or
- treatment recommendations.

Parameter refinement, retrospective calibration, and external validation are required before clinical interpretation or application.

---

# 1. Clinical and Mathematical Motivation

Diabetic kidney disease is commonly monitored using urinary albumin-to-creatinine ratio and eGFR. These measurements are clinically important but may not evolve in parallel.

Albuminuric injury may occur before substantial filtration loss, while creatinine-based filtration estimates may also be affected by differences in muscle mass. Cystatin C provides an additional filtration marker, and frailty represents physiological vulnerability that may influence treatment tolerance and interpretation of laboratory findings.

The problem considered here is therefore dynamic rather than purely associative.

Two individuals—or the same individual at different times—may have similar eGFR values while occupying different modelled physiological states because glycaemic burden, albuminuric injury, cystatin-C discordance, frailty, metabolic stress, and remaining reserve have evolved differently.

The framework uses concepts inspired by dynamical systems and financial engineering:

- a **dynamic risk surface** represents interacting burdens;
- **latent regimes** distinguish preserved, constrained, and depleted states;
- **uncertainty propagation** evaluates parameter risk; and
- an **optimal-stopping layer** formalises the mathematical comparison between continued observation and intervention under an assumed cost ordering.

The financial-engineering contribution concerns the architecture of the decision problem rather than treating medicine as finance.

---

# 2. Model State Variables

The model contains six dimensionless states.

| Symbol | State | Interpretation |
|---|---|---|
| `G(t)` | Glycaemic burden | Sustained glycaemic stress informed conceptually by HbA1c or glucose exposure |
| `A(t)` | Albuminuric injury | Albuminuria burden informed conceptually by UACR category or longitudinal UACR |
| `C(t)` | Cystatin-C discordance | Burden associated with cystatin-C elevation or discordance between creatinine- and cystatin-C-based eGFR |
| `F(t)` | Frailty burden | Physiological vulnerability informed conceptually by a Frailty Index, Clinical Frailty Scale, or compatible measure |
| `M(t)` | Metabolic-inflammatory stress | Latent shared burden connecting glycaemia and frailty with renal injury |
| `R(t)` | Renal reserve potential | Latent dimensionless capacity to withstand or recover from cumulative stress |

`R(t)` should not be interpreted as a directly measured renal functional reserve test.

---

# 3. Governing Dynamical System

The six-state system is:

```text
dG/dt = u_g(1-G) - k_g G

dA/dt = (a_g G + a_m M + a_f F)(1-A)
        - k_a A - e_rA R A

dC/dt = (c_a A + c_m M)(1-C)
        - k_c C - e_rC R C

dF/dt = (f_g G + f_c C + f_a A)(1-F)
        - k_f F - e_rF R F

dM/dt = (m_g G + m_f F)(1-M)
        - k_m M - e_rM R M

dR/dt = (rho + u_r)(1-R)
        - (d_g G + d_a A + d_c C + d_f F + d_m M)R
```

The glycaemic equation represents an input-clearance balance.

The remaining burden equations use bounded interactions in which input is moderated by remaining capacity. Renal reserve is restored toward one and depleted according to weighted state-dependent stress.

---

# 4. Renal Reserve Potential

The burden-adjusted renal reserve potential is defined as

```text
Q(t) = R(t)
       - 0.15G(t)
       - 0.25A(t)
       - 0.20C(t)
       - 0.18F(t)
       - 0.12M(t)
```

For computational summarisation, three provisional regimes are used:

| Reserve potential | Computational regime |
|---|---|
| `Q >= 0.45` | Preserved |
| `0.25 <= Q < 0.45` | Constrained |
| `Q < 0.25` | Depleted |

These labels are computational classifications and **not clinical staging thresholds**.

---

# 5. Dynamic Risk Surface

The nonlinear dynamic risk surface is

```text
Psi(t) =
    0.14G
  + 0.24A
  + 0.18C
  + 0.18F
  + 0.12M
  + 0.10AF
  + 0.08CF
  - 0.50R
```

A smooth transition hazard is defined by

```text
H(t) = {1 + exp[-8(Psi + 0.10)]}^(-1)
```

The interaction terms `AF` and `CF` allow vulnerability to depend on combinations of states rather than only additive contributions.

`H(t)` is a model-derived transition quantity. It is not a clinically validated probability of kidney failure or another patient outcome.

---

# 6. Optimal-Stopping Interpretation

A sequential decision layer compares the discounted cost of continued observation with the assumed cost of intervention.

The provisional stopping region is defined where

```text
H(t) > 0.60
```

or

```text
Q(t) < 0.25
```

This rule **does not prescribe treatment**.

It identifies mathematical states in which continued observation becomes less favourable under the assumed cost ordering.

Future clinical calibration would need to specify the intervention, treatment burden, adverse-event costs, relevant outcomes, and decision horizon.

---

# 7. Mathematical Analysis

The theoretical analysis considers:

- existence and uniqueness;
- positivity;
- boundedness;
- equilibrium structure;
- numerical Jacobian analysis;
- local stability;
- numerical continuation;
- two-parameter risk-surface geometry;
- finite-horizon basin geometry; and
- perturbation-based numerical stability.

Existence and uniqueness follow from continuous differentiability of the right-hand side on the closed unit hypercube.

Positivity and boundedness are supported by the inward-pointing structure of the system at the state boundaries.

A full Lyapunov proof is **not claimed**. Stability is instead evaluated using Jacobian eigenvalues, long-horizon numerical convergence, and perturbation experiments.

---

# 8. Parameter Strategy

Clinical variables and parameter directions were informed by published literature.

The numerical parameter values used in the present implementation are **provisional dimensionless reference values** rather than patient-fitted estimates.

The reference values include:

| Parameter | Role | Reference value |
|---|---|---:|
| `u_g` | Glycaemic input | 0.230 |
| `k_g` | Glycaemic clearance | 0.340 |
| `a_g` | Glycaemia-to-albuminuria coupling | 0.220 |
| `a_m` | Metabolic-to-albuminuria coupling | 0.080 |
| `a_f` | Frailty-to-albuminuria coupling | 0.050 |
| `k_a` | Albuminuria resolution | 0.280 |
| `e_rA` | Reserve protection of albuminuric state | 0.250 |
| `c_a` | Albuminuria-to-cystatin coupling | 0.180 |
| `c_m` | Metabolic-to-cystatin coupling | 0.100 |
| `k_c` | Cystatin-state resolution | 0.320 |
| `e_rC` | Reserve protection of cystatin state | 0.200 |
| `f_g` | Glycaemia-to-frailty coupling | 0.080 |
| `f_c` | Cystatin-to-frailty coupling | 0.120 |
| `f_a` | Albuminuria-to-frailty coupling | 0.090 |
| `k_f` | Frailty recovery | 0.240 |
| `e_rF` | Reserve protection of frailty | 0.160 |
| `m_g` | Glycaemia-to-metabolic stress | 0.100 |
| `m_f` | Frailty-to-metabolic stress | 0.110 |
| `k_m` | Metabolic stress resolution | 0.300 |
| `e_rM` | Reserve protection of metabolic stress | 0.180 |
| `rho` | Intrinsic reserve restoration | 0.200 |
| `d_g` | Glycaemic reserve damage | 0.100 |
| `d_a` | Albuminuric reserve damage | 0.220 |
| `d_c` | Cystatin reserve damage | 0.160 |
| `d_f` | Frailty reserve damage | 0.180 |
| `d_m` | Metabolic reserve damage | 0.120 |
| `u_r` | Intervention-associated reserve restoration | 0.000 |

These values should not be interpreted as empirically estimated biological rate constants.

---

# 9. Computational Analysis

The reference computational implementation uses:

| Setting | Value |
|---|---|
| Software | Python |
| ODE solver | `scipy.integrate.solve_ivp` |
| Integration method | RK45 |
| Simulation horizon | 12 model years |
| Relative tolerance | `1e-9` |
| Absolute tolerance | `1e-11` |
| State order | `[G, A, C, F, M, R]` |
| Initial state | `[0.28, 0.16, 0.10, 0.14, 0.12, 0.74]` |
| Continuation range for `u_g` | `0.08–0.90` |
| Two-parameter surface | `46 × 46` |
| Latin hypercube sets | `2,500` |
| Basin grid | `31 × 31 = 961` initial states |
| Local elasticity perturbation | `±1%` |
| Random seed | `20260727` |

The workflow includes deterministic trajectories, equilibrium computation, continuation, risk-surface analysis, uncertainty propagation, sensitivity analysis, elasticity analysis, robustness checks, intervention scenarios, and basin analysis.

---

# 10. Reference Equilibrium and Stability

Under the provisional reference parameter set, the computed equilibrium is approximately

```text
[G, A, C, F, M, R]
=
[0.404, 0.200, 0.100, 0.158, 0.125, 0.582]
```

All Jacobian eigenvalues have negative real parts.

The dominant real part is approximately

```text
-0.2160
```

and the fastest reported real part is approximately

```text
-0.5700
```

supporting local asymptotic stability under the reference parameterization.

No uncertainty simulation failed among the 2,500 valid Latin hypercube runs.

Tightening solver tolerances from `10^-6 / 10^-8` to `10^-10 / 10^-12` changed year-12 reserve potential by approximately `0.00000%`.

---

# 11. Reference Dynamics

At model year 12, the reference simulation produces approximately:

```text
G = 0.403
A = 0.197
C = 0.096
F = 0.153
M = 0.122
R = 0.594
```

The corresponding burden-adjusted reserve potential is

```text
Q = 0.423
```

placing the reference trajectory in the model-defined **constrained** regime.

The simulated transition hazard is

```text
H = 0.441
```

The difference between `R` and `Q` reflects the reduction of effective reserve by the weighted combination of glycaemic, albuminuric, cystatin-C, frailty, and metabolic burdens.

---

# 12. Continuation and Tipping Boundaries

Numerical continuation in dimensionless glycaemic input produces a progressive reduction in reserve potential.

The provisional depleted-state boundary

```text
Q = 0.25
```

occurs at approximately

```text
u_g = 0.556
```

The model-defined transition-hazard boundary

```text
H = 0.60
```

occurs earlier at approximately

```text
u_g = 0.393
```

Selected continuation values are:

| `u_g` | Reserve potential |
|---:|---:|
| 0.10 | 0.593 |
| 0.30 | 0.367 |
| 0.60 | 0.237 |

These are parameter-dependent model outputs and should **not** be translated directly into HbA1c or other clinical cut-offs.

---

# 13. Intervention Scenarios

The theoretical intervention scenarios produce:

| Scenario | Reserve potential | Hazard | Regime | Reserve gain | Hazard reduction |
|---|---:|---:|---|---:|---:|
| Reference | 0.423 | 0.441 | Constrained | 0.0% | 0.0% |
| Glycaemic input reduced 30% | 0.499 | 0.339 | Preserved | 17.9% | 23.1% |
| Albuminuric propagation reduced 35% | 0.473 | 0.372 | Preserved | 11.9% | 15.6% |
| Frailty coupling reduced 30% | 0.451 | 0.401 | Preserved | 6.6% | 9.0% |
| Reserve restoration added | 0.509 | 0.354 | Preserved | 20.2% | 19.7% |
| Combined multi-domain strategy | 0.644 | 0.203 | Preserved | 52.2% | 53.9% |

Within this theoretical parameterization, the combined strategy produces the largest increase in reserve potential and the largest reduction in transition hazard among the evaluated scenarios.

These percentages represent **simulation effects**, not clinical treatment-effect estimates.

---

# 14. Global Sensitivity and Uncertainty

Global uncertainty analysis uses 2,500 Latin hypercube parameter sets.

Selected PRCC results are:

| Parameter | PRCC with reserve | PRCC with hazard |
|---|---:|---:|
| `rho` | 0.947 | -0.920 |
| `u_g` | -0.925 | 0.938 |
| `k_g` | 0.802 | -0.830 |
| `a_g` | -0.794 | 0.824 |
| `d_a` | -0.613 | 0.526 |
| `k_a` | 0.592 | -0.635 |
| `f_g` | -0.557 | 0.610 |
| `d_f` | -0.516 | 0.437 |
| `k_f` | 0.473 | -0.525 |
| `c_a` | -0.450 | 0.481 |
| `d_c` | -0.323 | 0.261 |
| `k_c` | 0.269 | -0.294 |

The largest absolute PRCCs for reserve potential are observed for:

```text
rho = 0.947
u_g = -0.925
k_g = 0.802
a_g = -0.794
```

The simulated uncertainty distribution has:

```text
Median reserve potential = 0.427
5th–95th percentile      = 0.255–0.581
Depleted-regime frequency = 4.6%
H > 0.60 frequency        = 11.3%
```

These percentages describe uncertainty under the prespecified parameter ranges. They **do not represent population prevalence or clinical event rates**.

---

# 15. Basin Geometry

Finite-horizon basin analysis evaluates

```text
31 × 31 = 961
```

combinations of initial albuminuric injury and frailty burden.

At the evaluated horizon:

```text
Preserved   : 0 / 961   (0.0%)
Constrained : 961 / 961 (100.0%)
Depleted    : 0 / 961   (0.0%)
```

Higher initial albuminuric injury and frailty burden reduce finite-horizon reserve potential even though the trajectories approach the same local equilibrium over sufficiently long horizons.

This provides a distinction between long-run mathematical stability and finite-horizon recoverability.

---

# 16. Figures

The computational workflow produces five principal figures:

1. **Reserve-potential trajectories** under the reference and intervention scenarios.
2. **Numerical continuation** in dimensionless glycaemic input.
3. **Dynamic reserve surface** over glycaemic input and albuminuric propagation.
4. **Global PRCC sensitivity** of renal reserve potential.
5. **Finite-horizon basin geometry** over initial albuminuric injury and frailty burden.

The figures are stored in the `figures/` directory.

---

# 17. Repository Structure

```text
diabetic-kidney-latent-reserve-model/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
│
├── code/
│   └── renal_risk_surface_model.py
│
├── data/
│   └── outputs/
│       ├── scenario_results.csv
│       ├── continuation_results.csv
│       ├── global_sensitivity_prcc.csv
│       ├── local_elasticity.csv
│       ├── basin_results.csv
│       └── renal_risk_surface_results.csv
│
├── figures/
│   ├── figure1_reserve_trajectories.png
│   ├── figure2_continuation.png
│   ├── figure3_risk_surface.png
│   ├── figure4_prcc.png
│   └── figure5_basin.png
│
└── manuscript/
    └── Dynamic-Risk-Surface-Renal-Reserve-Full-Manuscript.pdf
```

---

# 18. Installation

Clone the repository:

```bash
git clone https://github.com/TitinPrihantini4/diabetic-kidney-latent-reserve-model.git
```

Enter the repository:

```bash
cd diabetic-kidney-latent-reserve-model
```

Install the required Python dependencies:

```bash
pip install -r requirements.txt
```

---

# 19. Reproducing the Analysis

Run the main computational script:

```bash
python code/renal_risk_surface_model.py
```

The accompanying Python implementation is intended to reproduce:

- model equations;
- reference trajectories;
- intervention scenarios;
- equilibrium calculations;
- continuation analysis;
- two-parameter risk surface;
- uncertainty analysis;
- PRCC sensitivity analysis;
- local elasticity analysis;
- solver robustness analysis;
- basin analysis;
- figures; and
- CSV computational outputs.

---

# 20. Data and Output Availability

No original patient-level dataset was used in this theoretical study.

The CSV files distributed with this repository are **model-generated computational outputs**, not clinical datasets.

They contain results from:

- reference simulations;
- intervention scenarios;
- continuation analysis;
- uncertainty analysis;
- global sensitivity analysis;
- local elasticity calculations; and
- finite-horizon basin analysis.

Published clinical literature informed variable selection and parameter direction, while numerical parameters were provisionally scaled and subjected to uncertainty and sensitivity analyses.

---

# 21. Limitations

The current implementation has several important limitations:

1. The parameter set was not fitted to a patient cohort.
2. Renal reserve potential is a theoretical latent construct and is not equivalent to stimulated renal functional reserve.
3. Cystatin-C discordance is represented as a single state although assay variation, inflammation, thyroid status, and corticosteroid exposure may affect cystatin C.
4. Frailty is represented as a single burden even though different frailty instruments may not be equivalent.
5. The metabolic-inflammatory state is latent and may not be identifiable without additional biomarkers.
6. The reserve-regime and stopping boundaries are computational labels rather than clinical thresholds.
7. Competing death, acute kidney injury, medication classes, blood pressure, cardiovascular disease, and dialysis are not explicitly modelled.
8. PRCC rankings are conditional on the selected parameter ranges and deterministic model structure.

---

# 22. Validation Roadmap

The model is intended for later clinical calibration rather than immediate clinical use.

A future validation programme should include:

1. protocol freezing and preregistration;
2. systematic parameter-provenance review;
3. selection of a longitudinal calibration cohort;
4. mapping of clinical measurements to dimensionless model states;
5. structural identifiability assessment;
6. parameter estimation using clinical longitudinal data;
7. practical identifiability analysis;
8. validation of latent-state transitions against prespecified kidney outcomes;
9. predictive-performance and calibration assessment;
10. expanded uncertainty analysis;
11. external temporal and geographical validation; and
12. reproducible archival of code, environment, parameter tables, and simulation outputs.

Potential future clinical measurements include longitudinal UACR, creatinine, cystatin C, frailty or functional measures, diabetes duration, medication exposure, and kidney outcomes.

---

# 23. Intended Research Use

The immediate purpose of this framework is not bedside prediction.

Potential research uses include:

- guiding variable selection;
- structuring longitudinal monitoring studies;
- identifying parameters requiring empirical estimation;
- comparing theoretical multi-domain intervention scenarios;
- studying stability and recoverability;
- examining nonlinear vulnerability interactions;
- defining candidate outcomes for retrospective calibration;
- studying intervention timing; and
- supporting collaboration between mathematical modellers, nephrologists, diabetologists, geriatricians, and clinical epidemiologists.

---

# 24. Reproducibility Record

```text
Software:
Python

ODE solver:
scipy.integrate.solve_ivp

Method:
RK45

Simulation horizon:
12 model years

State order:
G, A, C, F, M, R

Reference initial state:
[0.28, 0.16, 0.10, 0.14, 0.12, 0.74]

Relative tolerance:
1e-9

Absolute tolerance:
1e-11

Latin hypercube sets:
2,500

Basin grid:
31 × 31

Two-parameter surface:
46 × 46

Random seed:
20260727
```

---

# 25. Manuscript

The complete accompanying manuscript is:

**Dynamic Risk Surface Model Integrating Albuminuria, eGFR, Cystatin C, and Frailty to Infer Renal Reserve Potential in Diabetes**

The manuscript contains the full theoretical formulation, parameter strategy, computational analysis, results, discussion, limitations, validation roadmap, and references.

---

# 26. Citation

A `CITATION.cff` file will provide structured citation metadata for this repository.

A persistent DOI will be added following archival release.

If you use the mathematical framework, code, computational outputs, or figures in academic work, please cite the corresponding archived release.

---

# 27. License

The computational code in this repository is distributed under the **MIT License**.

See the [`LICENSE`](LICENSE) file for details.

The presence of the manuscript in this repository should not be interpreted as transferring copyright beyond the permissions explicitly stated for the associated research output.

---

# 28. Clinical Disclaimer

**This repository is intended for research and methodological development only.**

The model has not been clinically calibrated or externally validated.

The latent renal reserve potential, computational regimes, transition hazard, numerical boundaries, and simulated intervention effects must not be used for:

- diagnosis;
- individual risk prediction;
- treatment selection;
- medication adjustment;
- clinical staging; or
- direct patient-management decisions.

Clinical calibration and external validation are required before clinical use.

---

# 29. Author

**Prihantini**  
Financial Engineering, WorldQuant University, USA

---

## Research Areas

Mathematical Biology · Computational Medicine · Dynamical Systems · Diabetic Kidney Disease · Sensitivity Analysis · Uncertainty Quantification · Optimal Stopping
