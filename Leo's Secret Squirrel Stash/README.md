# TDA-Informed Orbital Refueling-Depot Placement

This packet is an educational baseline for selecting static orbital refueling-depot locations while balancing transfer cost and service redundancy. It is designed so a cadet can read every step, replace individual models, and verify intermediate results.

> **Scope warning:** The included transfer-cost function is a transparent surrogate, not a flight-quality trajectory model. Its output is useful for learning the pipeline and testing hypotheses, not for operational placement decisions.

## Research question

Given:

- a weighted set of customer orbits;
- a discrete set of candidate depot orbits;
- a fixed number of depots;
- a required redundancy level (q); and
- a maximum acceptable service cost (epsilon_*);

choose depot locations that minimize transfer cost while keeping the service network connected and redundantly covered.

The pipeline separates the problem into three layers:

1. **Astrodynamics:** Generate a customer-to-depot transfer-cost matrix (C=[c_{ij}]).
2. **Facility location:** Select a subset (S) of candidate depot sites.
3. **Topology and resilience:** Sweep the allowable transfer cost (epsilon), construct service graphs, and calculate Betti curves.

Keeping these layers separate is important. A higher-fidelity transfer model can replace the classroom surrogate without rewriting the optimizer or topology analysis.

## Why TDA belongs here

For a selected depot set, create a bipartite graph with customer vertices on one side and depot vertices on the other. Add edge ((i,j)) when

[
c_{ij}leqepsilon.
]

Sweeping (epsilon) produces a graph filtration.

- (eta_0(epsilon)) is the number of connected components. Large values indicate fragmented or isolated service regions.
- (eta_1(epsilon)=|E|-|V|+eta_0) is the number of independent graph cycles. In this bipartite service graph, cycles indicate alternate customer-depot connections and can serve as a redundancy signal.
- Per-customer degree measures explicit (q)-coverage.

This is a deliberately interpretable first use of TDA. It is not yet the relative-homology coverage certificate used in classical sensor-network theory. A later phase can construct Čech/Rips, nerve, or cubical complexes over a validated orbital-service domain.

## Baseline objective

For customer (i), sort its costs to selected depots:

[
c_{i(1)}leq c_{i(2)}leqcdots.
]

The baseline objective is

[
J(S)=
ar c_{(1)}
+w_bar c_{(q)}
+w_f c_{mathrm{single failure}}
+w_u D_q
+w_0[eta_0(epsilon_*)-1]_+
-w_1eta_1(epsilon_*).
]

Here:

- (ar c_{(1)}) is weighted mean primary-service delta-v;
- (ar c_{(q)}) is weighted mean cost to the (q)-th reachable depot;
- (c_{mathrm{single failure}}) is the worst mean cost after losing one selected depot;
- (D_q) is the weighted (q)-coverage deficit;
- the (eta_0) term penalizes fragmentation; and
- the bounded (eta_1) reward encourages alternate service paths.

The code leaves the terms visible rather than hiding them inside a solver. All coefficients are policy choices and must be justified by sensitivity analysis. In a research implementation, normalize the terms or formulate a constrained multiobjective/Pareto problem so incompatible units are not silently mixed.

## Run the example

Requirements:

- Python 3.10 or newer
- no third-party packages for the baseline

From this directory:

```powershell
python pipeline.py
```

To save the result:

```powershell
python pipeline.py --depots 3 --redundancy 2 --coverage-threshold 12 --output results.json
```

Run the tests:

```powershell
python -m unittest discover -s tests -v
```

If the command `python` is not available on Windows, try `py -3` after installing Python.

## Input files

Both CSV files use:

| Column | Meaning |
|---|---|
| `name` | Unique orbit or site label |
| `altitude_km` | Circular-orbit altitude above the reference Earth radius |
| `inclination_deg` | Inclination |
| `raan_deg` | Right ascension of the ascending node |
| `weight` | Relative customer demand; ignored for candidate depots |

The examples are synthetic and are not a recommended architecture.

## What the transfer surrogate does

`transfer_cost_km_s()` adds:

1. the two impulses of a coplanar Hohmann transfer between circular radii; and
2. a plane-change estimate at the larger radius.

It omits phasing, epoch, argument of latitude, rendezvous time, J2-driven RAAN evolution, finite thrust, eclipse constraints, boil-off, tanker capacity, docking compatibility, inventory, launch replenishment, and uncertainty. These omissions are explicit so the cadet knows exactly where fidelity must be added.

## Suggested cadet implementation sequence

### Milestone 1 — Reproduce the baseline

- Draw the selected bipartite service graph at three (epsilon) values.
- Plot (eta_0(epsilon)) and (eta_1(epsilon)).
- Explain why a cycle represents an alternate service path.
- Verify each unit test by hand.

### Milestone 2 — Sensitivity and validation

- Sweep depot count, redundancy (q), threshold (epsilon_*), demand weights, and objective weights.
- Report Pareto fronts for mean delta-v versus resilience.
- Compare the optimum with random placement and delta-v-only placement.
- Perform leave-one-depot-out and leave-one-orbital-plane-out tests.
- Check whether the choice is stable under cost-matrix perturbations.

### Milestone 3 — Improve orbital fidelity

Replace `transfer_cost_km_s()` with one of:

- impulsive Lambert rendezvous with explicit epochs and time windows;
- Edelbaum or another documented low-thrust approximation;
- a validated low-thrust optimal-control solver; or
- imported costs from mission-analysis software.

Do not call an asymmetric, time-dependent transfer-cost matrix a metric. A Vietoris–Rips complex requires a symmetric dissimilarity; any symmetrization must be stated and tested.

### Milestone 4 — Scale the optimizer

The exact combination search is intentionally simple and grows as (inom{n}{p}). For larger candidate sets, formulate a (p)-median or facility-location mixed-integer program. Continuous refinement of chosen orbital slots should come only after the discrete model and cost calculations have been validated.

### Milestone 5 — Full persistent homology

Use GUDHI or a similar validated library to analyze:

- Rips/Čech complexes when the service-space distance assumptions are defensible;
- a nerve complex of depot service regions;
- a cubical complex over a gridded (q)-th-nearest-depot cost field; or
- persistence diagrams under uncertainty and failure scenarios.

Long-lived topological features should be interpreted physically, not merely rewarded because they are persistent.

## Additional considerations beyond delta-v

A practical architecture should eventually include:

- depot deployment and replenishment cost;
- usable propellant delivered, not just vehicle delta-v;
- tank capacity, throughput, boil-off, and dormancy;
- client arrival rate and queueing;
- time of flight and launch-window availability;
- station-keeping and disposal cost;
- single-depot and common-mode failures;
- communications, navigation, custody, and docking constraints;
- conjunction risk and traffic management;
- uncertainty in state, demand, and maneuver execution; and
- lifecycle cost or equivalent mass to LEO.

The most defensible early result is a **Pareto frontier** rather than a single optimum: expected delta-v, worst-case delta-v, coverage deficit, deployment cost, and resilience should remain separately visible.

## Path to a research-grade pipeline

The classroom baseline is only the first rung. The full development plan is in [RESEARCH_ROADMAP.md](RESEARCH_ROADMAP.md). It covers:

- high-fidelity, time-dependent, directed transfer opportunities;
- mixed-integer facility location and continuous orbit refinement;
- depot inventory, replenishment, fleet state, and stochastic demand;
- directed time-expanded networks and rolling-horizon operation;
- (q)-coverage, disjoint temporal paths, and resource redundancy;
- zigzag persistence, vineyards, directed path homology, and multiparameter extensions;
- uncertainty quantification and statistical inference for topological features;
- solver scaling, surrogate validation, and transfer-cost caching;
- software architecture, provenance, verification, and publication-quality benchmarks.

The recommended progression is **static facility location first, dynamic logistics second, temporal TDA third**. Each layer must outperform or add information beyond simpler graph and optimization baselines.

## Sources and conceptual lineage

- Vin de Silva and Robert Ghrist, [Coverage in sensor networks via persistent homology](https://doi.org/10.2140/agt.2007.7.339), *Algebraic & Geometric Topology* 7 (2007), 339–358.
- Yuri Shimane, Nick Gollins, and Koki Ho, [Orbital Facility Location Problem for Satellite Constellation Servicing Depots](https://arxiv.org/abs/2302.12191) (2023).
- Ian Clark, [An Assessment of On-Orbit Cryogenic Refueling: Optimal Depot Orbits, Launch Vehicle Mass Savings, and Deep Space Mission Opportunities](https://ntrs.nasa.gov/citations/20210014171), NASA Technical Reports Server (2021).
- [GUDHI Python documentation](https://gudhi.inria.fr/python/latest/) for production persistent-homology extensions.
- Gunnar Carlsson and Vin de Silva, [Zigzag Persistence](https://doi.org/10.1007/s10208-010-9066-0), *Foundations of Computational Mathematics* 10 (2010), 367–405.
- Samir Chowdhury and Facundo Mémoli, [Persistent Path Homology of Directed Networks](https://arxiv.org/abs/1701.00565).
- Tristan Sarton du Jonchay, Hao Chen, Onalli Gunasekara, and Koki Ho, [Framework for Modeling and Optimization of On-Orbit Servicing Operations Under Demand Uncertainties](https://doi.org/10.2514/1.A34978), *Journal of Spacecraft and Rockets* (2021).
- Riccardo Apa, Jennifer Hudson, and Marcello Romano, [Low-Thrust Propulsion and Drift Orbit Optimization for Multi-Client Servicing Missions in Low Earth Orbit](https://doi.org/10.1007/s40295-026-00580-4), *The Journal of the Astronautical Sciences* 73, 36 (2026).

## Files

- `pipeline.py`: transfer surrogate, exact placement optimizer, redundancy metrics, and Betti curves
- `data/example_customers.csv`: synthetic demand orbits
- `data/example_candidates.csv`: synthetic depot candidates
- `tests/test_pipeline.py`: unit tests for the cost, topology, and optimizer logic
- `RESEARCH_ROADMAP.md`: staged plan for a dynamic, validated, publication-quality system
