# Research-Grade Development Roadmap

This roadmap turns the educational static-placement baseline into a reproducible research pipeline for designing and operating a resilient orbital refueling network.

## 1. Target scientific claim

A defensible final claim should be narrower than “TDA finds the best depot network.” A stronger formulation is:

> A coupled astrodynamics, logistics, optimization, and topological-analysis pipeline can identify depot architectures that trade propellant cost against service availability and remain robust across time, uncertainty, and failures.

TDA should provide multiscale structural information and robustness diagnostics. It should not replace transfer modeling, network flow, explicit redundancy constraints, or uncertainty analysis.

## 2. End-to-end architecture

The research pipeline should have seven independently testable layers:

```text
Scenario and demand data
        |
        v
Orbit propagation and transfer-opportunity generation
        |
        v
Versioned temporal cost database
        |
        v
Depot placement + fleet/inventory optimization
        |
        v
Dynamic service-network construction
        |
        v
TDA, graph resilience, and uncertainty analysis
        |
        v
Validated reports, figures, and reproducible artifacts
```

A clean interface between layers is essential. In particular, the optimizer should consume a cost/opportunity database rather than call a trajectory solver inside every optimization evaluation.

## 3. Replace the scalar cost matrix with a transfer-opportunity tensor

The classroom model uses one symmetric-looking value (c_{ij}). A research model should represent each feasible directed transfer as

[
a=(i,j,t_{mathrm{dep}},t_{mathrm{arr}},m_0,	ext{vehicle}),
]

with a cost vector such as

[
oldsymbol c_a =
left[
Delta v,,
m_{mathrm{prop}},,
t_{mathrm{flight}},,
P_{mathrm{success}},,
E_{mathrm{energy}},,
C_{mathrm{operations}}
ight].
]

The transfer relation is generally:

- directed;
- time dependent;
- vehicle and mass dependent;
- constrained by rendezvous geometry and docking windows; and
- nonmetric: symmetry and the triangle inequality should not be assumed.

Store the transfer opportunities in a versioned columnar format such as Parquet, with units and provenance attached to every field.

### Astrodynamics fidelity ladder

Build and validate one level at a time:

| Level | Transfer model | Intended use |
|---|---|---|
| F0 | Hohmann plus plane-change surrogate | Unit tests and classroom examples |
| F1 | Lambert or documented impulsive rendezvous model with epochs | Static-placement studies with time windows |
| F2 | Perturbed propagation including (J_2), drag where relevant, eclipse, and station keeping | Operationally meaningful LEO/MEO studies |
| F3 | Finite-thrust or low-thrust optimal control, mass depletion, drift-orbit RAAN correction | High-fidelity servicing and refueling |
| F4 | Cross-validation against an independent mission-analysis tool | Publication-quality validation |

At every level, report position, velocity, epoch, frame, time scale, force model, integrator tolerance, and stopping criteria.

## 4. Static placement as a rigorous facility-location problem

Before introducing network dynamics, replace exhaustive enumeration with a documented mathematical program.

Suggested decision variables:

- (y_jin{0,1}): candidate depot (j) is established;
- (x_{ij}in{0,1}): customer (i) is assigned to depot (j);
- (z_{ij}in{0,1}): depot (j) provides backup coverage to customer (i);
- continuous variables for propellant capacity, inventory, or continuous orbit refinement.

Core constraints should include:

[
x_{ij}leq y_j,
qquad
sum_j z_{ij}geq q,
qquad
z_{ij}=0;	ext{if};c_{ij}>epsilon_*.
]

Use a mixed-integer linear program when the transfer costs are precomputed. Use mixed-integer nonlinear or derivative-free optimization only when continuous orbital decisions require it. Compare every advanced solver against an exact small-instance benchmark.

### Do not collapse the science into one undocumented score

Prefer one of:

1. a constrained problem, such as minimizing expected propellant subject to (q)-coverage and reliability requirements;
2. an (epsilon)-constraint method that traces a Pareto frontier; or
3. an explicitly normalized multiobjective formulation.

At minimum, report these outputs separately:

- expected and worst-case propellant consumption;
- primary and (q)-th-best service cost;
- time of flight;
- depot deployment and replenishment mass;
- unserved demand;
- service delay;
- loss-of-one-depot and common-mode-failure performance;
- topology summaries; and
- lifecycle cost or equivalent mass to LEO.

## 5. Dynamic network formulation

Represent operations as a directed time-expanded multigraph

[
mathcal G=(mathcal V,mathcal A).
]

A vertex can encode orbital location, time, vehicle type, and resource state. Arcs represent:

- wait or loiter;
- customer-to-depot and depot-to-customer transfers;
- refueling and docking operations;
- tanker replenishment;
- launch from Earth;
- maintenance;
- disposal; and
- contingency actions.

For depot (j), a basic propellant balance is

[
I_{j,t+1}
=
I_{j,t}
+D_{j,t}
-S_{j,t}
-L_{j,t},
]

where (I) is inventory, (D) is delivered propellant, (S) is serviced demand, and (L) includes boil-off and other losses. Capacity, flow conservation, launch cadence, vehicle availability, and docking throughput become explicit constraints.

### Three meanings of redundancy

Keep these distinct:

1. **Geometric redundancy:** at least (q) depots are within a transfer threshold.
2. **Temporal-path redundancy:** at least (q) feasible, preferably node-disjoint, time-respecting service paths exist.
3. **Resource redundancy:** those paths have sufficient vehicle, docking, and propellant capacity when needed.

A network can satisfy the first definition while failing the other two. The final pipeline should report all three.

### Rolling-horizon operation

A dynamic implementation should repeatedly:

1. assimilate current ephemerides, inventory, demand, and health state;
2. generate or update feasible transfer arcs;
3. optimize over a finite planning horizon;
4. execute only the first decision interval;
5. update uncertainty and repeat.

This supports stochastic service requests and failures without assuming perfect long-term knowledge.

## 6. TDA progression

### T0 — Graph-filtration baseline

Retain the current (eta_0(epsilon)) and graph-cycle (eta_1(epsilon)) curves for interpretability. Add:

- area under the Betti curves;
- threshold at which all required customers enter one service component;
- threshold at which (q)-coverage is first achieved;
- leave-one-depot-out Betti curves; and
- comparisons against degree, articulation points, min cuts, and node-disjoint paths.

Do not reward raw (eta_1) without a cap or normalization; dense graphs can create many cycles that have little operational value.

### T1 — Service-region complexes

Construct topology from physically defined service regions:

- a nerve complex of depot service regions;
- a cubical complex over a gridded (q)-th-nearest transfer-cost field; or
- Čech/Rips complexes only when their distance assumptions are justified.

This is the appropriate point to connect the project more directly to sensor-coverage theory.

### T2 — Dynamic topology

At each time (t_k), construct a service complex (K_k). Edges and simplices will both appear and disappear as geometry, inventory, demand, and failures change. Standard persistence assumes nested inclusions, which generally do not exist across time.

Use:

- **snapshot persistence** for a simple baseline;
- **vineyards** or distances between consecutive persistence diagrams to track gradual evolution; and
- **zigzag persistence** when complexes gain and lose simplices over time.

Candidate dynamic outputs include:

- duration of full (q)-coverage;
- time spent in fragmented states;
- birth and death of isolated service regions;
- persistence of alternate-route structures;
- topology-change rate; and
- early warning indicators before a coverage loss.

### T3 — Directed and asymmetric topology

Transfer networks are directed. Symmetrizing (c_{ij}(t)) can erase meaningful asymmetry. Compare:

- persistent path homology;
- directed flag complexes;
- Dowker complexes; and
- conventional symmetrized Rips constructions as an ablation.

The directed method should earn its complexity by improving prediction or discrimination of operational failures.

### T4 — Multiple filtration parameters

Operational feasibility depends simultaneously on delta-v, time of flight, risk, and available inventory. Possible research extensions include:

- one-parameter slices through a two-parameter filtration;
- rank invariants or fibered barcodes;
- sensitivity maps over ((Delta v,t_{mathrm{flight}}));
- topology conditioned on failure probability or inventory margin.

Multiparameter persistence is substantially harder to interpret and compute. It should follow, not precede, a validated one-parameter analysis.

## 7. Optimization and TDA coupling

Evaluate three coupling strategies:

### Post hoc diagnostic

Optimize the logistics model first, then use TDA to compare candidate architectures. This is the most interpretable and should be the initial scientific baseline.

### Constraint or screening rule

Require topological and graph-resilience conditions, for example:

[
eta_0(epsilon_*)=1,
qquad
deg(i)geq q,
qquad
kappa(i,	ext{depot set})geq q,
]

where (kappa) is an appropriate node- or edge-connectivity measure.

### Topology-aware objective

Add stable, normalized topological summaries only after demonstrating that they predict meaningful outcomes such as service loss, delay, or propellant shortfall. Use ablation studies to show whether the topology term improves decisions beyond ordinary graph metrics.

## 8. Uncertainty and statistical validity

Model at least:

- orbit-determination covariance;
- maneuver-execution error;
- transfer-model error;
- service-demand time and quantity;
- launch delay;
- depot and servicer failure;
- docking or transfer failure;
- boil-off and propellant gauging uncertainty; and
- uncertain operating cost.

Recommended methods include:

- Monte Carlo ensembles;
- scenario-based stochastic programming;
- chance constraints;
- distributionally robust optimization;
- conditional value at risk;
- robust rolling-horizon optimization; and
- global sensitivity analysis.

For TDA, quantify whether apparent features survive uncertainty:

- compute persistence summaries across the ensemble;
- use bottleneck or Wasserstein distances;
- report confidence bands for Betti curves or persistence landscapes;
- bootstrap when its assumptions are appropriate; and
- distinguish physical signal from discretization and sampling artifacts.

## 9. Computational scaling

High-fidelity transfer generation will dominate runtime. Use:

- cached and versioned transfer solutions;
- parallel opportunity generation;
- sparse graph and matrix representations;
- adaptive sampling of departure time and vehicle mass;
- surrogate models with held-out error bounds;
- active learning near optimizer-relevant boundaries;
- decomposition between placement, routing, and inventory;
- warm starts for rolling-horizon solves; and
- reproducible random seeds.

Every surrogate must expose its domain of validity and prediction uncertainty.

## 10. Software architecture

A research implementation should evolve toward:

```text
ssm/
  config/            scenario schemas and unit definitions
  astrodynamics/     propagation and transfer solvers
  opportunities/     temporal cost database and caching
  logistics/         inventory, demand, and vehicle models
  optimization/      static, stochastic, and rolling-horizon solvers
  networks/          static and time-expanded graph builders
  tda/               filtrations, persistence, and topology summaries
  uncertainty/       sampling, covariance, and scenario generation
  validation/        benchmarks and independent cross-checks
  reporting/         figures, tables, and provenance manifests
  cli/               reproducible command-line workflows
tests/
examples/
docs/
```

Engineering requirements:

- typed configuration schemas;
- explicit units;
- structured logs;
- deterministic seeds;
- immutable scenario identifiers;
- data and model versioning;
- unit, property, integration, and regression tests;
- continuous integration;
- environment lock file;
- archived configuration and provenance for every figure;
- no silent fallback from a failed high-fidelity solver to a surrogate.

## 11. Verification and validation ladder

A result should not advance to the next gate until the current one passes.

| Gate | Evidence |
|---|---|
| V0 — Mathematical correctness | Hand-computed toy cases and dimensional checks |
| V1 — Software correctness | Unit, property, regression, and integration tests |
| V2 — Astrodynamics validation | Comparison with analytic cases and an independent propagator |
| V3 — Optimization validation | Exact small-instance optimum, bounds, and convergence evidence |
| V4 — Topological validity | Known complexes, stability tests, and graph-metric ablations |
| V5 — Uncertainty validity | Calibration, coverage checks, and convergence of ensemble statistics |
| V6 — External validity | Realistic scenarios, expert review, and independent reproduction |

## 12. Benchmark experiment suite

Create fixed, versioned benchmarks:

1. identical coplanar circular orbits;
2. pure altitude change;
3. pure plane change;
4. RAAN-separated constellation with and without (J_2) drift;
5. clustered demand where the static optimum is visually obvious;
6. deliberately disconnected demand regions;
7. single-depot failure and common-plane failure;
8. stochastic demand with inventory depletion;
9. time-varying accessibility that creates and removes service components;
10. high-fidelity comparison case from published literature.

For each benchmark, compare:

- delta-v-only placement;
- conventional (p)-median or facility location;
- explicit (q)-coverage;
- graph-resilience constraints;
- TDA diagnostic only; and
- topology-aware optimization.

## 13. Recommended implementation sequence

### Research-grade static release

1. Define schemas, units, frames, epochs, and scenario provenance.
2. Replace F0 costs with at least F1 transfer opportunities.
3. Implement a mixed-integer (p)-median/facility-location baseline.
4. Implement explicit (q)-coverage and (N-1) failure constraints.
5. Add topology as a post hoc diagnostic.
6. Produce Pareto fronts and uncertainty/sensitivity results.
7. Validate on exact small instances and one published orbital case.

### Dynamic release

1. Introduce time-indexed transfer opportunities.
2. Add depot inventory, tanker replenishment, vehicle state, and service demand.
3. Construct the directed time-expanded network.
4. Implement rolling-horizon optimization.
5. Add snapshot topology and conventional temporal-network metrics.
6. Add zigzag or vineyard analysis.
7. Compare directed persistent path homology against symmetrized baselines.
8. Demonstrate whether TDA supplies information beyond min-cut, path diversity, and degree.

## 14. Publication-quality questions

The strongest paper questions are likely to be:

- When does topology-aware placement select architectures different from delta-v-only facility location?
- Which topological summaries predict service loss under depot or orbital-plane failures?
- How stable are those summaries under orbit, demand, and transfer-model uncertainty?
- Does temporal TDA detect impending loss of service earlier than standard graph metrics?
- When does directed topology materially outperform symmetrized analysis?
- How much additional propellant or deployment mass buys a specified improvement in resilience?
- Are the conclusions robust across LEO, MEO, GEO, and cislunar regimes?

Negative results are valuable. If conventional graph connectivity and (q)-coverage explain all useful behavior, that should be reported rather than forcing a TDA contribution.

## References for the research-grade extension

- Gunnar Carlsson and Vin de Silva, [Zigzag Persistence](https://doi.org/10.1007/s10208-010-9066-0), *Foundations of Computational Mathematics* 10 (2010), 367–405.
- Samir Chowdhury and Facundo Mémoli, [Persistent Path Homology of Directed Networks](https://arxiv.org/abs/1701.00565).
- Peter Bubenik, [Statistical Topological Data Analysis using Persistence Landscapes](https://jmlr.org/papers/v16/bubenik15a.html), *Journal of Machine Learning Research* 16 (2015), 77–102.
- Tristan Sarton du Jonchay, Hao Chen, Onalli Gunasekara, and Koki Ho, [Framework for Modeling and Optimization of On-Orbit Servicing Operations Under Demand Uncertainties](https://doi.org/10.2514/1.A34978), *Journal of Spacecraft and Rockets* (2021).
- Riccardo Apa, Jennifer Hudson, and Marcello Romano, [Low-Thrust Propulsion and Drift Orbit Optimization for Multi-Client Servicing Missions in Low Earth Orbit](https://doi.org/10.1007/s40295-026-00580-4), *The Journal of the Astronautical Sciences* 73, 36 (2026).
- Euihyeon Choi and Koki Ho, [Orbital Depot Location Optimization for Satellite Constellation Servicing with Low-Thrust Transfers](https://hdl.handle.net/1853/80712).
