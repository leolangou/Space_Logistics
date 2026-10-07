#DESCRIPTION
'''
This optimization was done on TDA distance data.
We minimize the total round-trip Δv summed over all satellites going to their assigned depots.

1. Every satellite is assigned to exactly one depot.
2. Satellites can only be assigned to open depots.
3. Exactly p depots are open, and the model is re-solved for each p in a sweep. (SUM Y_j = p)
4. The total demand of the satellites assigned to a depot cannot exceed that depot's inventory (2,500 kg).
5. A satellite can only be assigned to a depot if its round-trip Δv is within the 0.70 km/s sortie limit (F_ij = 1 if it is, 0 if not). (Not seen by solver)
6. X_ij and Y_j are binary.
'''


# %% Step A — imports
# Run this in VS Code's Interactive Window (Shift+Enter per cell), from the
# project root, with data_pipeline.py / oflp_geo.py / oflp_fa.py /
# oflp_fa_bridge.py all sitting next to this file.

import data_pipeline as dp
from oflp_fa import FAParams
from oflp_fa_bridge import (
    clients_from_pipeline,
    params_from_pipeline_config,
    make_slots_for_clients,
    run_full_sweep,
)

# %% Step B — run the partner pipeline (STEP 1-4: load, filter, propagate, demand)
targets = dp.run_pipeline()
targets.head()

# %% Step C — map to the solver's client schema
clients = clients_from_pipeline(targets)
clients.describe()

# %% Step D — pull mission parameters straight from data_pipeline.py's config
# (station inventory, sortie dv limit, objective mode, etc. — nothing retyped)
params = params_from_pipeline_config()
print(params)

# %% Step E — build candidate station locations
# Tune lon_step_deg / inc_step_deg here — coarser grid = faster, worse locations.
slots = make_slots_for_clients(clients, lon_step_deg=10.0, inc_step_deg=1.0)
print(f"{len(slots)} candidate slots")

# %% Step F — sweep station count and solve each one to exact optimality
summary, solutions = run_full_sweep(clients, params, slots, time_limit_s=120.0)
summary

# %% Step G — inspect a chosen station count
# Read the marginal_savings column in `summary`: pick the p where adding one
# more station stops being worth it for your purposes, then look at it here.
p_choice = int(summary.loc[summary["feasible"], "p_stations"].min())  # smallest feasible, by default
sol = solutions[p_choice]

print(f"p_stations = {p_choice}   objective = {sol.objective:.3f}")
sol.station_summary

# %% Step H — full client-to-station assignment
sol.assignment

# %%
