"""
FA — capacitated depot location & assignment, matching the partner pipeline's
mission model (data_pipeline.py):

  * No depot establishment / launch cost (unlike oflp_geo's OFLP, which
    carries the paper's EMLEO facility term). This pipeline never defines a
    depot dry mass or launch vehicle, so there is nothing to price there.
  * Capacity is a station's total propellant INVENTORY in kg
    (STATION_INVENTORY_KG), not a launch-mass cap. The constraint is
    sum_i demand_i * X_ij <= STATION_INVENTORY_KG * Y_j.
  * Demand is PER-CLIENT (DEMAND_KG), not a fixed payload per trip.
  * Objective is either round-trip Delta-v ("delta_v") or a rocket-equation
    propellant estimate ("propellant_kg"), selected by OBJECTIVE_MODE.
  * Station count p is not fixed: the pipeline searches for a station count
    (up to MAX_STATIONS_TO_TEST) via a metaheuristic (GSA). This module
    replaces that search with an EXACT MILP solve at each candidate p, since
    an exact solver is already available -- no metaheuristic tuning needed.

This is the "FA" formulation from the thesis-formulation discussion:
capacitated p-median, Delta-v (or propellant) objective, satellites travel to
a fixed-location depot and back, no facility cost.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

from oflp_geo import A_GEO, MU_EARTH, G0, make_candidate_slots  # reuse constants/slots


# ----------------------------------------------------------------------------
# Phasing Delta-v: k-revolution drift-orbit model
# ----------------------------------------------------------------------------

def phasing_delta_v_kms(delta_lambda_rad: np.ndarray, a_km: np.ndarray,
                        n_revs: int) -> np.ndarray:
    """
    Round-trip Delta-v (km/s) for a phasing maneuver that closes an angular
    separation delta_lambda using n_revs revolutions of a drift orbit.

    Standard small-Delta-a phasing-maneuver approximation (Vallado / Curtis):
    a drift orbit differing from the circular target orbit by Delta-a closes
    delta_lambda after n_revs revolutions when

        |Delta-a / a| = delta_lambda / (3 pi n_revs)

    and a single small tangential burn to reach that Delta-a costs
    dv = (v_circ / 2) * (Delta-a / a); two such burns (enter + return) give

        dv_phase = v_circ * delta_lambda / (3 pi n_revs)

    More revolutions (n_revs) trade transfer TIME for lower Delta-v, matching
    the pipeline's stated model. n_revs = 1 is the fastest / most expensive
    option.
    """
    v_circ = np.sqrt(MU_EARTH / a_km)
    dlam = np.abs(delta_lambda_rad)
    return v_circ * dlam / (3.0 * np.pi * max(n_revs, 1))


def round_trip_delta_v_matrix(clients: pd.DataFrame, slots: pd.DataFrame,
                              n_revs: int) -> np.ndarray:
    """
    Round-trip Delta-v (km/s), shape (n_clients, n_slots): altitude change +
    plane change + phasing, each paid out-and-back.
    """
    a1 = slots["a_km"].to_numpy(float)[None, :]     # depot/station
    a2 = clients["a_km"].to_numpy(float)[:, None]   # client

    # altitude change (Hohmann, out and back)
    at = 0.5 * (a1 + a2)
    dv_alt = (np.abs(np.sqrt(MU_EARTH * (2.0 / a1 - 1.0 / at))
                     - np.sqrt(MU_EARTH / a1))
              + np.abs(np.sqrt(MU_EARTH / a2)
                       - np.sqrt(MU_EARTH * (2.0 / a2 - 1.0 / at))))
    dv_alt = np.where(np.abs(a2 - a1) < 1e-9, 0.0, dv_alt)
    dv_alt_roundtrip = 2.0 * dv_alt

    # plane change (combined inclination difference, paid out and back)
    di = np.radians(np.abs(clients["inc_deg"].to_numpy(float)[:, None]
                           - slots["inc_deg"].to_numpy(float)[None, :]))
    v_pc = np.sqrt(MU_EARTH / np.maximum(a1, a2))
    dv_plane_roundtrip = 2.0 * (2.0 * v_pc * np.sin(0.5 * di))

    # phasing (longitude), already a round-trip quantity by construction
    dlon = np.abs(clients["lon_deg"].to_numpy(float)[:, None]
                  - slots["lon_deg"].to_numpy(float)[None, :]) % 360.0
    dlon = np.minimum(dlon, 360.0 - dlon)
    dv_phase_roundtrip = phasing_delta_v_kms(np.radians(dlon), a1, n_revs)

    return dv_alt_roundtrip + dv_plane_roundtrip + dv_phase_roundtrip


# ----------------------------------------------------------------------------
# Propellant-mass objective (rocket equation, variable demand)
# ----------------------------------------------------------------------------

def round_trip_propellant_kg(dv_round_kms: np.ndarray, demand_kg: np.ndarray,
                             ssc_dry_mass_kg: float,
                             ssc_isp_s: float) -> np.ndarray:
    """
    Round-trip SSC propellant mass, generalizing oflp_geo's fixed-payload
    Eq.(22)-(24) bookkeeping to a PER-CLIENT demand_kg (variable payload).

    Splits the round-trip Delta-v evenly between outbound and inbound legs
    (symmetric geometry assumption, as in oflp_geo):
        m_2+  = ssc_dry_mass_kg * exp(dv_in  / (g0 Isp))     leaves client
        m_2-  = m_2+ + demand_i                              arrives at client
        m_1   = m_2- * exp(dv_out / (g0 Isp))                leaves station
        prop  = m_1 - ssc_dry_mass_kg - demand_i
    """
    dv_leg = dv_round_kms / 2.0
    m2p = ssc_dry_mass_kg * np.exp(dv_leg / (G0 * ssc_isp_s))
    m2m = m2p + demand_kg[:, None]
    m1 = m2m * np.exp(dv_leg / (G0 * ssc_isp_s))
    return m1 - ssc_dry_mass_kg - demand_kg[:, None]


# ----------------------------------------------------------------------------
# Cost matrix + feasibility, per the pipeline's OBJECTIVE_MODE
# ----------------------------------------------------------------------------

@dataclass
class FAParams:
    station_inventory_kg: float = 2500.0
    ssc_payload_capacity_kg: float = 200.0
    max_sortie_dv_kms: float | None = 0.70
    phasing_revolutions: int = 1
    objective_mode: str = "delta_v"     # "delta_v" | "propellant_kg"
    ssc_dry_mass_kg: float = 1000.0
    ssc_isp_s: float = 320.0
    max_stations_to_test: int = 10


def build_cost_matrix(clients: pd.DataFrame, slots: pd.DataFrame, p: FAParams):
    """
    Returns (cost, dv_round, feasible):
      cost      -- (n_clients, n_slots) objective coefficients, per p.objective_mode
      dv_round  -- (n_clients, n_slots) round-trip Delta-v, km/s (always computed,
                   used for the feasibility screen regardless of objective_mode)
      feasible  -- boolean mask, True where dv_round <= p.max_sortie_dv_kms
    """
    dv_round = round_trip_delta_v_matrix(clients, slots, p.phasing_revolutions)

    if p.objective_mode == "delta_v":
        cost = dv_round.copy()
    elif p.objective_mode == "propellant_kg":
        demand = clients["demand_kg"].to_numpy(float)
        cost = round_trip_propellant_kg(dv_round, demand, p.ssc_dry_mass_kg,
                                        p.ssc_isp_s)
    else:
        raise ValueError(f"Unknown objective_mode: {p.objective_mode!r}")

    if p.max_sortie_dv_kms is not None:
        feasible = dv_round <= p.max_sortie_dv_kms
    else:
        feasible = np.ones_like(dv_round, dtype=bool)

    return cost, dv_round, feasible


# ----------------------------------------------------------------------------
# Exact MILP: FA at a fixed station count p
# ----------------------------------------------------------------------------

@dataclass
class FASolution:
    status: str
    feasible: bool
    objective: float
    p: int
    open_slots: list
    assignment: pd.DataFrame
    station_summary: pd.DataFrame
    solve_time_s: float = 0.0


def solve_fa(clients: pd.DataFrame, slots: pd.DataFrame, cost: np.ndarray,
            feasible_mask: np.ndarray, demand_kg: np.ndarray,
            station_inventory_kg: float, p_stations: int,
            time_limit_s: float = 60.0, verbose: bool = False) -> FASolution:
    """
    Exact solve of:

        min_{X,Y}   sum_i sum_j  cost_ij X_ij                       [no f_j term]
        s.t.        sum_j X_ij = 1                    for all i
                    X_ij <= Y_j                        for all i, j
                    sum_j Y_j = p_stations
                    sum_i demand_i X_ij <= station_inventory_kg Y_j  for all j
                    X_ij, Y_j in {0,1}
    """
    t0 = time.time()
    n_i, n_j = cost.shape

    reachable_slot = feasible_mask.any(axis=0)
    if not reachable_slot.any():
        raise RuntimeError("No candidate slot reaches any client under the "
                           "sortie Delta-v limit.")
    slot_idx = np.where(reachable_slot)[0]
    slots_k = slots.iloc[slot_idx].reset_index(drop=True)
    cost_k = cost[:, slot_idx]
    feas_k = feasible_mask[:, slot_idx]
    n_jk = len(slot_idx)

    if n_jk < p_stations:
        return FASolution("infeasible: fewer reachable slots than p_stations",
                          False, np.inf, p_stations, [], pd.DataFrame(),
                          pd.DataFrame(), time.time() - t0)

    unreachable = ~feas_k.any(axis=1)
    if unreachable.any():
        bad = clients.loc[unreachable, "norad_cat_id"].tolist()
        return FASolution(f"infeasible: unreachable clients {bad}", False,
                          np.inf, p_stations, [], pd.DataFrame(),
                          pd.DataFrame(), time.time() - t0)

    nX = n_i * n_jk
    nvars = nX + n_jk

    ii, jj = np.nonzero(feas_k)
    xcols = ii * n_jk + jj
    n_pairs = xcols.size

    obj = np.zeros(nvars)
    obj[xcols] = cost_k[ii, jj]

    rows_l, cols_l, vals_l, lb_l, ub_l = [], [], [], [], []
    r = 0

    # each client assigned exactly once
    rows_l.append(ii); cols_l.append(xcols); vals_l.append(np.ones(n_pairs))
    lb_l.append(np.ones(n_i)); ub_l.append(np.ones(n_i))
    r += n_i

    # X_ij <= Y_j
    pr = r + np.arange(n_pairs)
    rows_l.append(pr); cols_l.append(xcols); vals_l.append(np.ones(n_pairs))
    rows_l.append(pr); cols_l.append(nX + jj); vals_l.append(-np.ones(n_pairs))
    lb_l.append(np.full(n_pairs, -np.inf)); ub_l.append(np.zeros(n_pairs))
    r += n_pairs

    # sum_j Y_j = p_stations
    rows_l.append(np.full(n_jk, r)); cols_l.append(nX + np.arange(n_jk))
    vals_l.append(np.ones(n_jk))
    lb_l.append(np.array([float(p_stations)]))
    ub_l.append(np.array([float(p_stations)]))
    r += 1

    # sum_i demand_i X_ij <= station_inventory_kg * Y_j
    cap_rows = r + jj
    rows_l.append(cap_rows); cols_l.append(xcols)
    vals_l.append(demand_kg[ii])
    rows_l.append(r + np.arange(n_jk)); cols_l.append(nX + np.arange(n_jk))
    vals_l.append(-station_inventory_kg * np.ones(n_jk))
    lb_l.append(np.full(n_jk, -np.inf)); ub_l.append(np.zeros(n_jk))
    r += n_jk

    A = coo_matrix((np.concatenate(vals_l),
                    (np.concatenate(rows_l), np.concatenate(cols_l))),
                   shape=(r, nvars)).tocsc()
    cons = LinearConstraint(A, np.concatenate(lb_l), np.concatenate(ub_l))

    up = np.ones(nvars)
    ni, nj = np.nonzero(~feas_k)
    up[ni * n_jk + nj] = 0.0

    res = milp(c=obj, constraints=cons, integrality=np.ones(nvars),
              bounds=Bounds(np.zeros(nvars), up),
              options={"time_limit": time_limit_s, "disp": verbose})

    if res.x is None or not res.success:
        return FASolution(res.message, False, np.inf, p_stations, [],
                          pd.DataFrame(), pd.DataFrame(), time.time() - t0)

    X = res.x[:nX].reshape(n_i, n_jk)
    Y = res.x[nX:]
    open_j = np.where(Y > 0.5)[0]

    rows = []
    for i in range(n_i):
        j = int(np.argmax(X[i]))
        rows.append({
            "norad_cat_id": clients.iloc[i].get("norad_cat_id"),
            "object_name": clients.iloc[i].get("object_name", ""),
            "demand_kg": float(demand_kg[i]),
            "station_slot": int(slots_k.iloc[j]["slot_id"]),
            "station_lon_deg": float(slots_k.iloc[j]["lon_deg"]),
            "station_inc_deg": float(slots_k.iloc[j]["inc_deg"]),
            "cost": float(cost_k[i, j]),
        })
    assignment = pd.DataFrame(rows)

    summ = []
    for j in open_j:
        members = assignment[assignment["station_slot"]
                             == int(slots_k.iloc[j]["slot_id"])]
        summ.append({
            "station_slot": int(slots_k.iloc[j]["slot_id"]),
            "lon_deg": float(slots_k.iloc[j]["lon_deg"]),
            "inc_deg": float(slots_k.iloc[j]["inc_deg"]),
            "n_clients": int(len(members)),
            "demand_used_kg": float(members["demand_kg"].sum()),
            "inventory_kg": station_inventory_kg,
            "utilization": float(members["demand_kg"].sum()) / station_inventory_kg,
        })
    station_summary = (pd.DataFrame(summ).sort_values("lon_deg")
                       .reset_index(drop=True))

    return FASolution(res.message, True, float(res.fun), p_stations,
                      [int(slots_k.iloc[j]["slot_id"]) for j in open_j],
                      assignment, station_summary, time.time() - t0)


# ----------------------------------------------------------------------------
# Sweep over station count p (replaces the GSA search with exact solves)
# ----------------------------------------------------------------------------

def min_feasible_p(demand_kg: np.ndarray, station_inventory_kg: float) -> int:
    """Hard lower bound on station count from total demand vs. inventory alone
    (ignores reachability -- the true minimum feasible p can be higher)."""
    return int(np.ceil(demand_kg.sum() / station_inventory_kg))


def sweep_stations(clients: pd.DataFrame, slots: pd.DataFrame, p: FAParams,
                   p_min: int | None = None, time_limit_s: float = 60.0,
                   verbose: bool = False) -> tuple[pd.DataFrame, dict[int, FASolution]]:
    """
    Exact solve at every station count from p_min (or the inventory-derived
    lower bound) to p.max_stations_to_test. Returns a summary table and the
    dict of full solutions, so you can pick the count you want by inspecting
    the marginal-savings curve rather than trusting a single heuristic answer.
    """
    demand_kg = clients["demand_kg"].to_numpy(float)
    cost, dv_round, feasible = build_cost_matrix(clients, slots, p)

    lo = p_min if p_min is not None else min_feasible_p(
        demand_kg, p.station_inventory_kg)
    lo = max(lo, 1)

    rows = []
    solutions: dict[int, FASolution] = {}
    prev_obj = None
    for p_stations in range(lo, p.max_stations_to_test + 1):
        sol = solve_fa(clients, slots, cost, feasible, demand_kg,
                       p.station_inventory_kg, p_stations,
                       time_limit_s=time_limit_s, verbose=verbose)
        solutions[p_stations] = sol
        marginal = (prev_obj - sol.objective) if (sol.feasible and prev_obj
                                                  is not None) else np.nan
        rows.append({
            "p_stations": p_stations,
            "feasible": sol.feasible,
            "objective": sol.objective if sol.feasible else np.nan,
            "marginal_savings": marginal,
            "solve_time_s": round(sol.solve_time_s, 2),
        })
        if sol.feasible:
            prev_obj = sol.objective

    summary = pd.DataFrame(rows)
    return summary, solutions
