"""
Orbital Facility Location Problem (OFLP) adapted to a GEO client set.

Formulation follows Shimane, Gollins & Ho (2024), "Orbital Facility Location
Problem for Satellite Constellation Servicing Depots", J. Spacecraft & Rockets
61(3), Eqs. (15)-(24).

Adaptation notes (see README section in the accompanying report):
  * The paper targets MEO constellations spread over many orbital planes and
    uses a Q-law low-thrust controller to price each transfer. GEO clients are
    near-coplanar and near-circular, so transfers here are priced with an
    impulsive Delta-v budget (phasing + plane change + altitude change) rather
    than Q-law. The paper itself notes GEO servicing needs only a phasing
    maneuver in the coplanar limit.
  * The EMLEO facility-usage cost (Eqs. 15-19) carries over unchanged.
  * The backward-in-time round-trip propellant accounting (Eqs. 22-24) carries
    over unchanged, with the rocket equation replacing Q-law integration.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from scipy.optimize import milp, LinearConstraint, Bounds

MU_EARTH = 398600.4418          # km^3/s^2
R_EARTH = 6378.137              # km
G0 = 9.80665e-3                 # km/s^2
OMEGA_EARTH = 7.2921150e-5      # rad/s
A_GEO = 42164.17                # km
V_GEO = np.sqrt(MU_EARTH / A_GEO)   # km/s, ~3.0747


# ----------------------------------------------------------------------------
# Parameters (Table 1 of the paper, with GEO-appropriate defaults)
# ----------------------------------------------------------------------------

@dataclass
class OFLPParams:
    # --- launch / depot establishment (paper Table 1) ---
    m_l_max: float = 12950.0        # kg, max launch mass (Ariane 64 sub-GEO)
    r0: float = 6578.0              # km, LV parking orbit radius
    isp_l: float = 457.0            # s, launch vehicle
    isp_d: float = 320.0            # s, depot insertion stage
    m_d_dry: float = 1500.0         # kg, depot dry mass

    # --- servicer ---
    m_s_dry: float = 500.0          # kg, servicer dry mass
    m_s_L: float = 100.0            # kg, payload delivered per trip
    isp_s: float = 1790.0           # s, servicer
    D: int = 1                      # trips per client over depot life

    # --- GEO transfer model ---
    t_transfer_max_days: float = 30.0   # time allowed for a phasing maneuver
    phasing_mode: str = "charged"       # "charged" | "free"
    dv_margin: float = 1.0              # multiplier on computed Delta-v

    # --- feasibility screening ---
    max_alloc_dv_kms: float = 1.5   # drop (i,j) pairs above this round-trip dv


# ----------------------------------------------------------------------------
# Client element extraction
# ----------------------------------------------------------------------------

def elements_from_gp(df: pd.DataFrame) -> pd.DataFrame:
    """Preferred path: pull mean elements straight from GP / TLE records.

    Expects Space-Track GP column names. Returns one row per satellite using
    the most recent epoch available for that satellite.
    """
    cols = {c.upper(): c for c in df.columns}
    need = ["NORAD_CAT_ID", "MEAN_MOTION", "ECCENTRICITY", "INCLINATION",
            "RA_OF_ASC_NODE"]
    missing = [c for c in need if c not in cols]
    if missing:
        raise KeyError(f"GP frame missing columns: {missing}")

    d = df.copy()
    if "EPOCH" in cols:
        d[cols["EPOCH"]] = pd.to_datetime(d[cols["EPOCH"]], utc=True,
                                          errors="coerce")
        d = d.sort_values(cols["EPOCH"]).groupby(cols["NORAD_CAT_ID"]).tail(1)

    n_rev_day = pd.to_numeric(d[cols["MEAN_MOTION"]], errors="coerce")
    n_rad_s = n_rev_day * 2.0 * np.pi / 86400.0
    a_km = (MU_EARTH / n_rad_s**2) ** (1.0 / 3.0)

    out = pd.DataFrame({
        "norad_cat_id": pd.to_numeric(d[cols["NORAD_CAT_ID"]],
                                      errors="coerce").astype("Int64"),
        "a_km": a_km.to_numpy(),
        "ecc": pd.to_numeric(d[cols["ECCENTRICITY"]],
                             errors="coerce").to_numpy(),
        "inc_deg": pd.to_numeric(d[cols["INCLINATION"]],
                                 errors="coerce").to_numpy(),
        "raan_deg": pd.to_numeric(d[cols["RA_OF_ASC_NODE"]],
                                  errors="coerce").to_numpy(),
    })
    if "OBJECT_NAME" in cols:
        out["object_name"] = d[cols["OBJECT_NAME"]].to_numpy()

    # Sub-satellite longitude at epoch is not recoverable from mean elements
    # alone without the epoch GMST; callers who need it should use the TEME
    # path below. Left as NaN here and filled by merge if available.
    out["lon_deg"] = np.nan
    return out.reset_index(drop=True)


def elements_from_teme_timeseries(good: pd.DataFrame) -> pd.DataFrame:
    """Fallback path: derive working elements from the propagated TEME series.

    Uses the fact that a near-circular GEO orbit's out-of-plane excursion over
    a full day has amplitude r*sin(i), so inclination can be recovered from the
    z-amplitude of the timeseries even though a single position cannot give it.
    Longitude is taken at the first timestep.
    """
    d = good.copy()
    d["time_utc"] = pd.to_datetime(d["time_utc"], utc=True)
    r = np.sqrt(d["x_teme_km"]**2 + d["y_teme_km"]**2 + d["z_teme_km"]**2)
    d["radius_km"] = r
    d["z_over_r"] = d["z_teme_km"] / r

    grp = d.groupby("norad_cat_id")
    a_km = grp["radius_km"].mean()
    # sin(i) ~ max |z|/r sampled over the span
    sin_i = grp["z_over_r"].apply(lambda s: np.nanmax(np.abs(s.to_numpy())))
    inc_deg = np.degrees(np.arcsin(np.clip(sin_i, 0.0, 1.0)))
    # crude eccentricity proxy from radius spread
    ecc = (grp["radius_km"].max() - grp["radius_km"].min()) / (2.0 * a_km)

    t0 = d["time_utc"].min()
    snap = d[d["time_utc"] == t0].set_index("norad_cat_id")
    lon = np.degrees(np.arctan2(snap["y_teme_km"], snap["x_teme_km"])) % 360.0

    out = pd.DataFrame({
        "norad_cat_id": a_km.index.to_numpy(),
        "a_km": a_km.to_numpy(),
        "ecc": ecc.reindex(a_km.index).to_numpy(),
        "inc_deg": inc_deg.reindex(a_km.index).to_numpy(),
        "lon_deg": lon.reindex(a_km.index).to_numpy(),
    })
    if "object_name" in d.columns:
        names = grp["object_name"].first()
        out["object_name"] = names.reindex(a_km.index).to_numpy()
    out["raan_deg"] = np.nan
    return out.reset_index(drop=True)


# ----------------------------------------------------------------------------
# Candidate depot slots (paper Sec. III.C — discretization)
# ----------------------------------------------------------------------------

def make_candidate_slots(
    lon_deg=np.arange(0.0, 360.0, 10.0),
    dalt_km=(-300.0, 0.0, 300.0, 800.0),
    inc_deg=(0.0, 1.0, 3.0),
    ecc=(0.0,),
) -> pd.DataFrame:
    """Discretized XGEO/GEO depot slots.

    dalt_km is the altitude offset from the GEO radius; negative is sub-GEO,
    positive is the XGEO / near-graveyard band.
    """
    rows = []
    for lo in lon_deg:
        for da in dalt_km:
            for ic in inc_deg:
                for e in ecc:
                    rows.append({
                        "slot_id": len(rows),
                        "lon_deg": float(lo),
                        "dalt_km": float(da),
                        "a_km": A_GEO + float(da),
                        "inc_deg": float(ic),
                        "ecc": float(e),
                    })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Facility usage cost: EMLEO mass ratios (paper Eqs. 15-19)
# ----------------------------------------------------------------------------

def emleo_mass_ratios(a_km, ecc, p: OFLPParams):
    """Return (phi, phi_d) per slot.

    phi   = phi_d * phi_l : full EMLEO ratio, used in the objective.
    phi_d = depot's own insertion burn ratio, used in the launch-mass constraint.
    Follows Eqs. (15)-(19): a coplanar Hohmann from circular LEO r0 to either
    the slot's perigee or apogee, whichever gives the smaller ratio.
    """
    a = np.asarray(a_km, dtype=float)
    e = np.asarray(ecc, dtype=float)
    r0 = p.r0
    rp = a * (1.0 - e)
    ra = a * (1.0 + e)

    def _pair(r_target):
        # burn 1 (launch vehicle): circular r0 -> transfer ellipse r0 x r_target
        dv1 = (np.sqrt(MU_EARTH * (2.0 / r0 - 2.0 / (r0 + r_target)))
               - np.sqrt(MU_EARTH / r0))
        # burn 2 (depot): transfer ellipse -> target orbit at r_target
        dv2 = (np.sqrt(MU_EARTH * (2.0 / r_target - 1.0 / a))
               - np.sqrt(MU_EARTH * (2.0 / r_target - 2.0 / (r0 + r_target))))
        dv1 = np.abs(dv1)
        dv2 = np.abs(dv2)
        phi_l = np.exp(dv1 / (G0 * p.isp_l))
        phi_d = np.exp(dv2 / (G0 * p.isp_d))
        return phi_d * phi_l, phi_d

    phi_p, phid_p = _pair(rp)
    phi_a, phid_a = _pair(ra)
    use_p = phi_p <= phi_a
    phi = np.where(use_p, phi_p, phi_a)
    phi_d = np.where(use_p, phid_p, phid_a)
    return phi, phi_d


# ----------------------------------------------------------------------------
# Allocation cost: GEO round-trip Delta-v and propellant (paper Eqs. 22-24)
# ----------------------------------------------------------------------------

def geo_transfer_dv(slot, client, p: OFLPParams):
    """One-way Delta-v (km/s) from a depot slot to a client orbit.

    Three additive contributions, all impulsive:
      1. altitude change   : Hohmann between the two semimajor axes
      2. plane change      : combined inclination difference at GEO speed
      3. longitude phasing : drift-orbit phasing within t_transfer_max
    """
    a1 = float(slot["a_km"])
    a2 = float(client["a_km"])

    # 1. Hohmann altitude change
    if abs(a2 - a1) < 1e-9:
        dv_alt = 0.0
    else:
        at = 0.5 * (a1 + a2)
        dv_alt = (abs(np.sqrt(MU_EARTH * (2.0 / a1 - 1.0 / at))
                      - np.sqrt(MU_EARTH / a1))
                  + abs(np.sqrt(MU_EARTH / a2)
                        - np.sqrt(MU_EARTH * (2.0 / a2 - 1.0 / at))))

    # 2. plane change at the higher of the two radii (cheaper there)
    di = np.radians(abs(float(client["inc_deg"]) - float(slot["inc_deg"])))
    v_pc = np.sqrt(MU_EARTH / max(a1, a2))
    dv_plane = 2.0 * v_pc * np.sin(0.5 * di)

    # 3. longitude phasing
    if p.phasing_mode == "free":
        dv_phase = 0.0
    else:
        dlon = abs(float(client["lon_deg"]) - float(slot["lon_deg"])) % 360.0
        dlon = min(dlon, 360.0 - dlon)
        T = p.t_transfer_max_days * 86400.0
        dv_phase = 2.0 * A_GEO * np.radians(dlon) / (3.0 * T)

    return p.dv_margin * (dv_alt + dv_plane + dv_phase)


def transfer_dv_matrix(slots: pd.DataFrame, clients: pd.DataFrame,
                       p: OFLPParams) -> np.ndarray:
    """Vectorized one-way Delta-v, shape (n_clients, n_slots). Same model as
    geo_transfer_dv, which is retained for readability / single-pair use."""
    a1 = slots["a_km"].to_numpy(float)[None, :]      # depot
    a2 = clients["a_km"].to_numpy(float)[:, None]    # client

    at = 0.5 * (a1 + a2)
    dv_alt = (np.abs(np.sqrt(MU_EARTH * (2.0 / a1 - 1.0 / at))
                     - np.sqrt(MU_EARTH / a1))
              + np.abs(np.sqrt(MU_EARTH / a2)
                       - np.sqrt(MU_EARTH * (2.0 / a2 - 1.0 / at))))
    dv_alt = np.where(np.abs(a2 - a1) < 1e-9, 0.0, dv_alt)

    di = np.radians(np.abs(clients["inc_deg"].to_numpy(float)[:, None]
                           - slots["inc_deg"].to_numpy(float)[None, :]))
    v_pc = np.sqrt(MU_EARTH / np.maximum(a1, a2))
    dv_plane = 2.0 * v_pc * np.sin(0.5 * di)

    if p.phasing_mode == "free":
        dv_phase = 0.0
    else:
        dlon = np.abs(clients["lon_deg"].to_numpy(float)[:, None]
                      - slots["lon_deg"].to_numpy(float)[None, :]) % 360.0
        dlon = np.minimum(dlon, 360.0 - dlon)
        T = p.t_transfer_max_days * 86400.0
        dv_phase = 2.0 * A_GEO * np.radians(dlon) / (3.0 * T)

    return p.dv_margin * (dv_alt + dv_plane + dv_phase)


def allocation_cost_matrix(slots: pd.DataFrame, clients: pd.DataFrame,
                           p: OFLPParams):
    """Round-trip servicer propellant mass c_tilde[i, j] in kg.

    Backward-in-time accounting exactly as in Eqs. (22)-(24):
        m_s,2+ = m_s,dry * exp(dv_in  / (g0 Isp))     leaves client
        m_s,2- = m_s,2+ + m_s,L                       arrives at client
        m_s,1  = m_s,2- * exp(dv_out / (g0 Isp))      leaves depot
        c      = m_s,1 - m_s,dry - m_s,L
    Returns (c_tilde, dv_round, feasible_mask).
    """
    dv_out = transfer_dv_matrix(slots, clients, p)
    dv_in = dv_out          # symmetric geometry; mass differs, not dv
    dvr = dv_out + dv_in

    m2p = p.m_s_dry * np.exp(dv_in / (G0 * p.isp_s))
    m2m = m2p + p.m_s_L
    m1 = m2m * np.exp(dv_out / (G0 * p.isp_s))
    c = m1 - p.m_s_dry - p.m_s_L

    feasible = dvr <= p.max_alloc_dv_kms
    return c, dvr, feasible


# ----------------------------------------------------------------------------
# The binary linear program (paper Eqs. 20-21)
# ----------------------------------------------------------------------------

@dataclass
class OFLPSolution:
    status: str
    objective_emleo_kg: float
    open_slots: list
    assignment: pd.DataFrame
    depot_summary: pd.DataFrame
    solve_time_s: float = 0.0
    n_vars: int = 0


def solve_oflp(slots: pd.DataFrame, clients: pd.DataFrame, p: OFLPParams,
               time_limit_s: float = 300.0, verbose: bool = True) -> OFLPSolution:
    """Solve the OFLP.

        min  sum_j m_d,dry phi_j Y_j
           + sum_i sum_j D (c_ij + m_s,L) phi_j X_ij
        s.t. sum_j X_ij = 1                      for all clients i
             X_ij <= Y_j                         for all i, j
             m_d,dry phi_d_j Y_j
               + sum_i D (c_ij + m_s,L) phi_d_j X_ij <= m_l,max   for all j
             X, Y binary
    """
    import time
    t_start = time.time()

    phi, phi_d = emleo_mass_ratios(slots["a_km"].to_numpy(),
                                   slots["ecc"].to_numpy(), p)
    c_til, dv_round, feasible = allocation_cost_matrix(slots, clients, p)

    n_i, n_j = c_til.shape

    # Prune slots that cannot reach at least one client, and prune infeasible
    # (i, j) pairs by freezing X_ij = 0 (paper Sec. III.D).
    reachable_slot = feasible.any(axis=0)
    if not reachable_slot.any():
        raise RuntimeError("No candidate slot can reach any client under "
                           "max_alloc_dv_kms; loosen the screen.")
    slot_idx = np.where(reachable_slot)[0]
    slots_k = slots.iloc[slot_idx].reset_index(drop=True)
    phi_k, phid_k = phi[slot_idx], phi_d[slot_idx]
    c_k = c_til[:, slot_idx]
    feas_k = feasible[:, slot_idx]
    dv_k = dv_round[:, slot_idx]
    n_jk = len(slot_idx)

    unreachable_clients = ~feas_k.any(axis=1)
    if unreachable_clients.any():
        bad = clients.loc[unreachable_clients, "norad_cat_id"].tolist()
        raise RuntimeError(f"Clients unreachable from every slot: {bad}")

    # Variable layout: [X (n_i * n_jk) | Y (n_jk)]
    nX = n_i * n_jk
    nvars = nX + n_jk

    def xi(i, j):
        return i * n_jk + j

    # --- objective ---
    alloc_cost = p.D * (c_k + p.m_s_L) * phi_k[None, :]   # kg EMLEO per (i,j)
    obj = np.concatenate([alloc_cost.ravel(), p.m_d_dry * phi_k])

    # --- constraints (assembled vectorized) ---
    ii, jj = np.nonzero(feas_k)                 # feasible (client, slot) pairs
    xcols = ii * n_jk + jj
    n_pairs = xcols.size

    rows_l, cols_l, vals_l, lb_l, ub_l = [], [], [], [], []
    r = 0

    # (20b) each client assigned exactly once
    rows_l.append(ii)
    cols_l.append(xcols)
    vals_l.append(np.ones(n_pairs))
    lb_l.append(np.ones(n_i)); ub_l.append(np.ones(n_i))
    r += n_i

    # (20c) X_ij - Y_j <= 0, one row per feasible pair
    pair_rows = r + np.arange(n_pairs)
    rows_l.append(pair_rows);   cols_l.append(xcols)
    vals_l.append(np.ones(n_pairs))
    rows_l.append(pair_rows);   cols_l.append(nX + jj)
    vals_l.append(-np.ones(n_pairs))
    lb_l.append(np.full(n_pairs, -np.inf)); ub_l.append(np.zeros(n_pairs))
    r += n_pairs

    # (21b) launch mass cap per depot
    cap_rows = r + np.arange(n_jk)
    rows_l.append(cap_rows); cols_l.append(nX + np.arange(n_jk))
    vals_l.append(p.m_d_dry * phid_k)
    rows_l.append(r + jj);   cols_l.append(xcols)
    vals_l.append(p.D * (c_k[ii, jj] + p.m_s_L) * phid_k[jj])
    lb_l.append(np.full(n_jk, -np.inf)); ub_l.append(np.full(n_jk, p.m_l_max))
    r += n_jk

    from scipy.sparse import coo_matrix
    A = coo_matrix((np.concatenate(vals_l),
                    (np.concatenate(rows_l), np.concatenate(cols_l))),
                   shape=(r, nvars)).tocsc()
    cons = LinearConstraint(A, np.concatenate(lb_l), np.concatenate(ub_l))

    # frozen infeasible pairs
    up = np.ones(nvars)
    infeas_i, infeas_j = np.nonzero(~feas_k)
    up[infeas_i * n_jk + infeas_j] = 0.0

    res = milp(
        c=obj,
        constraints=cons,
        integrality=np.ones(nvars),
        bounds=Bounds(np.zeros(nvars), up),
        options={"time_limit": time_limit_s, "disp": verbose},
    )

    if res.x is None:
        raise RuntimeError(f"Solver returned no solution: {res.message}")

    X = res.x[:nX].reshape(n_i, n_jk)
    Y = res.x[nX:]
    open_j = np.where(Y > 0.5)[0]

    assign_rows = []
    for i in range(n_i):
        j = int(np.argmax(X[i]))
        assign_rows.append({
            "norad_cat_id": clients.iloc[i].get("norad_cat_id"),
            "object_name": clients.iloc[i].get("object_name", ""),
            "client_lon_deg": clients.iloc[i].get("lon_deg"),
            "client_inc_deg": clients.iloc[i].get("inc_deg"),
            "depot_slot": int(slots_k.iloc[j]["slot_id"]),
            "depot_lon_deg": float(slots_k.iloc[j]["lon_deg"]),
            "depot_dalt_km": float(slots_k.iloc[j]["dalt_km"]),
            "depot_inc_deg": float(slots_k.iloc[j]["inc_deg"]),
            "roundtrip_dv_kms": float(dv_k[i, j]),
            "roundtrip_prop_kg": float(c_k[i, j]),
        })
    assignment = pd.DataFrame(assign_rows)

    dep_rows = []
    for j in open_j:
        members = assignment[assignment["depot_slot"]
                             == int(slots_k.iloc[j]["slot_id"])]
        wet = (p.m_d_dry
               + p.D * (c_k[:, j] * X[:, j]).sum()
               + p.D * p.m_s_L * X[:, j].sum())
        dep_rows.append({
            "depot_slot": int(slots_k.iloc[j]["slot_id"]),
            "lon_deg": float(slots_k.iloc[j]["lon_deg"]),
            "dalt_km": float(slots_k.iloc[j]["dalt_km"]),
            "a_km": float(slots_k.iloc[j]["a_km"]),
            "inc_deg": float(slots_k.iloc[j]["inc_deg"]),
            "n_clients": int(len(members)),
            "depot_wet_kg": float(wet),
            "depot_emleo_kg": float(
                p.m_d_dry * phi_k[j]
                + (p.D * (c_k[:, j] + p.m_s_L) * X[:, j]).sum() * phi_k[j]),
        })
    depot_summary = (pd.DataFrame(dep_rows)
                     .sort_values("lon_deg")
                     .reset_index(drop=True))

    return OFLPSolution(
        status=res.message,
        objective_emleo_kg=float(res.fun),
        open_slots=[int(slots_k.iloc[j]["slot_id"]) for j in open_j],
        assignment=assignment,
        depot_summary=depot_summary,
        solve_time_s=time.time() - t_start,
        n_vars=nvars,
    )
