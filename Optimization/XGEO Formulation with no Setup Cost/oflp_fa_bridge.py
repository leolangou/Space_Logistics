"""
Bridge from data_pipeline.py (the partner's loader/propagator/demand pipeline)
into oflp_fa.py (the exact FA solver matching that pipeline's mission model).

Typical use, once STATUS_FILENAME/ACTIVE_CATALOG_FILENAME/HISTORY_FILENAME
point at real files in DATA_DIRECTORY:

    import data_pipeline as dp
    from oflp_fa import FAParams
    from oflp_fa_bridge import clients_from_pipeline, run_full_sweep

    targets = dp.run_pipeline()
    clients = clients_from_pipeline(targets)
    params = params_from_pipeline_config()
    summary, solutions = run_full_sweep(clients, params)
    display(summary)

    best_p = ...                      # pick from the marginal-savings curve
    display(solutions[best_p].station_summary)
    display(solutions[best_p].assignment)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import data_pipeline as dp
from oflp_geo import make_candidate_slots, A_GEO
from oflp_fa import FAParams, sweep_stations


def clients_from_pipeline(targets: pd.DataFrame) -> pd.DataFrame:
    """
    Map data_pipeline.py's output columns onto the schema oflp_fa expects:
    norad_cat_id, object_name, a_km, ecc, inc_deg, lon_deg, demand_kg.

    Requires prepare_target_locations() and attach_demands() to have already
    run (i.e. call dp.run_pipeline(), or the individual STEP functions, first).
    """
    required = ["NORAD_CAT_ID", "SEMIMAJOR_AXIS", "ECCENTRICITY",
               "INCLINATION", "LONGITUDE_DEG", "DEMAND_KG"]
    missing = [c for c in required if c not in targets.columns]
    if missing:
        raise KeyError(
            f"targets is missing {missing}. Did you run prepare_target_locations() "
            "and attach_demands() before calling clients_from_pipeline()?")

    out = pd.DataFrame({
        "norad_cat_id": targets["NORAD_CAT_ID"],
        "object_name": targets.get("OBJECT_NAME", ""),
        "a_km": pd.to_numeric(targets["SEMIMAJOR_AXIS"]),
        "ecc": pd.to_numeric(targets["ECCENTRICITY"]),
        "inc_deg": pd.to_numeric(targets["INCLINATION"]),
        "lon_deg": pd.to_numeric(targets["LONGITUDE_DEG"]) % 360.0,
        "demand_kg": pd.to_numeric(targets["DEMAND_KG"]),
    })
    n0 = len(out)
    out = out.dropna(subset=["a_km", "inc_deg", "lon_deg", "demand_kg"])
    print(f"[oflp_fa_bridge] clients: {len(out)} kept of {n0}")
    return out.reset_index(drop=True)


def params_from_pipeline_config() -> FAParams:
    """Read the mission-assumption constants directly out of data_pipeline.py's
    STEP 0 configuration block, so the solver always matches whatever the
    notebook is currently set to -- no numbers re-typed by hand."""
    return FAParams(
        station_inventory_kg=dp.STATION_INVENTORY_KG,
        ssc_payload_capacity_kg=dp.SSC_PAYLOAD_CAPACITY_KG,
        max_sortie_dv_kms=dp.MAX_SORTIE_DV_KM_S,
        phasing_revolutions=dp.PHASING_REVOLUTIONS,
        objective_mode=dp.OBJECTIVE_MODE,
        ssc_dry_mass_kg=dp.SSC_DRY_MASS_KG,
        ssc_isp_s=dp.SSC_ISP_SECONDS,
        max_stations_to_test=dp.MAX_STATIONS_TO_TEST,
    )


def make_slots_for_clients(clients: pd.DataFrame, lon_step_deg: float = 10.0,
                           dalt_km=(0.0,), inc_pad_deg: float = 1.0,
                           inc_step_deg: float = 1.0) -> pd.DataFrame:
    """
    Candidate station slots. The pipeline doesn't specify a discretization
    (it searches station location in continuous space via GSA), so this
    builds a grid spanning the client set: full longitude ring, and
    inclination from 0 up to the client maximum (+ padding), stepped by
    inc_step_deg. Tune this before trusting results -- a coarser grid finds
    worse locations, a finer one costs more solve time per p.
    """
    max_inc = float(clients["inc_deg"].max()) + inc_pad_deg
    inc_deg = tuple(np.arange(0.0, max_inc + inc_step_deg, inc_step_deg))
    return make_candidate_slots(
        lon_deg=np.arange(0.0, 360.0, lon_step_deg),
        dalt_km=dalt_km,
        inc_deg=inc_deg,
    )


def run_full_sweep(clients: pd.DataFrame, params: FAParams | None = None,
                   slots: pd.DataFrame | None = None,
                   time_limit_s: float = 60.0):
    """Build slots (if not supplied), sweep station count, return results."""
    p = params or params_from_pipeline_config()
    s = slots if slots is not None else make_slots_for_clients(clients)

    total_demand = clients["demand_kg"].sum()
    print(f"[oflp_fa_bridge] clients={len(clients)} total_demand={total_demand:.0f} kg "
         f"| station_inventory={p.station_inventory_kg:.0f} kg "
         f"| objective_mode={p.objective_mode!r} | phasing_revs={p.phasing_revolutions} "
         f"| max_sortie_dv={p.max_sortie_dv_kms} km/s | candidate slots={len(s)}")

    summary, solutions = sweep_stations(clients, s, p, time_limit_s=time_limit_s)
    print(summary.to_string(index=False))
    return summary, solutions
