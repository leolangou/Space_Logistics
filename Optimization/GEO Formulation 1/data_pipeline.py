# %% STEP 0 — Imports and user configuration

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# -------------------------------
# Input and output configuration
# -------------------------------

STATUS_FILENAME = "us_geo_catalog_status.csv"
ACTIVE_CATALOG_FILENAME = "celestrak_active_gpz_catalog.csv"
HISTORY_FILENAME = "gp_history_20260801_20260901_3cbd5060be.json"

# Leave as None to let the script search ./upload, the script directory, and cwd.
DATA_DIRECTORY = Path("data/geo_us")
OUTPUT_DIRECTORY = Path("geo_fuel_station_outputs")

# Leave as None to use the latest epoch found in the supplied history.
ANALYSIS_EPOCH_UTC: str | None = None


# ------------------------------------------
# Mission assumptions that YOU must replace
# ------------------------------------------

DEMO_MODE = True

# Option A: provide a CSV containing NORAD_CAT_ID and DEMAND_KG.
DEMAND_CSV: Path | None = None
DEMAND_ID_COLUMN = "NORAD_CAT_ID"
DEMAND_VALUE_COLUMN = "DEMAND_KG"

# Option B: demonstration-only uniform demand used when DEMAND_CSV is None.
UNIFORM_TARGET_DEMAND_KG = 100.0

# Capacity available to the reusable SSC associated with each station.
# The capacity constraint is sum(q_i assigned to station j) <= C.
# Total propellant stored at each fuel station over the planning horizon.
STATION_INVENTORY_KG = 2_500.0

# Propellant the SSC can carry and deliver on one station -> target -> station sortie.
SSC_PAYLOAD_CAPACITY_KG = 200.0

# If not None, a target is accessible only when its round-trip sortie delta-v
# is at or below this limit. Set this from the SSC fuel budget.
MAX_SORTIE_DV_KM_S: float | None = 0.70

# Circular GEO phasing model. One phasing revolution is the fastest/highest-dv
# option. More revolutions reduce delta-v but increase transfer time.
PHASING_REVOLUTIONS = 1

# The solver minimizes either total round-trip delta-v or a rocket-equation
# propellant estimate. "propellant_kg" requires dry mass and Isp assumptions.
OBJECTIVE_MODE = "delta_v"  # "delta_v" or "propellant_kg"
SSC_DRY_MASS_KG = 1_000.0
SSC_ISP_SECONDS = 320.0


# ------------------------
# Optimization configuration
# ------------------------

RANDOM_SEED = 27
MAX_STATIONS_TO_TEST = 10

# Fast first pass used to identify the smallest feasible station count.
SCREEN_POPULATION = 16
SCREEN_ITERATIONS = 40
SCREEN_RESTARTS = 1

# Longer second pass used to minimize cost for the selected station count.
FINAL_POPULATION = 24
FINAL_ITERATIONS = 80
FINAL_RESTARTS = 2

# Greedy allocation is repeated with slightly different target orders; the best
# capacity-feasible allocation is kept.
GREEDY_ALLOCATION_ATTEMPTS = 6

# These defaults are sized for an interactive first run. For a publication
# trade study, increase them (for example 40/200/5 and 60/500/10) and confirm
# that independent seeds converge to the same layouts and objective values.

# GSA constants.
GSA_G0 = 100.0
GSA_ALPHA = 20.0

# Large lexicographic penalties: feasibility dominates cost.
UNASSIGNED_TARGET_PENALTY = 1_000_000.0
CAPACITY_EXCESS_PENALTY = 100_000.0


# %% STEP 1 — Locate and load the three supplied files


def locate_input(filename: str) -> Path:
    """Find an input file without hard-coding the temporary workspace path."""
    candidates: list[Path] = []
    if DATA_DIRECTORY is not None:
        candidates.append(Path(DATA_DIRECTORY) / filename)

    script_dir = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    candidates.extend(
        [
            Path.cwd() / "upload" / filename,
            script_dir / "upload" / filename,
            Path.cwd() / filename,
            script_dir / filename,
        ]
    )
    for path in candidates:
        if path.exists():
            return path.resolve()

    searched = "\n".join(f"  - {p}" for p in candidates)
    raise FileNotFoundError(f"Could not find {filename}. Searched:\n{searched}")


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    status_path = locate_input(STATUS_FILENAME)
    active_path = locate_input(ACTIVE_CATALOG_FILENAME)
    history_path = locate_input(HISTORY_FILENAME)

    status = pd.read_csv(status_path)
    active = pd.read_csv(active_path)
    with history_path.open("r", encoding="utf-8") as handle:
        history = pd.DataFrame(json.load(handle))

    status["NORAD_CAT_ID"] = status["NORAD_CAT_ID"].astype(str)
    active["NORAD_CAT_ID"] = active["NORAD_CAT_ID"].astype(str)
    history["NORAD_CAT_ID"] = history["NORAD_CAT_ID"].astype(str)
    history["EPOCH"] = pd.to_datetime(history["EPOCH"], utc=True)

    print(f"Status targets: {len(status):,}")
    print(f"Active GPZ catalog objects: {len(active):,}")
    print(f"History records: {len(history):,}")
    print(f"History satellites: {history['NORAD_CAT_ID'].nunique():,}")
    return status, active, history


# %% STEP 2 — Keep the latest public orbital state for every catalog target


ORBITAL_NUMERIC_COLUMNS = [
    "MEAN_MOTION",
    "ECCENTRICITY",
    "INCLINATION",
    "RA_OF_ASC_NODE",
    "ARG_OF_PERICENTER",
    "MEAN_ANOMALY",
    "SEMIMAJOR_AXIS",
]


def latest_states_for_catalog(
    status: pd.DataFrame, active: pd.DataFrame, history: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    status = status.copy()
    active_ids = set(active["NORAD_CAT_ID"])
    status["IN_ACTIVE_GPZ_CATALOG"] = status["NORAD_CAT_ID"].isin(active_ids)
    print(
        "Status targets also present in active GPZ catalog: "
        f"{status['IN_ACTIVE_GPZ_CATALOG'].sum():,}/{len(status):,}"
    )

    catalog_ids = set(status["NORAD_CAT_ID"])
    matching = history[history["NORAD_CAT_ID"].isin(catalog_ids)].copy()

    for column in ORBITAL_NUMERIC_COLUMNS:
        matching[column] = pd.to_numeric(matching[column], errors="coerce")

    required = [
        "EPOCH",
        "MEAN_MOTION",
        "ECCENTRICITY",
        "INCLINATION",
        "RA_OF_ASC_NODE",
        "ARG_OF_PERICENTER",
        "MEAN_ANOMALY",
        "SEMIMAJOR_AXIS",
    ]
    matching = matching.dropna(subset=required)

    latest = (
        matching.sort_values("EPOCH")
        .groupby("NORAD_CAT_ID", as_index=False)
        .tail(1)
        .copy()
    )

    # Add current operating status without overwriting the OMM orbital fields.
    catalog_metadata = status[
        ["NORAD_CAT_ID", "OPS_STATUS_CODE", "IN_ACTIVE_GPZ_CATALOG"]
    ].drop_duplicates("NORAD_CAT_ID")
    latest = latest.merge(catalog_metadata, on="NORAD_CAT_ID", how="left")

    included_ids = set(latest["NORAD_CAT_ID"])
    excluded = status[~status["NORAD_CAT_ID"].isin(included_ids)].copy()
    excluded["EXCLUSION_REASON"] = "No complete public GP/OMM state in history file"

    print(f"Usable targets: {len(latest):,}")
    print(f"Excluded targets: {len(excluded):,}")
    return latest.reset_index(drop=True), excluded.reset_index(drop=True)


# %% STEP 3 — Propagate mean elements to one common epoch


MU_EARTH_KM3_S2 = 398_600.4418
GEO_RADIUS_KM = 42_164.0
SIDEREAL_DAY_S = 86_164.0905
STANDARD_GRAVITY_M_S2 = 9.80665


def wrap_to_pi(angle_rad: np.ndarray | float) -> np.ndarray:
    """Wrap radians to [-pi, pi)."""
    return (np.asarray(angle_rad) + np.pi) % (2.0 * np.pi) - np.pi


def solve_kepler(mean_anomaly_rad: np.ndarray, eccentricity: np.ndarray) -> np.ndarray:
    """Vectorized Newton solve for E - e sin(E) = M."""
    mean_anomaly_rad = np.mod(mean_anomaly_rad, 2.0 * np.pi)
    eccentric_anomaly = mean_anomaly_rad.copy()
    for _ in range(20):
        residual = eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly_rad
        derivative = 1.0 - eccentricity * np.cos(eccentric_anomaly)
        step = residual / derivative
        eccentric_anomaly -= step
        if np.max(np.abs(step)) < 1e-13:
            break
    return eccentric_anomaly


def julian_date(timestamp: pd.Timestamp) -> float:
    """UTC timestamp to Julian date; adequate for the GMST approximation below."""
    return timestamp.timestamp() / 86_400.0 + 2_440_587.5


def greenwich_mean_sidereal_time(timestamp: pd.Timestamp) -> float:
    """Return GMST in radians using the standard low-order expression."""
    jd = julian_date(timestamp)
    centuries = (jd - 2_451_545.0) / 36_525.0
    gmst_deg = (
        280.46061837
        + 360.98564736629 * (jd - 2_451_545.0)
        + 0.000387933 * centuries**2
        - centuries**3 / 38_710_000.0
    )
    return np.deg2rad(gmst_deg % 360.0)


def propagate_mean_elements_to_eci(
    states: pd.DataFrame, analysis_epoch: pd.Timestamp
) -> np.ndarray:
    """
    Propagate two-body mean elements and return ECI position vectors in km.

    TLE/OMM elements are SGP4 mean elements, so this is a transparent research
    prototype rather than a replacement for an operational SGP4 propagator.
    Every target is nevertheless compared at exactly the same epoch.
    """
    dt_seconds = (analysis_epoch - states["EPOCH"]).dt.total_seconds().to_numpy()
    mean_motion_rad_s = states["MEAN_MOTION"].to_numpy() * 2.0 * np.pi / 86_400.0
    mean_anomaly = np.deg2rad(states["MEAN_ANOMALY"].to_numpy()) + mean_motion_rad_s * dt_seconds

    eccentricity = states["ECCENTRICITY"].to_numpy()
    eccentric_anomaly = solve_kepler(mean_anomaly, eccentricity)
    true_anomaly = 2.0 * np.arctan2(
        np.sqrt(1.0 + eccentricity) * np.sin(eccentric_anomaly / 2.0),
        np.sqrt(1.0 - eccentricity) * np.cos(eccentric_anomaly / 2.0),
    )

    semimajor_axis = states["SEMIMAJOR_AXIS"].to_numpy()
    radius = semimajor_axis * (1.0 - eccentricity * np.cos(eccentric_anomaly))
    arg_latitude = np.deg2rad(states["ARG_OF_PERICENTER"].to_numpy()) + true_anomaly
    inclination = np.deg2rad(states["INCLINATION"].to_numpy())
    raan = np.deg2rad(states["RA_OF_ASC_NODE"].to_numpy())

    cos_u, sin_u = np.cos(arg_latitude), np.sin(arg_latitude)
    cos_o, sin_o = np.cos(raan), np.sin(raan)
    cos_i, sin_i = np.cos(inclination), np.sin(inclination)

    x = radius * (cos_o * cos_u - sin_o * sin_u * cos_i)
    y = radius * (sin_o * cos_u + cos_o * sin_u * cos_i)
    z = radius * (sin_u * sin_i)
    return np.column_stack([x, y, z])


def eci_to_ecef(eci_km: np.ndarray, analysis_epoch: pd.Timestamp) -> np.ndarray:
    theta = greenwich_mean_sidereal_time(analysis_epoch)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    x = cos_t * eci_km[:, 0] + sin_t * eci_km[:, 1]
    y = -sin_t * eci_km[:, 0] + cos_t * eci_km[:, 1]
    return np.column_stack([x, y, eci_km[:, 2]])


def prepare_target_locations(states: pd.DataFrame) -> tuple[pd.DataFrame, pd.Timestamp]:
    if ANALYSIS_EPOCH_UTC is None:
        analysis_epoch = states["EPOCH"].max()
    else:
        analysis_epoch = pd.Timestamp(ANALYSIS_EPOCH_UTC)
        if analysis_epoch.tzinfo is None:
            analysis_epoch = analysis_epoch.tz_localize("UTC")
        else:
            analysis_epoch = analysis_epoch.tz_convert("UTC")

    eci = propagate_mean_elements_to_eci(states, analysis_epoch)
    ecef = eci_to_ecef(eci, analysis_epoch)

    result = states.copy()
    result["ANALYSIS_EPOCH_UTC"] = analysis_epoch.isoformat()
    result["ECI_X_KM"] = eci[:, 0]
    result["ECI_Y_KM"] = eci[:, 1]
    result["ECI_Z_KM"] = eci[:, 2]
    result["ECEF_X_KM"] = ecef[:, 0]
    result["ECEF_Y_KM"] = ecef[:, 1]
    result["ECEF_Z_KM"] = ecef[:, 2]
    result["LONGITUDE_RAD"] = np.arctan2(ecef[:, 1], ecef[:, 0])
    result["LONGITUDE_DEG"] = np.rad2deg(result["LONGITUDE_RAD"])
    result["GEOCENTRIC_LATITUDE_DEG"] = np.rad2deg(
        np.arctan2(ecef[:, 2], np.hypot(ecef[:, 0], ecef[:, 1]))
    )

    print(f"Common analysis epoch: {analysis_epoch.isoformat()}")
    return result, analysis_epoch


# %% STEP 4 — Attach each target's refueling demand q_i


def attach_demands(targets: pd.DataFrame) -> pd.DataFrame:
    result = targets.copy()

    if DEMAND_CSV is not None:
        demand_path = Path(DEMAND_CSV)
        demand = pd.read_csv(demand_path)
        demand[DEMAND_ID_COLUMN] = demand[DEMAND_ID_COLUMN].astype(str)
        demand[DEMAND_VALUE_COLUMN] = pd.to_numeric(
            demand[DEMAND_VALUE_COLUMN], errors="raise"
        )
        result = result.merge(
            demand[[DEMAND_ID_COLUMN, DEMAND_VALUE_COLUMN]],
            left_on="NORAD_CAT_ID",
            right_on=DEMAND_ID_COLUMN,
            how="left",
        )
        result["DEMAND_KG"] = result[DEMAND_VALUE_COLUMN]
        if result["DEMAND_KG"].isna().any():
            missing = result.loc[result["DEMAND_KG"].isna(), "NORAD_CAT_ID"].tolist()
            raise ValueError(f"Demand file is missing NORAD IDs: {missing}")
    else:
        result["DEMAND_KG"] = float(UNIFORM_TARGET_DEMAND_KG)
        if DEMO_MODE:
            warnings.warn(
                "DEMO_MODE: using a uniform placeholder demand for every target. "
                "Replace this before reporting mission-design results.",
                stacklevel=2,
            )

    if (result["DEMAND_KG"] <= 0).any():
        raise ValueError("Every target demand must be positive.")
    if (result["DEMAND_KG"] > SSC_PAYLOAD_CAPACITY_KG).any():
        bad = result.loc[
            result["DEMAND_KG"] > SSC_PAYLOAD_CAPACITY_KG,
            ["NORAD_CAT_ID", "OBJECT_NAME", "DEMAND_KG"],
        ]
        raise ValueError(
            "At least one target requires more propellant than the SSC can carry "
            "on one sortie:\n" + bad.to_string(index=False)
        )
    return result


def run_pipeline() -> pd.DataFrame:
    """Convenience wrapper chaining STEP 1-4 exactly as the notebook cells do."""
    status, active, history = load_inputs()
    latest, excluded = latest_states_for_catalog(status, active, history)
    located, epoch = prepare_target_locations(latest)
    with_demand = attach_demands(located)
    return with_demand

# %%
