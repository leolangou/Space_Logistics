from __future__ import annotations

import hashlib
import io
import json
import os
import time
from getpass import getpass
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import spacetrack.operators as op
from sgp4 import omm
from sgp4.api import SGP4_ERRORS, Satrec, jday
from spacetrack import SpaceTrackClient

# -----------------------------
# USER SETTINGS
# -----------------------------
START_UTC = "2026-08-01T00:00:00Z"
END_UTC = "2026-09-01T00:00:00Z"
CADENCE = "6h"                 # e.g. "1h", "6h", "1D"
ELEMENT_SELECTION = "nearest"  # "nearest" or "past"
MAX_ELEMENT_AGE_DAYS = 7
US_OWNER_CODES = {"US"}        # add "USBZ" if you want joint US/Brazil objects

DATA_DIR = Path("data/geo_us")
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Avoid writing a literal full URL in one string; the pieces are easier to swap/test.
CELESTRAK_BASE = "https://" + "celestrak.org"
SATCAT_PATH = "/satcat/records.php"


def get_current_us_geo_catalog(cache_hours: float = 24.0) -> pd.DataFrame:
    """Return current active US GEO/geosynchronous payloads from CelesTrak SATCAT."""
    cache_file = DATA_DIR / "celestrak_active_gpz_catalog.csv"

    if cache_file.exists():
        age_hours = (time.time() - cache_file.stat().st_mtime) / 3600.0
        if age_hours < cache_hours:
            df = pd.read_csv(cache_file)
        else:
            df = _download_catalog(cache_file)
    else:
        df = _download_catalog(cache_file)

    # CelesTrak SATCAT uses OWNER="US" for United States.
    df = df[df["OWNER"].isin(US_OWNER_CODES)].copy()
    df = df[df["OBJECT_TYPE"].eq("PAY")].copy()

    # Space-Track's geosynchronous report convention is 1430--1450 minutes.
    # SPECIAL=GPZ already narrows to the GEO protected zone; this is an extra guard.
    if "PERIOD" in df.columns:
        df["PERIOD"] = pd.to_numeric(df["PERIOD"], errors="coerce")
        df = df[df["PERIOD"].between(1430.0, 1450.0, inclusive="both")].copy()

    df["NORAD_CAT_ID"] = pd.to_numeric(df["NORAD_CAT_ID"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["NORAD_CAT_ID"]).copy()
    df["NORAD_CAT_ID"] = df["NORAD_CAT_ID"].astype(int)

    return df.sort_values("NORAD_CAT_ID").reset_index(drop=True)


def _download_catalog(cache_file: Path) -> pd.DataFrame:
    params = {
        "SPECIAL": "gpz",
        "FORMAT": "CSV",
        "ACTIVE": "1",
    }
    response = requests.get(CELESTRAK_BASE + SATCAT_PATH, params=params, timeout=60)
    response.raise_for_status()
    df = pd.read_csv(io.StringIO(response.text))
    df.to_csv(cache_file, index=False)
    return df


def get_spacetrack_credentials() -> tuple[str, str]:
    username = os.getenv("SPACETRACK_USER")
    password = os.getenv("SPACETRACK_PASSWORD")

    if not username:
        username = input("Space-Track username/email: ").strip()
    if not password:
        password = getpass("Space-Track password: ")

    return username, password


def history_cache_path(ids: list[int], start: pd.Timestamp, end: pd.Timestamp) -> Path:
    digest = hashlib.sha1(
        ",".join(map(str, sorted(ids))).encode("utf-8")
    ).hexdigest()[:10]
    s = start.strftime("%Y%m%d")
    e = end.strftime("%Y%m%d")
    return DATA_DIR / f"gp_history_{s}_{e}_{digest}.json"


def download_gp_history(
    ids: list[int],
    start: pd.Timestamp,
    end: pd.Timestamp,
    buffer_days: int = 7,
) -> list[dict]:
    """Download GP_HISTORY once and cache it locally."""
    cache_file = history_cache_path(ids, start, end)
    if cache_file.exists():
        with cache_file.open("r", encoding="utf-8") as f:
            return json.load(f)

    username, password = get_spacetrack_credentials()

    query_start = (start - pd.Timedelta(days=buffer_days)).to_pydatetime().replace(tzinfo=None)
    query_end = (end + pd.Timedelta(days=buffer_days)).to_pydatetime().replace(tzinfo=None)

    with SpaceTrackClient(identity=username, password=password) as st:
        # No explicit format="json": the client then returns parsed JSON.
        records = st.gp_history(
            norad_cat_id=ids,
            epoch=op.inclusive_range(query_start, query_end),
            orderby=["NORAD_CAT_ID", "EPOCH"],
            timeout=120.0,
        )

    with cache_file.open("w", encoding="utf-8") as f:
        json.dump(records, f, indent=2)

    return records


def index_history(records: list[dict]) -> dict[int, tuple[np.ndarray, list[dict]]]:
    grouped: dict[int, list[tuple[pd.Timestamp, dict]]] = {}

    for rec in records:
        try:
            sat_id = int(rec["NORAD_CAT_ID"])
            epoch = pd.to_datetime(rec["EPOCH"], utc=True)
        except (KeyError, TypeError, ValueError):
            continue
        grouped.setdefault(sat_id, []).append((epoch, rec))

    indexed = {}
    for sat_id, items in grouped.items():
        items.sort(key=lambda x: x[0])
        epochs_ns = np.array([epoch.value for epoch, _ in items], dtype=np.int64)
        recs = [rec for _, rec in items]
        indexed[sat_id] = (epochs_ns, recs)

    return indexed


def select_element(
    epochs_ns: np.ndarray,
    records: list[dict],
    target: pd.Timestamp,
    mode: str,
) -> tuple[dict, pd.Timestamp] | tuple[None, None]:
    target_ns = target.value
    pos = int(np.searchsorted(epochs_ns, target_ns, side="left"))

    if mode == "past":
        idx = pos - 1 if pos > 0 else None
    elif mode == "nearest":
        candidates = []
        if pos > 0:
            candidates.append(pos - 1)
        if pos < len(epochs_ns):
            candidates.append(pos)
        if not candidates:
            return None, None
        idx = min(candidates, key=lambda i: abs(int(epochs_ns[i]) - target_ns))
    else:
        raise ValueError("ELEMENT_SELECTION must be 'nearest' or 'past'.")

    if idx is None:
        return None, None

    rec = records[idx]
    epoch = pd.to_datetime(rec["EPOCH"], utc=True)
    return rec, epoch


def propagate_omm(record: dict, target: pd.Timestamp):
    sat = Satrec()
    omm.initialize(sat, record)

    sec = target.second + target.microsecond / 1e6 + target.nanosecond / 1e9
    jd, fr = jday(
        target.year,
        target.month,
        target.day,
        target.hour,
        target.minute,
        sec,
    )
    error, r_km, v_km_s = sat.sgp4(jd, fr)
    return error, r_km, v_km_s


def build_uniform_timeseries(
    catalog: pd.DataFrame,
    records: list[dict],
    start: pd.Timestamp,
    end: pd.Timestamp,
    cadence: str,
) -> pd.DataFrame:
    history = index_history(records)
    targets = pd.date_range(start=start, end=end, freq=cadence, tz="UTC")

    name_map = catalog.set_index("NORAD_CAT_ID")["OBJECT_NAME"].to_dict()
    rows = []

    for sat_id in catalog["NORAD_CAT_ID"].tolist():
        if sat_id not in history:
            continue

        epochs_ns, sat_records = history[sat_id]

        for target in targets:
            rec, element_epoch = select_element(
                epochs_ns, sat_records, target, ELEMENT_SELECTION
            )
            if rec is None:
                continue

            age_hours = abs((target - element_epoch).total_seconds()) / 3600.0
            if age_hours > MAX_ELEMENT_AGE_DAYS * 24.0:
                continue

            try:
                error, r_km, v_km_s = propagate_omm(rec, target)
            except (KeyError, TypeError, ValueError) as exc:
                rows.append(
                    {
                        "time_utc": target,
                        "norad_cat_id": sat_id,
                        "object_name": name_map.get(sat_id, rec.get("OBJECT_NAME")),
                        "element_epoch_utc": element_epoch,
                        "element_age_hours": age_hours,
                        "sgp4_error": -1,
                        "sgp4_error_text": str(exc),
                    }
                )
                continue

            row = {
                "time_utc": target,
                "norad_cat_id": sat_id,
                "object_name": name_map.get(sat_id, rec.get("OBJECT_NAME")),
                "element_epoch_utc": element_epoch,
                "element_age_hours": age_hours,
                "sgp4_error": int(error),
                "sgp4_error_text": "" if error == 0 else SGP4_ERRORS.get(error, "unknown"),
                "x_teme_km": np.nan,
                "y_teme_km": np.nan,
                "z_teme_km": np.nan,
                "vx_teme_km_s": np.nan,
                "vy_teme_km_s": np.nan,
                "vz_teme_km_s": np.nan,
            }

            if error == 0:
                row.update(
                    {
                        "x_teme_km": r_km[0],
                        "y_teme_km": r_km[1],
                        "z_teme_km": r_km[2],
                        "vx_teme_km_s": v_km_s[0],
                        "vy_teme_km_s": v_km_s[1],
                        "vz_teme_km_s": v_km_s[2],
                    }
                )

            rows.append(row)

    return pd.DataFrame(rows).sort_values(["time_utc", "norad_cat_id"]).reset_index(drop=True)


def main():
    start = pd.Timestamp(START_UTC)
    end = pd.Timestamp(END_UTC)

    print("1) Loading current active US GEO catalog...")
    catalog = get_current_us_geo_catalog()
    print(f"   Selected {len(catalog)} satellites")
    print(catalog[["NORAD_CAT_ID", "OBJECT_NAME", "OWNER", "PERIOD"]].head(20).to_string(index=False))

    ids = catalog["NORAD_CAT_ID"].tolist()
    if not ids:
        raise RuntimeError("No satellites matched the US GEO filters.")

    print("\n2) Loading historical GP/OMM element sets...")
    records = download_gp_history(ids, start, end)
    print(f"   Loaded {len(records):,} historical element sets")

    print("\n3) Propagating to common timestamps...")
    ts = build_uniform_timeseries(catalog, records, start, end, CADENCE)

    out_csv = DATA_DIR / "us_geo_timeseries_teme.csv"
    ts.to_csv(out_csv, index=False)

    good = ts[ts["sgp4_error"].eq(0)].copy()
    n_times = good["time_utc"].nunique() if not good.empty else 0
    n_sats = good["norad_cat_id"].nunique() if not good.empty else 0

    print(f"   Valid propagated rows: {len(good):,}")
    print(f"   Time steps: {n_times:,}")
    print(f"   Satellites represented: {n_sats:,}")
    print(f"\nSaved: {out_csv}")

    # Example point cloud for one timestamp:
    if not good.empty:
        t0 = good["time_utc"].iloc[0]
        cloud = good.loc[
            good["time_utc"].eq(t0),
            ["norad_cat_id", "object_name", "x_teme_km", "y_teme_km", "z_teme_km"],
        ]
        print(f"\nExample point cloud at {t0}:")
        print(cloud.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
