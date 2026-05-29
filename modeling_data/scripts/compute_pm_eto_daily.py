"""
compute_pm_eto_daily.py — Download NASA POWER daily climate, compute FAO-56 PM ETo.
Replicates the R nasapower_PM_ETo.R formula exactly (same parameters, same algebra).
Saves: modeling_data/data/nasapower_ETo_daily_py.csv  (State, Year, DOY, ETo_mm)
Run from /home/user/LT-trait:
    python3 modeling_data/scripts/compute_pm_eto_daily.py
"""

import requests
import pandas as pd
import numpy as np
import os, time

SITES = {
    "AR": {"lon": -94.17, "lat": 36.08, "elev_m": 381},
    "MO": {"lon": -92.33, "lat": 38.95, "elev_m": 274},
    "NE": {"lon": -96.70, "lat": 40.82, "elev_m": 357},
}
START = "19900101"
END   = "20191231"
PARS  = "T2M_MAX,T2M_MIN,RH2M,WS2M,ALLSKY_SFC_SW_DWN"
URL   = "https://power.larc.nasa.gov/api/temporal/daily/point"
OUT   = "modeling_data/data/nasapower_ETo_daily_py.csv"


def download_site(state, lon, lat):
    print(f"  Downloading {state} ({lon}, {lat}) 1990–2019 …", flush=True)
    for attempt in range(4):
        try:
            r = requests.get(URL, params={
                "parameters": PARS, "community": "AG",
                "longitude": lon, "latitude": lat,
                "start": START, "end": END, "format": "JSON",
            }, timeout=180)
            r.raise_for_status()
            break
        except Exception as e:
            if attempt == 3:
                raise
            wait = 2 ** attempt
            print(f"    retry in {wait}s ({e})")
            time.sleep(wait)

    data = r.json()["properties"]["parameter"]
    dates = sorted(data["T2M_MAX"].keys())
    df = pd.DataFrame({"date": dates})
    for par in PARS.split(","):
        df[par] = df["date"].map(data[par])
    df["date"]  = pd.to_datetime(df["date"], format="%Y%m%d")
    df["Year"]  = df["date"].dt.year
    df["DOY"]   = df["date"].dt.dayofyear
    df["State"] = state
    # Replace NASA POWER fill value with NaN
    num_cols = PARS.split(",")
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
        df.loc[df[c] < -900, c] = np.nan
    return df


def calc_eto(df, lat, elev_m):
    LAT = lat
    ALT = elev_m
    tmed = (df["T2M_MAX"] + df["T2M_MIN"]) / 2.0
    WS   = df["WS2M"]
    Rs   = df["ALLSKY_SFC_SW_DWN"]
    RH   = df["RH2M"]
    DOY  = df["DOY"].values.astype(float)

    # Solar geometry (FAO-56 / R script)
    decl = 23.45 * np.sin(((360 * (DOY - 80)) / 365) * (np.pi / 180))
    hn   = np.degrees(np.arccos(
               np.clip(-np.tan(LAT * np.pi/180) * np.tan(decl * np.pi/180), -1, 1)
           ))
    dD   = 1 + 0.033 * np.cos((DOY * 360 / 365) * (np.pi / 180))
    Qo   = 37.6 * dD * (
        (np.pi/180) * hn * np.sin(LAT*np.pi/180) * np.sin(decl*np.pi/180) +
        np.cos(LAT*np.pi/180) * np.cos(decl*np.pi/180) * np.sin(hn*np.pi/180)
    )

    # Vapour pressures
    es  = 0.6108 * np.exp(17.27 * tmed / (tmed + 237.3))
    ea  = es * (RH / 100.0)
    vpd = es - ea

    # Net radiation
    BOC = Rs * (1 - 0.23)
    Qgo = (0.75 + 2e-5 * ALT) * Qo
    BOL = -4.903e-9 * (
        (tmed + 273)**4 *
        (0.34 - 0.14 * np.sqrt(ea.clip(lower=0))) *
        (1.35 * (Rs / np.where(Qgo == 0, np.nan, Qgo)) - 0.35)
    )
    Rn  = BOC + BOL

    # Slope of saturation vapour pressure curve
    s_v = (4098 * es) / (tmed + 237.3)**2

    # FAO-56 PM ETo (γ = 0.062 kPa °C⁻¹, fixed)
    ETo = (0.408 * Rn * s_v + 0.062 * 900 * WS * vpd / (tmed + 273)) / \
          (s_v + 0.062 * (1 + 0.34 * WS))
    ETo = np.maximum(0, ETo)

    out = df[["State", "Year", "DOY"]].copy()
    out["ETo_mm"] = ETo.values
    return out


if __name__ == "__main__":
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    all_dfs = []
    for state, site in SITES.items():
        raw = download_site(state, site["lon"], site["lat"])
        eto = calc_eto(raw, site["lat"], site["elev_m"])
        all_dfs.append(eto)

    result = pd.concat(all_dfs, ignore_index=True)
    result = result[result["Year"].between(1990, 2019)].dropna(subset=["ETo_mm"])
    result = result.sort_values(["State", "Year", "DOY"]).reset_index(drop=True)

    result.to_csv(OUT, index=False)
    print(f"\nSaved: {OUT}  ({len(result):,} rows)")
    print(result.groupby("State")["ETo_mm"].describe().round(2))
