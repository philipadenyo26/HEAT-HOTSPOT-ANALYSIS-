"""
Step 1 (Python). ERA5 downloads from the Copernicus Climate Data Store.
One request per month with valid dates only, sent one at a time (CDS rule).
Files already downloaded are reused. Failed requests print the full CDS error.

Produces:
  data/raw/era5_tcwv_by_scene.csv       -> upload to GEE as table asset 'era5_tcwv_scenes'
  data/raw/era5_tmax_dry_2020_2023.tif  -> used in R for the LST validation
Run from the thesis_analysis folder:  python python/01_era5_download.py
"""
import os, glob
import cdsapi
import pandas as pd
import xarray as xr
import rioxarray  # noqa: F401

AREA = [6.25, -0.65, 5.35, 0.80]          # North, West, South, East
DRY = ["01", "02", "03", "11", "12"]
DAYS = [f"{d:02d}" for d in range(1, 32)]
os.makedirs("data/raw/era5", exist_ok=True)

scenes = pd.read_csv("data/gee/landsat_scene_list.csv", dtype={"hour_utc": str})
scenes["date"] = pd.to_datetime(scenes["date"])
hours = sorted({f"{int(h):02d}:00" for h in scenes["hour_utc"]})

import calendar

def valid_days(year, month):
    return [f"{d:02d}" for d in range(1, calendar.monthrange(year, int(month))[1] + 1)]

def get_tcwv(year, month):
    target = f"data/raw/era5/tcwv_{year}-{month}.nc"
    if os.path.exists(target):
        return f"tcwv {year}-{month}: already done"
    days = sorted(scenes.loc[(scenes["date"].dt.year == year) & (scenes["date"].dt.strftime("%m") == month),
                             "date"].dt.strftime("%d").unique().tolist())
    try:
        cdsapi.Client().retrieve("reanalysis-era5-single-levels", {
            "product_type": ["reanalysis"], "variable": ["total_column_water_vapour"],
            "year": [str(year)], "month": [month], "day": days, "time": hours,
            "area": AREA, "data_format": "netcdf", "download_format": "unarchived"}, target)
        return f"tcwv {year}-{month}: downloaded"
    except Exception as e:
        return f"tcwv {year}-{month}: FAILED -> {e}"

def get_tmax(year, month):
    target = f"data/raw/era5/tmax_{year}_{month}.nc"
    if os.path.exists(target):
        return f"tmax {year}_{month}: already done"
    try:
        cdsapi.Client().retrieve("derived-era5-single-levels-daily-statistics", {
            "product_type": "reanalysis", "variable": ["2m_temperature"],
            "year": str(year), "month": [month], "day": valid_days(year, month),
            "daily_statistic": "daily_maximum", "time_zone": "utc+00:00",
            "frequency": "1_hourly", "area": AREA}, target)
        return f"tmax {year}_{month}: downloaded"
    except Exception as e:
        return f"tmax {year}_{month}: FAILED -> {e}"

ym = sorted({(d.year, d.strftime("%m")) for d in scenes["date"]})
jobs = [(get_tcwv, y, m) for y, m in ym] + [(get_tmax, y, m) for y in range(2020, 2024) for m in DRY]
import time
print(len(jobs), "monthly requests to check or download, one at a time", flush=True)
for fn, y, m in jobs:                                      # sequential, as the CDS requires
    for attempt in range(20):
        msg = fn(y, m)
        print(msg, flush=True)
        if "Maximum number queued" in msg:
            print("  CDS queue full: waiting 5 minutes before retrying", flush=True)
            time.sleep(300)
            continue
        break

def load_series(pattern, var):
    parts = []
    for f in sorted(glob.glob(pattern)):
        ds = xr.open_dataset(f)
        t = "valid_time" if "valid_time" in ds.dims else "time"
        s = ds[var].mean(dim=["latitude", "longitude"]).to_series()
        s.index = pd.to_datetime(s.index.get_level_values(t) if isinstance(s.index, pd.MultiIndex) else s.index)
        parts.append(s)
    s = pd.concat(parts).sort_index()
    return s[~s.index.duplicated()]

# 1. Water vapour per scene (g/cm2)
tcwv = load_series("data/raw/era5/tcwv_*.nc", "tcwv")
rows = []
for _, s in scenes.iterrows():
    stamp = s["date"] + pd.Timedelta(hours=int(s["hour_utc"]))
    nearest = tcwv.index[abs(tcwv.index - stamp).argmin()]
    rows.append({"scene_id": s["scene_id"], "date": s["date"].strftime("%Y-%m-%d"),
                 "w": round(float(tcwv.loc[nearest]) / 10.0, 4)})
pd.DataFrame(rows).to_csv("data/raw/era5_tcwv_by_scene.csv", index=False)
print("Water vapour written for", len(rows), "scenes")

# 2. Dry-season mean daily maximum 2 m temperature (deg C), 2020 to 2023
grids = []
for f in sorted(glob.glob("data/raw/era5/tmax_*.nc")):
    ds = xr.open_dataset(f)
    t = "valid_time" if "valid_time" in ds.dims else "time"
    grids.append(ds["t2m"].rename({t: "time"}))
t2m = xr.concat(grids, dim="time")
t2m = t2m.sel(time=~t2m.get_index("time").duplicated())
dry_mean = (t2m.mean(dim="time") - 273.15).rename("tmax_c")
dry_mean = dry_mean.rio.write_crs("EPSG:4326").rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
dry_mean.rio.to_raster("data/raw/era5_tmax_dry_2020_2023.tif")
print("Dry-season mean daily maximum written")
