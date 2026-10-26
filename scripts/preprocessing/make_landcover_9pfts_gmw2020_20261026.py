#!/usr/bin/env python3
# =============================================================================
#  MAKE_LANDCOVER_9PFTS_GMW2020  --  9-PFT land-cover files (no pastures, no crops) with the FIXED GMW-2020 mangrove area
#  Dr Sri Nandini-Weiss, ZMT Bremen.                             2026-10-26
#  Runs on Levante. Replaces the 9-PFT files of make_landcover_8_9pfts_20261023.py, which used the CGMD 1984/2023 maps (area changing
#  with time, made for the transient 1984-2023 experiment). The 8-PFT files are not affected and are not rewritten here.
#  Inputs  /pool/data/ICON/grids/public/mpim/0049/land/r0103/bc_land_frac_12pfts_YYYY.nc   (1850-2024)
#          /work/bm1548/u241304/preprocessing/PFT13/input_boundary_files/gmw_mng_2020_v4019_gtiff/bc_land_frac_13pfts_2020.nc
#          (GMW 2020 mangrove; fract_pft13 = fraction of the vegetated tile, identical in every year: same area in all years)
#  Rule    8 natural PFTs as in make_landcover_8_9pfts_20261023.py (C3 pasture + C3 crop -> PFT07; C4 pasture + C4 crop -> PFT08);
#          mangrove area taken from PFT01 and PFT02 in proportion to their cover, from all of PFT01-08 where those are too small.
#  Output  <outdir>/bc_land_frac_9pfts_YYYY.nc ; checks printed: PFT sum = 1 in every land cell, no negative fraction,
#          mangrove physical area (fract_pft13 * fract_land * fract_veg * cell fraction) identical in all years.
# =============================================================================
import os, sys
import numpy as np
from netCDF4 import Dataset
POOL = "/pool/data/ICON/grids/public/mpim/0049/land/r0103"
GMW = "/work/bm1548/u241304/preprocessing/PFT13/input_boundary_files/gmw_mng_2020_v4019_gtiff/bc_land_frac_13pfts_2020.nc"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
def write(src, dst, F, keep):
    with Dataset(src) as s, Dataset(dst, "w", format=s.data_model) as o:
        o.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
        o.history = "2026-10-26 make_landcover_9pfts_gmw2020_20261026.py: pastures and crops removed (area to the grass of the pathway), fixed GMW-2020 mangrove; " + getattr(s, "history", "")
        for n, dim in s.dimensions.items(): o.createDimension(n, len(dim))
        for n, v in s.variables.items():
            if n.startswith("fract_pft") and n not in keep: continue
            w = o.createVariable(n, v.dtype, v.dimensions, fill_value=getattr(v, "_FillValue", None))
            w.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"})
            w[:] = F[n] if n in F else v[:]
        if "fract_pft13" in keep and "fract_pft13" not in s.variables:
            v = s.variables["fract_pft01"]
            w = o.createVariable("fract_pft13", v.dtype, v.dimensions, fill_value=getattr(v, "_FillValue", None))
            w.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"}); w.long_name = "mangrove (GMW 2020)"
            w[:] = F["fract_pft13"]
with Dataset(GMW) as m:
    m.set_auto_mask(False); f13 = np.clip(np.array(m["fract_pft13"][:], dtype=float), 0, 1)
bad = 0; area0 = None; amax = 0.0
for y in range(1850, 2025):
    src = f"{POOL}/bc_land_frac_12pfts_{y}.nc"
    with Dataset(src) as d:
        d.set_auto_mask(False)
        f = {i: np.array(d[f"fract_pft{i:02d}"][:], dtype=float) for i in range(1, 13)}
        land = np.array(d["notsea"][:]) > 0
        cap = np.array(d["fract_land"][:], dtype=float) * np.array(d["fract_veg"][:], dtype=float)
    g = {i: f[i].copy() for i in range(1, 9)}
    g[7] += f[9] + f[11]; g[8] += f[10] + f[12]
    h = {i: g[i].copy() for i in range(1, 9)}
    p12 = h[1] + h[2]; take12 = np.minimum(f13, p12)
    for i in (1, 2): h[i] -= np.where(p12 > 0, take12 * g[i] / np.where(p12 > 0, p12, 1), 0.0)
    rest = f13 - take12; tot = sum(h.values())
    for i in range(1, 9): h[i] -= np.where(tot > 0, rest * h[i] / np.where(tot > 0, tot, 1), 0.0)
    for i in range(1, 9): h[i] = np.maximum(h[i], 0.0)
    s9 = sum(h.values()) + f13
    e9 = np.abs(s9 - 1)[land].max(); neg = min(v[land].min() for v in h.values())
    if e9 > 1e-9 or neg < -1e-12: bad += 1
    area = f13 * cap
    if area0 is None: area0 = area.copy()
    amax = max(amax, float(np.abs(area - area0).max()))
    F9 = {f"fract_pft{i:02d}": h[i] for i in range(1, 9)}; F9["fract_pft13"] = f13
    write(src, f"{OUT}/bc_land_frac_9pfts_{y}.nc", F9, {f"fract_pft{i:02d}" for i in range(1, 9)} | {"fract_pft13"})
    if y in (1850, 1979, 2024): print(f"{y}: max |sum-1| {e9:.1e}; smallest fraction {neg:.1e}; mangrove cells {(f13[land] > 0).sum()}")
print("years failing a check:", bad, "| max change of mangrove area (fraction of cell) between years:", amax)
