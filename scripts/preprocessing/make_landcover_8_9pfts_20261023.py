#!/usr/bin/env python3
# =============================================================================
#  MAKE_LANDCOVER_8_9PFTS  --  land-cover files without pastures and crops, for a dedicated 8-PFT and 9-PFT use case
#  Dr Sri Nandini-Weiss, ZMT Bremen.                             2026-10-23
#  Runs on the Mac on COPIES (input/), writes output/8pfts and output/9pfts; touches nothing on Levante.
#  Inputs  input/pool_r0103/bc_land_frac_12pfts_YYYY.nc  (pool 0049 r0103, 1850-2024: natural PFTs ESA CCI 2022, crops/pasture LUH3)
#          input/mangrove_13pfts/bc_land_frac_13pfts_YYYY.nc (CGMD mangrove, 1984-2023; fract_pft13 only is used)
#  Rule 8 PFTs (revised 26 Oct 2026, decision of the author): all pasture and crop area goes to grass of its pathway:
#          C3 pasture + C3 crop (PFT09 + PFT11) -> C3 grass PFT07; C4 pasture + C4 crop (PFT10 + PFT12) -> C4 grass PFT08.
#  Rule 9 PFTs: the 8-PFT map plus fract_pft13 of the mangrove file (1850-1983: 1984 map, 2024: 2023 map); the mangrove area is taken
#          from PFT01 and PFT02 in proportion to their cover, from all of PFT01-08 where those two are too small.
#  Checks printed per year: PFT fractions sum to 1 in every land cell; all fractions >= 0.
# =============================================================================
import os, shutil, sys
import numpy as np
from netCDF4 import Dataset
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C3NAT = [1, 2, 3, 4, 5, 6, 7]
def fr(d, i): return d[f"fract_pft{i:02d}"][:].astype(float).filled(0) if hasattr(d[f"fract_pft{i:02d}"][:], "filled") else d[f"fract_pft{i:02d}"][:].astype(float)
def write(src, dst, F, keep):
    with Dataset(src) as s, Dataset(dst, "w", format=s.data_model) as o:
        o.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
        o.history = "2026-10-23 make_landcover_8_9pfts_20261023.py: pastures and crops removed (area to the natural PFTs of the pathway); " + getattr(s, "history", "")
        for n, dim in s.dimensions.items(): o.createDimension(n, len(dim))
        for n, v in s.variables.items():
            if n.startswith("fract_pft") and n not in keep: continue
            w = o.createVariable(n, v.dtype, v.dimensions, fill_value=getattr(v, "_FillValue", None))
            w.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"})
            w[:] = F[n] if n in F else v[:]
        if "fract_pft13" in keep and "fract_pft13" not in s.variables:
            v = s.variables["fract_pft01"]
            w = o.createVariable("fract_pft13", v.dtype, v.dimensions, fill_value=getattr(v, "_FillValue", None))
            w.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"}); w.long_name = "mangrove (CGMD)"
            w[:] = F["fract_pft13"]
bad = 0
for y in range(1850, 2025):
    src = f"{B}/input/pool_r0103/bc_land_frac_12pfts_{y}.nc"
    with Dataset(src) as d:
        d.set_auto_mask(False)
        f = {i: d[f"fract_pft{i:02d}"][:].astype(float) for i in range(1, 13)}
        land = d["notsea"][:] > 0
    c3m, c4m = f[9] + f[11], f[10] + f[12]
    g = {i: f[i].copy() for i in range(1, 9)}
    g[7] += c3m
    g[8] += c4m
    s8 = sum(g.values())
    # 9 PFTs
    ym = min(max(y, 1984), 2023)
    with Dataset(f"{B}/input/mangrove_13pfts/bc_land_frac_13pfts_{ym}.nc") as m:
        m.set_auto_mask(False); f13 = np.clip(m["fract_pft13"][:].astype(float), 0, 1)
    h = {i: g[i].copy() for i in range(1, 9)}
    p12 = h[1] + h[2]
    take12 = np.minimum(f13, p12)
    for i in (1, 2):
        h[i] -= np.where(p12 > 0, take12 * g[i] / np.where(p12 > 0, p12, 1), 0.0)
    rest = f13 - take12
    tot = sum(h.values())
    for i in range(1, 9):
        h[i] -= np.where(tot > 0, rest * h[i] / np.where(tot > 0, tot, 1), 0.0)
    for i in range(1, 9):
        h[i] = np.maximum(h[i], 0.0)          # rounding of the subtraction can leave -1e-19
    h13 = f13
    s9 = sum(h.values()) + h13
    e8 = np.abs(s8 - 1)[land].max(); e9 = np.abs(s9 - 1)[land].max()
    neg = min(min(v[land].min() for v in g.values()), min(v[land].min() for v in h.values()))
    if e8 > 1e-9 or e9 > 1e-9 or neg < -1e-12: bad += 1
    write(src, f"{B}/output/8pfts/bc_land_frac_8pfts_{y}.nc", {f"fract_pft{i:02d}": g[i] for i in range(1, 9)}, {f"fract_pft{i:02d}" for i in range(1, 9)})
    F9 = {f"fract_pft{i:02d}": h[i] for i in range(1, 9)}; F9["fract_pft13"] = h13
    write(src, f"{B}/output/9pfts/bc_land_frac_9pfts_{y}.nc", F9, {f"fract_pft{i:02d}" for i in range(1, 9)} | {"fract_pft13"})
    if y in (1850, 1900, 1950, 2000, 2020, 2024):
        print(f"{y}: managed share moved {np.average((c3m+c4m)[land]):.3f} (cell mean); max |sum-1| 8 PFTs {e8:.1e}, 9 PFTs {e9:.1e}; smallest fraction {neg:.1e}; mangrove cells {(h13[land]>0).sum()}")
print("years failing a check:", bad)
