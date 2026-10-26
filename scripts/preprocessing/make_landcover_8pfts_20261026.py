#!/usr/bin/env python3
# =============================================================================
#  MAKE_LANDCOVER_8PFTS  --  land-cover files without pastures and crops for the dedicated 8-PFT use case
#  Dr Sri Nandini-Weiss, ZMT Bremen.                                                      2026-10-26
#  RUN     python3 make_landcover_8pfts_20261026.py <outdir>          (Levante, login node is enough: 175 small files)
#  Input   /pool/data/ICON/grids/public/mpim/0049/land/r0103/bc_land_frac_12pfts_YYYY.nc      YYYY = 1850 ... 2024
#  Rule    PFT01-PFT06 unchanged. PFT07 = f7 + f9 + f11 (C3 pasture + C3 crop to C3 grass). PFT08 = f8 + f10 + f12 (C4 pasture + C4 crop to C4 grass).
#          notsea, fract_veg, fract_land, veg_ratio_max and every other variable are copied. Only fract_pft01 ... fract_pft08 are kept.
#  Output  <outdir>/bc_land_frac_8pfts_YYYY.nc
#  Same numbers as the 8-PFT part of make_landcover_8_9pfts_20261023.py (Mac, inputs copied to input/), which is no longer needed.
# =============================================================================
import os, sys
import numpy as np
from netCDF4 import Dataset
POOL = "/pool/data/ICON/grids/public/mpim/0049/land/r0103"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
KEEP = {f"fract_pft{i:02d}" for i in range(1, 9)}
bad = 0
for y in range(1850, 2025):
    src = f"{POOL}/bc_land_frac_12pfts_{y}.nc"
    with Dataset(src) as s:
        s.set_auto_mask(False)
        f = {i: np.array(s[f"fract_pft{i:02d}"][:], dtype=float) for i in range(1, 13)}
        land = np.array(s["notsea"][:]) > 0
        g = {i: f[i].copy() for i in range(1, 9)}
        g[7] += f[9] + f[11]
        g[8] += f[10] + f[12]
        err = np.abs(sum(g.values()) - 1)[land].max()
        neg = min(v[land].min() for v in g.values())
        if err > 1e-9 or neg < -1e-12:
            bad += 1
        with Dataset(f"{OUT}/bc_land_frac_8pfts_{y}.nc", "w", format=s.data_model) as o:
            o.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
            o.history = "2026-10-26 make_landcover_8pfts_20261026.py: pastures and crops removed (area to the grass of the pathway); " + getattr(s, "history", "")
            for n, dim in s.dimensions.items():
                o.createDimension(n, len(dim))
            for n, v in s.variables.items():
                if n.startswith("fract_pft") and n not in KEEP:
                    continue
                w = o.createVariable(n, v.dtype, v.dimensions, fill_value=getattr(v, "_FillValue", None))
                w.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"})
                w[:] = g[int(n[-2:])] if n in KEEP else v[:]
    if y in (1850, 1979, 2020, 2024):
        print(f"{y}: max |sum-1| {err:.1e}; smallest fraction {neg:.1e}")
print("years failing a check:", bad)
