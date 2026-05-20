from pathlib import Path
import rasterio

base = Path("/work/bm1548/u241304/preprocessing/pft13_traits/raw_climate")

for pattern in ["*bio01*.tif", "*bio12*.tif"]:
    path = next(base.glob(pattern))
    with rasterio.open(path) as src:
        print("\nFILE:", path.name)
        print("CRS:", src.crs)
        print("Width, height:", src.width, src.height)
        print("Bounds:", src.bounds)
        print("Resolution:", src.res)
        print("NoData:", src.nodata)
        print("Dtype:", src.dtypes[0])
