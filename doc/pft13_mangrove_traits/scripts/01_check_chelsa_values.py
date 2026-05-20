from pathlib import Path
import rasterio
import numpy as np

base = Path("/work/bm1548/u241304/preprocessing/pft13_traits/raw_climate")

files = {
    "bio01_temperature": next(base.glob("*bio01*.tif")),
    "bio12_precipitation": next(base.glob("*bio12*.tif")),
}

for name, path in files.items():
    with rasterio.open(path) as src:
        window = rasterio.windows.Window(10000, 5000, 2000, 2000)
        data = src.read(1, window=window, masked=True)
        data = data.compressed()

        print("\nFILE:", path.name)
        print("VARIABLE:", name)
        print("MIN raw:", np.nanmin(data))
        print("MEAN raw:", np.nanmean(data))
        print("MAX raw:", np.nanmax(data))

        if "bio01" in path.name:
            temp_c = data / 10.0 - 273.15
            print("MIN temp degC:", np.nanmin(temp_c))
            print("MEAN temp degC:", np.nanmean(temp_c))
            print("MAX temp degC:", np.nanmax(temp_c))

        if "bio12" in path.name:
            print("BIO12 is annual precipitation in mm/year")
