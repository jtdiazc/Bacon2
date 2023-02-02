import rasterio
import pandas as pd
import os
import numpy as np

fom_0 = rasterio.open(r'\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\OM_QAQC\fom_0_v2.tif')
fomund_0 = rasterio.open(r'\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\OM_QAQC\fomund_0_v2.tif')

fom_0_np= fom_0.read(1)
fomund_0_np= fomund_0.read(1)

np.save(r"C:\Projects\5630\20221122\Base\fom_0.npy",fom_0_np)
np.save(r"C:\Projects\5630\20221122\Base\fomund_0.npy",fomund_0_np)