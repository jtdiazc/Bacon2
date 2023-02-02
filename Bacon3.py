import pyhf
import os
import pandas as pd
import flopy
import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
#from moviepy.editor import VideoClipPM
from moviepy.video.io.bindings import mplfig_to_npimage
import shutil

#Version of simulation
sim_vers="20221122"

#We set directory
wdir=os.path.join(r"C:\Projects\5630",sim_vers)

np_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\Numpy",sim_vers)

shp_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector",sim_vers)

ras_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster",sim_vers)

csv_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV",sim_vers)

SLR=pd.read_csv(r"C:\Projects\5630\20221122\SLR.csv",index_col=0)


ml = flopy.modflow.Modflow.load('Base/BAU/MF_inputs/Bacon.nam')

toedrains_in=pd.read_csv(r"ToeDrains_Index_Flopy.csv")

Transects=pd.read_csv("Transects.csv")

rice_df=pd.read_csv(os.path.join("Base","WLR","Rice.csv"))

isam=pyhf.isa(wdir,np_dir,shp_dir,ras_dir,csv_dir,SLR,ml,toedrains_in, Transects, rice_df)

#Organic matter to raster
if False:
    isam.init_SUBCALC()
    flopy.export.utils.export_array(isam.grid, os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\20221122\Base\BAU", "fom_0.tif"),
                                    isam.SC_Input["fom"])
    flopy.export.utils.export_array(isam.grid, os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\20221122\Base\BAU", "fomund_0.tif"),
                                    isam.SC_Input["fomund"])
    isam.run_BAU("Base",2018,2018)


sensitivities=["Base","UB","LB"]

Start_Year=2018
End_Year=2070



#sens="Base"
#isam.run_WLR(sens,Start_Year,End_Year)
#isam.run_BAU(sens,Start_Year,End_Year)
#sens="LB"
#isam.run_WLR(sens,Start_Year,End_Year)
#isam.run_BAU(sens,Start_Year,End_Year)
sens="UB"
isam.run_WLR(sens,Start_Year,End_Year)
#isam.run_BAU(sens,Start_Year,End_Year)

#Density plots for RPF
out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_12\RPF"
csv_path=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV\20221122"



#Let's move and rename RPF files
#Base
exp_prefix="Base"
suffix='RPF_'
sensitivities=["Base","LB","UB"]
names=["BAU","WLR"]
for sens in sensitivities:
    for run in names:

        #Let's loop through years
        for year in range(Start_Year,End_Year+1):
            source=os.path.join(wdir,sens,run,"Output","RPF","Transects_"+str(year)+".csv")
            target=os.path.join(csv_path,sens+"_"+run+"_"+suffix+str(year)+".csv")
            shutil.copy(source, target)
wdir=r"C:\Projects\5630\20221122"
sens="Base"
year=2018
run="BAU"

min, max=pyhf.utils.density_plot(csv_path,names,out_path,video=True,
                                 calc_x=False,
                                 kernel='linear',
                                 bandwidth=0.1,
                                 bands=True,
                                 alpha=0.2,
                                 min_x=0,
                                 valkey="RPF_Total",
                                 max_x=1,
                                 suffix=suffix,
                                 End_Year=End_Year)

min2, max2=pyhf.utils.ecdf(csv_path,names,out_path,video=True,
                                 calc_x=False,
                                 kernel='linear',
                                 bandwidth=0.1,
                                 bands=True,
                                 alpha=0.2,
                                 min_x=0,
                                 valkey="RPF_Total",
                                 max_x=1,
                                 suffix=suffix,
                                 End_Year=End_Year)

out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_12\Gradients"
min, max=pyhf.utils.density_plot(csv_path,names,out_path,video=True,calc_x=False,kernel='linear',bandwidth=0.1,bands=True,alpha=1.0)

min2, max2=pyhf.utils.ecdf(csv_path,names,out_path,video=True,End_Year=End_Year,bands=True)




