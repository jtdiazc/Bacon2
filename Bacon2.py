import pyhf
import os
import pandas as pd
import flopy
import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
#from moviepy.editor import VideoClipPM
from moviepy.video.io.bindings import mplfig_to_npimage

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

if True:
    isam.init_SUBCALC()

sensitivities=["Base","UB","LB"]

Start_Year=2018
End_Year=2070

run_BAU=False

sens="Base"
#isam.run_WLR(sens,Start_Year,End_Year)
isam.run_BAU(sens,Start_Year,End_Year)
#sens="LB"
#isam.run_WLR(sens,Start_Year,End_Year)
#isam.run_BAU(sens,Start_Year,End_Year)
#sens="UB"
#isam.run_WLR(sens,Start_Year,End_Year)
#isam.run_BAU(sens,Start_Year,End_Year)


#isam.run_BAU(sens,Start_Year,2045)
#isam.run_BAU(sens,2045,2045,resume_run=True)
#isam.run_BAU(sens,2046,End_Year,resume_run=True)








#for sens in sensitivities:



#    if run_BAU:
#        isam.run_BAU(sens,Start_Year,End_Year)
        #isam.run_BAU(sens, 2050, End_Year,resume_run=True)

        # Let's create band plots
#        pyhf.utils.band_plot(Start_Year-1, End_Year, isam.RPF_BAU_out,leg_y=0.2,leg_x=1)

        # Let's export shapefile of cross sections for first and last time step
#        pyhf.utils.shp_df_join(isam.RPF_BAU_out, ["Transects_" + str(Start_Year - 1) + ".csv",
#                                                  "Transects_" + str(End_Year) + ".csv"]

#                               , os.path.join(shp_dir, "RPF", "Transect_Cells_v3.shp"),
#                               os.path.join(shp_dir, "RPF"),
#                               sens + "_BAU_"
#                               )

#    isam.run_WLR(sens,Start_Year,End_Year)

    # Let's create band plots
#    pyhf.utils.band_plot(Start_Year-1, End_Year, isam.RPF_WLR_out)



#    pyhf.utils.shp_df_join(isam.RPF_WLR_out,["Transects_"+str(Start_Year-1)+".csv",
#                                        "Transects_" + str(End_Year) + ".csv"]
#                           ,os.path.join(shp_dir,"RPF","Transect_Cells_v3.shp"),
#                           os.path.join(shp_dir, "RPF"),
#                           sens + "_WLR_"
#                           )

###Let's plot an histogram of the peat bottom elevations for the Base scenario, BAU and 2018 to assess the baseline
###RPF


#Base_2018_BAU_dir=os.path.join(wdir,"Base", "BAU", "Output", "RPF")
#RPF_2018_df=pd.read_csv(os.path.join(Base_2018_BAU_dir,"Transects_2017.csv"),index_col=0)
#out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_11\RPF"

#pyhf.utils.histogram(RPF_2018_df,
#                     "PT_bot",
#                     20,
#                     out_path,
#                     "BAU_2018_PT_bot_hist",
#                     "Peat Bottom (ft)",
#                     "Density")

#Now, let's plot an histogram of the surface elevations at levee toes
#pyhf.utils.histogram(RPF_2018_df,
#                     'PT_top',
#                     20,
#                     out_path,
#                     "BAU_2018_PT_top_hist",
#                     "Peat Top (ft)",
#                     "Density")

#Now, let's plot an histogram of levee crests
#pyhf.utils.histogram(RPF_2018_df,
#                     'Z',
#                     20,
#                     out_path,
#                     "BAU_2018_Levee_top_hist",
#                     "Levee Crest Elevation (ft)",
#                     "Density")

#Now, let's plot a histogram of RPFseep
pyhf.utils.histogram(RPF_2018_df,
                     'RPF_Seep',
                     20,
                     out_path,
                     "BAU_2018_RPF_Seep_hist",
                     "RPF seepage",
                     "Density")

#Now, let's plot a histogram of RPFslope
pyhf.utils.histogram(RPF_2018_df,
                     'RPF_Slope',
                     20,
                     out_path,
                     "BAU_2018_RPF_Slope_hist",
                     "RPF slope",
                     "Density")

#Now, let's plot a histogram of total RPF

pyhf.utils.histogram(RPF_2018_df,
                     'RPF_Total',
                     20,
                     out_path,
                     "BAU_2018_RPF_Total_hist",
                     "RPF total",
                     "Density")

#Subsidence statistics for Base scenario
ras_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\20221122\Base\BAU"
prefix="Subs_ft_"
years=np.array(range(2018,2071))
BAU_sub_stats=pyhf.utils.ras_stats(ras_dir,prefix,years)

#Subsidence statistics for LB scenario
ras_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\20221122\LB\BAU"
prefix="Subs_ft_"
years=np.array(range(2018,2071))
BAU_LB_sub_stats=pyhf.utils.ras_stats(ras_dir,prefix,years)

#Subsidence statistics for UB scenario
ras_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\20221122\UB\BAU"
prefix="Subs_ft_"
years=np.array(range(2018,2071))
BAU_UB_sub_stats=pyhf.utils.ras_stats(ras_dir,prefix,years)



#Path to csv files
csv_path=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV\20221122"

#Names of series
names=["BAU","WLR"]

out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_12\Gradients"

#min, max=pyhf.utils.density_plot(csv_path,names,out_path,video=True,calc_x=False,kernel='linear',bandwidth=0.1,bands=True,alpha=1.0)

#min2, max2=pyhf.utils.ecdf(csv_path,names,out_path,video=True,End_Year=End_Year,bands=True)

#Density plots for RPF

out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_12\RPF"

#Let's move and rename RPF files
#Base

