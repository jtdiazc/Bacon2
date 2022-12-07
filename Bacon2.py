import pyhf
import os
import pandas as pd
import flopy

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

sensitivities=["Base","UB","LB"]

Start_Year=2018
End_Year=2019

for sens in sensitivities:




    isam.run_BAU(sens,Start_Year,End_Year)

    # Let's create band plots
    pyhf.utils.band_plot(Start_Year-1, End_Year, isam.RPF_BAU_out,leg_y=0.2,leg_x=1)

    isam.run_WLR(sens,Start_Year,End_Year)

    # Let's create band plots
    pyhf.utils.band_plot(Start_Year-1, End_Year, isam.RPF_WLR_out)

    #Let's export shapefile of cross sections for first and last time step
    pyhf.utils.shp_df_join(isam.RPF_BAU_out,["Transects_"+str(Start_Year-1)+".csv",
                                    "Transects_" + str(End_Year) + ".csv"]
                       ,os.path.join(shp_dir,"RPF","Transect_Cells_v3.shp"),
                       os.path.join(shp_dir, "RPF"),
                       sens + "_BAU_"
                       )

    pyhf.utils.shp_df_join(isam.RPF_WLR_out,["Transects_"+str(Start_Year-1)+".csv",
                                        "Transects_" + str(End_Year) + ".csv"]
                           ,os.path.join(shp_dir,"RPF","Transect_Cells_v3.shp"),
                           os.path.join(shp_dir, "RPF"),
                           sens + "_WLR_"
                           )

###Let's plot an histogram of the peat bottom elevations for the Base scenario, BAU and 2018 to assess the baseline
###RPF


Base_2018_BAU_dir=os.path.join(wdir,"Base", "BAU", "Output", "RPF")
RPF_2018_df=pd.read_csv(os.path.join(Base_2018_BAU_dir,"Transects_2017.csv"),index_col=0)
out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_11\RPF"

pyhf.utils.histogram(RPF_2018_df,
                     "PT_bot",
                     20,
                     out_path,
                     "BAU_2018_PT_bot_hist",
                     "Peat Bottom (ft)",
                     "Density")

#Now, let's plot an histogram of the surface elevations at levee toes
pyhf.utils.histogram(RPF_2018_df,
                     'PT_top',
                     20,
                     out_path,
                     "BAU_2018_PT_top_hist",
                     "Peat Top (ft)",
                     "Density")

#Now, let's plot an histogram of levee crests
pyhf.utils.histogram(RPF_2018_df,
                     'Z',
                     20,
                     out_path,
                     "BAU_2018_Levee_top_hist",
                     "Levee Crest Elevation (ft)",
                     "Density")

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
