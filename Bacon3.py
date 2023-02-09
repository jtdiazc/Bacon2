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
from statsmodels.distributions.empirical_distribution import ECDF

#Version of simulation
sim_vers="20221122"

#We set directory
wdir=os.path.join(r"C:\Projects\5630",sim_vers)

np_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\Numpy",sim_vers)

shp_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector",sim_vers)

ras_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster",sim_vers)

csv_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV",sim_vers)

SLR=pd.read_csv(r"C:\Projects\5630\20221122\SLR.csv",index_col=0)


ml = flopy.modflow.Modflow.load(os.path.join(wdir,r'Base\BAU\MF_inputs\Bacon.nam'))

toedrains_in=pd.read_csv(os.path.join(wdir,r"ToeDrains_Index_Flopy.csv"))

Transects=pd.read_csv(os.path.join(wdir,"Transects.csv"))

rice_df=pd.read_csv(os.path.join(wdir,"Base","WLR","Rice.csv"))

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
#sens="UB"
#isam.run_WLR(sens,Start_Year,End_Year)
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

#Let's plot histogram of peat bottom elevations for 2018

Base_2018_BAU_path=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV\20221122\Base_BAU_RPF_2018.csv"
RPF_2018_df=pd.read_csv(Base_2018_BAU_path,index_col=0)
out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2023_01\RPF"



#Let's plot histogram of marsh deposit elevations for 2018
pyhf.utils.histogram(RPF_2018_df,
                     "PT_bot",
                     20,
                     out_path,
                     "BAU_2018_Marsh_bot_hist",
                     "Marsh Deposits Bottom (ft)",
                     "Density")

#Let's plot histogram of  RPF due to seepage  for 2018
pyhf.utils.histogram(RPF_2018_df,
                     "RPF_Seep",
                     20,
                     out_path,
                     "BAU_2018_RPF_Seep_hist",
                     "RPF Seepage",
                     "Density")

#Let's plot histogram of  slope RPF for 2018
pyhf.utils.histogram(RPF_2018_df,
                     "RPF_Slope",
                     20,
                     out_path,
                     "BAU_2018_RPF_Slope_hist",
                     "RPF Slope",
                     "Density")

#Let's plot histogram of total RPF for 2018
pyhf.utils.histogram(RPF_2018_df,
                     "RPF_Total",
                     20,
                     out_path,
                     "BAU_2018_RPF_Total_hist",
                     "RPF total",
                     "Density")

#Let's plot ECDF of total gradients for various years of the BAU scenario

years=[2018,2030,2040,2050,2060,2070]

# We initiate figure
fig, ax = plt.subplots()
fig.set_figwidth(7)
fig.set_figheight(7)
nbins=1000
min=0
max=1
x=np.linspace(min,max,nbins)
line_width=1
leg_x=1
leg_y=0.8

#Due to seepage
#Let's loop through years
for year in years:
    #Let's import dataframe
    df_dum=pd.read_csv(os.path.join(csv_dir,"Base_BAU_RPF_"+str(year)+".csv"))
    #Let's calculate ECDF
    ecdf_dum=ECDF(df_dum['RPF_Seep'])
    ax.plot(x, ecdf_dum(x), label=str(year), linewidth=line_width)

ax.legend(bbox_to_anchor=(leg_x, leg_y))
ax.set_xlabel('RPF Seepage')
plt.savefig(os.path.join(out_path, "ecdf_Base_BAU_Seepage.png"))
plt.close()

#Due to slope

fig, ax = plt.subplots()
fig.set_figwidth(7)
fig.set_figheight(7)

for year in years:
    #Let's import dataframe
    df_dum=pd.read_csv(os.path.join(csv_dir,"Base_BAU_RPF_"+str(year)+".csv"))
    #Let's calculate ECDF
    ecdf_dum=ECDF(df_dum['RPF_Slope'])
    ax.plot(x, ecdf_dum(x), label=str(year), linewidth=line_width)

ax.legend(bbox_to_anchor=(leg_x, leg_y))
ax.set_xlabel('RPF Slope')
plt.savefig(os.path.join(out_path, "ecdf_Base_BAU_Slope.png"))
plt.close()

#Total

fig, ax = plt.subplots()
fig.set_figwidth(7)
fig.set_figheight(7)

for year in years:
    #Let's import dataframe
    df_dum=pd.read_csv(os.path.join(csv_dir,"Base_BAU_RPF_"+str(year)+".csv"))
    #Let's calculate ECDF
    ecdf_dum=ECDF(df_dum['RPF_Total'])
    ax.plot(x, ecdf_dum(x), label=str(year), linewidth=line_width)

ax.legend(bbox_to_anchor=(leg_x, leg_y))
ax.set_xlabel('RPF Total')
plt.savefig(os.path.join(out_path, "ecdf_Base_BAU_Total.png"))
plt.close()

##Let's calculate percentiles for 2070
RPF_2070_df=pd.read_csv(os.path.join(csv_dir,"Base_BAU_RPF_"+str(2070)+".csv"))

##Due to seepage
perc_seep_2070=np.percentile(RPF_2070_df.RPF_Seep,[25,50,75])

##Due to static slope failure
perc_slope_2070=np.percentile(RPF_2070_df.RPF_Slope,[25,50,75])

##Total RPFs
perc_slope_2070=np.percentile(RPF_2070_df.RPF_Total,[25,50,75])


#Let's do WLR scenario now
#Let's plot ECDF of total gradients for various years of the WLR scenario

years=[2018,2030,2040,2050,2060,2070]

# We initiate figure
fig, ax = plt.subplots()
fig.set_figwidth(7)
fig.set_figheight(7)
nbins=1000
min=0
max=1
x=np.linspace(min,max,nbins)
line_width=1
leg_x=1
leg_y=0.8

#Due to seepage
#Let's loop through years
for year in years:
    #Let's import dataframe
    df_dum=pd.read_csv(os.path.join(csv_dir,"Base_WLR_RPF_"+str(year)+".csv"))
    #Let's calculate ECDF
    ecdf_dum=ECDF(df_dum['RPF_Seep'])
    ax.plot(x, ecdf_dum(x), label=str(year), linewidth=line_width)

ax.legend(bbox_to_anchor=(leg_x, leg_y))
ax.set_xlabel('RPF Seepage')
plt.savefig(os.path.join(out_path, "ecdf_Base_WLR_Seepage.png"))
plt.close()

#Due to slope

fig, ax = plt.subplots()
fig.set_figwidth(7)
fig.set_figheight(7)

for year in years:
    #Let's import dataframe
    df_dum=pd.read_csv(os.path.join(csv_dir,"Base_WLR_RPF_"+str(year)+".csv"))
    #Let's calculate ECDF
    ecdf_dum=ECDF(df_dum['RPF_Slope'])
    ax.plot(x, ecdf_dum(x), label=str(year), linewidth=line_width)

ax.legend(bbox_to_anchor=(leg_x, leg_y))
ax.set_xlabel('RPF Slope')
plt.savefig(os.path.join(out_path, "ecdf_Base_WLR_Slope.png"))
plt.close()

#Total

fig, ax = plt.subplots()
fig.set_figwidth(7)
fig.set_figheight(7)

for year in years:
    #Let's import dataframe
    df_dum=pd.read_csv(os.path.join(csv_dir,"Base_WLR_RPF_"+str(year)+".csv"))
    #Let's calculate ECDF
    ecdf_dum=ECDF(df_dum['RPF_Total'])
    ax.plot(x, ecdf_dum(x), label=str(year), linewidth=line_width)

ax.legend(bbox_to_anchor=(leg_x, leg_y))
ax.set_xlabel('RPF Total')
plt.savefig(os.path.join(out_path, "ecdf_Base_WLR_Total.png"))
plt.close()

##Let's calculate percentiles for 2070
RPF_2070_df_WLR=pd.read_csv(os.path.join(csv_dir,"Base_WLR_RPF_"+str(2070)+".csv"))

##Due to seepage
perc_seep_2070_WLR=np.percentile(RPF_2070_df_WLR.RPF_Seep,[25,50,75])
perc_seep_2018=np.percentile(RPF_2018_df.RPF_Seep,[25,50,75])

##Due to static slope failure
perc_slope_2070_WLR=np.percentile(RPF_2070_df_WLR.RPF_Slope,[25,50,75])
perc_slope_2018=np.percentile(RPF_2018_df.RPF_Slope,[25,50,75])

##Total RPFs
perc_total_2070_WLR=np.percentile(RPF_2070_df_WLR.RPF_Total,[25,50,75])
perc_total_2018=np.percentile(RPF_2018_df.RPF_Total,[25,50,75])