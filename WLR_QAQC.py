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

isam.sens="Base"

isam.sim_years = isam.End_Year - isam.Start_Year + 2

if isam.sens == "Base":
    isam.sedcalc_dum = pyhf.SEDCALC.sedcalc(h2oin=pd.read_csv(os.path.join("Base", "WLR", "h2oin.csv")),
                                            minin=pd.read_csv(os.path.join("Base", "WLR", "minin.csv")),
                                            orgin=pd.read_csv(os.path.join("Base", "WLR", "orgin.csv")),
                                            porelim=pd.read_csv(os.path.join("Base", "WLR", "porelim.csv")),
                                            endtim=isam.sim_years)
elif isam.sens == "LB":
    isam.sedcalc_dum = pyhf.SEDCALC.sedcalc(h2oin=pd.read_csv(os.path.join("Base", "WLR", "h2oin.csv")),
                                            minin=pd.read_csv(os.path.join("Base", "WLR", "minin.csv")),
                                            orgin=pd.read_csv(os.path.join("Base", "WLR", "orgin.csv")),
                                            porelim=pd.read_csv(os.path.join("Base", "WLR", "porelim_LB.csv")),
                                            endtim=isam.sim_years,
                                            k1=0.5)
elif isam.sens == "UB":
    isam.sedcalc_dum = pyhf.SEDCALC.sedcalc(h2oin=pd.read_csv(os.path.join("Base", "WLR", "h2oin.csv")),
                                            minin=pd.read_csv(os.path.join("Base", "WLR", "minin.csv")),
                                            orgin=pd.read_csv(os.path.join("Base", "WLR", "orgin.csv")),
                                            porelim=pd.read_csv(os.path.join("Base", "WLR", "porelim_UB.csv")),
                                            endtim=isam.sim_years,
                                            k1=12.5)

isam.relelv_dum = isam.sedcalc_dum["relelv"].values
isam.depth_dum = isam.sedcalc_dum["depth"].values
isam.porg_dum = isam.sedcalc_dum["porg"].values
isam.bulkd_dum = isam.sedcalc_dum["bulkd"].values
isam.acc_rate_dum = pd.DataFrame(columns=["Year", "Yearly Accretion (cm)", "Yearly Accretion (ft)"])
isam.totC_dum = pd.DataFrame(columns=["Year", "Yearly Accretion (cm)", "gC/cm3", "gC/cm2"])
year = isam.Start_Year

for i in range(isam.sim_years - 1):
    # Index for the backwards sum
    i_bw = isam.sim_years - i - 1

    #        else:
    acc_rate_dum_dum = pd.DataFrame({"Year": [year],
                                     "Yearly Accretion (cm)": [(isam.depth_dum[i_bw] - isam.depth_dum[i_bw - 1])],
                                     "Yearly Accretion (ft)": [
                                         (isam.depth_dum[i_bw] - isam.depth_dum[i_bw - 1]) * 0.0328084]})
    totC_dum_dum = pd.DataFrame({"Year": [year],
                                 "Yearly Accretion (cm)": [(isam.depth_dum[i_bw] - isam.depth_dum[i_bw - 1])],
                                 "gC/cm3": [isam.porg_dum[i_bw] * isam.bulkd_dum[i_bw] / 2],
                                 "gC/cm2": [
                                     (isam.depth_dum[i_bw] - isam.depth_dum[i_bw - 1]) * isam.porg_dum[i_bw] *
                                     isam.bulkd_dum[
                                         i_bw] / 2]})

    isam.acc_rate_dum = isam.acc_rate_dum.append(acc_rate_dum_dum, ignore_index=True)
    isam.totC_dum = isam.totC_dum.append(totC_dum_dum, ignore_index=True)
    year = year + 1


isam.acc_rate_dum["Yearly Accretion (ft)"] = isam.acc_rate_dum["Yearly Accretion (cm)"]*0.0328084
# Let's export the total carbon table as a csv file
isam.totC_dum.to_csv(os.path.join(isam.sens, "WLR", "SEDCALC_totC_output.csv"), index=False)

# Let's export the SEDCALC accretion rates
isam.acc_rate_dum.to_csv(os.path.join(isam.sens, "SEDCALC_Accretion.csv"), index=False)
isam.sedcalc_ts = isam.acc_rate_dum

isam.WLR_ras_dir = os.path.join(isam.ras_dir, isam.sens, "WLR")
if not os.path.exists(isam.WLR_ras_dir):
    os.makedirs(isam.WLR_ras_dir)
# Directory for shapefiles
isam.shp_WLR_path = os.path.join(isam.shp_dir, isam.sens, "WLR")
if not os.path.exists(isam.shp_WLR_path):
    os.makedirs(isam.shp_WLR_path)

isam.RPF_WLR_out = os.path.join(isam.wdir, isam.sens, "WLR", "Output", "RPF")
if not os.path.exists(isam.RPF_WLR_out):
    os.makedirs(isam.RPF_WLR_out)


#Let's loop through years now
for year in range(2018, 2018):
    if year == isam.Start_Year:
        isam.ml = flopy.modflow.Modflow.load('Base/WLR/MF_inputs/Bacon.nam')
        isam.ml.change_model_ws(new_pth=os.path.join(isam.sens, "WLR", str(year)))
        isam.ml.pcg.rclose = 864
        isam.ml.exe_name = isam.mf_exe_path
        isam.ml.modelgrid.set_coord_info(xoff=isam.xoff, yoff=isam.yoff, epsg=isam.epsg)
        # Let's set rice and wetland to constant head
        isam.ml.bas6.ibound[0][isam.rice_mask] = -1
        isam.ml.bas6.ibound[0][isam.wetland_mask] = -1
        # Let's add constant heads for rice
        isam.ml.bas6.strt[0][isam.rice_mask] = isam.ml.dis.top[isam.rice_mask]

        isam.ml.write_input()
        isam.ml.run_model()
        pyhf.RPF.rpf(isam.Transects, isam.ml, year - 1, isam.SLR, isam.RPF_WLR_out)

    # Let´s add wetlands accretion to land surface
    isam.ml.dis.top[isam.wetland_mask] = np.minimum(isam.ml.dis.top[isam.wetland_mask] +
                                                    float(isam.sedcalc_ts.loc[isam.sedcalc_ts.Year == year,
                                                    "Yearly Accretion (ft)"]),isam.SLR.loc[year, "2_ft"])
    # Let´s add wetlands constant head values
    isam.ml.bas6.strt[0][isam.wetland_mask] = np.minimum(isam.ml.dis.top[isam.wetland_mask] + 1,
                                                         isam.SLR.loc[year, "2_ft"])
    # Peat thickness
    isam.PT_thck = isam.ml.dis.gettop()[0] - isam.ml.dis.getbotm()[0]

    # Depth to groundwater
    isam.get_DTW()

    # We sample elevations
    isam.drns['top'] = isam.ml.dis.top[isam.drns['i'], isam.drns['j']]

    # Let's update constant heads

    for layer in range(3):
        isam.ml.bas6.strt[layer][isam.CH] = isam.SLR.loc[year, "2_ft"]

    isam.ml.change_model_ws(new_pth=os.path.join(os.path.join(isam.sens, "WLR", str(year))))
    isam.ml.write_input()
    # Let's run MODFLOW
    isam.ml.run_model()

    isam.h = flopy.utils.HeadFile(os.path.join(isam.ml.model_ws, "Bacon.hds"), model=isam.ml)
    isam.heads = isam.h.get_data()
    isam.drns['h_PT'] = isam.heads[0][isam.drns['i'], isam.drns['j']]
    isam.drns['h_TM'] = isam.heads[1][isam.drns['i'], isam.drns['j']]
    isam.drns['h_SP'] = isam.heads[2][isam.drns['i'], isam.drns['j']]

    # Let's calculate gradients
    isam.calc_grads()

    flopy.export.shapefile_utils.recarray2shp(isam.drns,
                                              geoms=isam.polygons,
                                              shpname=os.path.join(isam.shp_WLR_path,
                                                                   "WLR_DRN_" + str(year) + ".shp"),
                                              epsg=isam.grid.epsg)
    # Let's export heads
    flopy.export.utils.export_array(isam.grid, os.path.join(isam.WLR_ras_dir, "H_Pt_ft_" + str(year) + ".tif"), isam.heads[0])
    flopy.export.utils.export_array(isam.grid, os.path.join(isam.WLR_ras_dir, "H_TM_ft_" + str(year) + ".tif"), isam.heads[1])
    flopy.export.utils.export_array(isam.grid, os.path.join(isam.WLR_ras_dir, "H_SP_ft_" + str(year) + ".tif"), isam.heads[2])
    isam.drns_pd = pd.DataFrame(isam.drns)

    # Let's subset to toedrains
    isam.toedrains_dum = pd.merge(isam.drns_pd, isam.toedrains_in, how="right", on=["i", "j"])
    # Let's convert to recarray
    isam.toedrains_dum_rec = isam.toedrains_dum.to_records()

    flopy.export.shapefile_utils.recarray2shp(isam.toedrains_dum_rec,
                                              geoms=isam.polygons_toedrn,
                                              shpname=os.path.join(isam.shp_WLR_path,
                                                                   "WLR_TOEDRNS_" + str(year) + ".shp"),
                                              epsg=isam.grid.epsg)
    isam.drns_pd["Year"] = year
    isam.toedrains_dum["Year"] = year
    isam.drns_pd.to_csv(os.path.join(isam.csv_dir,isam.sens + "_WLR_DRNS" + str(year) + ".csv"), index=False)
    isam.toedrains_dum.to_csv(os.path.join(isam.csv_dir, isam.sens +"_WLR_TOEDRNS" + str(year) + ".csv"), index=False)

    # Let's export top elevation
    flopy.export.utils.export_array(isam.grid, os.path.join(isam.WLR_ras_dir, "Top_ft_" + str(year) + ".tif"),
                                    isam.ml.dis.top[:])
    pyhf.RPF.rpf(isam.Transects, isam.ml, year, isam.SLR, isam.RPF_WLR_out)



isam.rice_df
isam.rice_mask = list(zip(isam.rice_df.row.values, isam.rice_df.column_lef.values))
isam.rice_mask = tuple(np.array(isam.rice_mask).T)

# Let´s create mask of wetlands
isam.wetland_df = pd.merge(isam.ACNL_df, isam.rice_df, how="left", left_on=["row", "col"],
                           right_on=["row", "column_lef"],
                           indicator=True)

isam.wetland_df = isam.wetland_df[isam.wetland_df._merge == "left_only"].reset_index()

isam.wetland_mask = list(zip(isam.wetland_df.row.values, isam.wetland_df.col.values))