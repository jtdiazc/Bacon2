import sys
#path to the hydrofocus functions module
#sys.path.insert(0, r'\\hydro-nas\Team\Projects\5630_DSC\Codes')
import pyhf
#Path to flopy module
#sys.path.insert(0, r'\\hydro-nas\Team\Projects\5630_DSC\Codes\flopy')
import flopy
import os
import pandas as pd
import numpy as np
import rasterio
import subprocess
import datetime
import geopandas as gpd

#Version of simulation
sim_vers="20221122"

#We set directory
wdir=os.path.join(r"C:\Projects\5630",sim_vers)

if not os.path.exists(wdir):
    os.makedirs(wdir)
os.chdir(wdir)

#Numpy arrays
np_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\Numpy",sim_vers)
if not os.path.exists(np_dir):
    os.makedirs(np_dir)


#Shapefiles
shp_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector",sim_vers)
if not os.path.exists(shp_dir):
    os.makedirs(shp_dir)

#Rasters
ras_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster",sim_vers)
if not os.path.exists(ras_dir):
    os.makedirs(ras_dir)

#CSVs
csv_dir=os.path.join(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV",sim_vers)
if not os.path.exists(csv_dir):
    os.makedirs(csv_dir)


#MODFLOW executable path
mf_exe_path=r"C:\wrdapp\MF2005.1_12\bin\mf2005.exe"


##Sea Level Rise timeseries
SLR=pd.read_csv("SLR.csv",index_col=0)

#First year of simulation
Start_Year=2018

#Last year of simulation
End_Year=2070

#Let's load model
ml = flopy.modflow.Modflow.load('Base/BAU/MF_inputs/Bacon.nam')



ml.modelgrid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)



nrg=ml.nrow
ncg=ml.ncol

#Let's get active cells
#bas = flopy.modflow.ModflowBas.load('MF_inputs/Bacon.bas', ml)
bas=ml.bas6

ibound=np.array(bas.ibound[0][:])



Active_cells=np.where(ibound!=0)
Inactive_cells=np.where(ibound==0)

grid = ml.modelgrid

subsidence=np.zeros(shape=(nrg,ncg))

#Let's get mask of levees
levees=np.where(ml.lpf.hk[0][:]==min(np.unique(ml.lpf.hk[0][:])))

#Let's get indexes of toedrains
toedrains_in=pd.read_csv(r"ToeDrains_Index_Flopy.csv")
toedrains_in=toedrains_in.drop(['k','elev', 'cond', 'iface', 'top', 'PT_bot', 'TM_bot','h_PT', 'h_TM', 'h_SP', 'Grad'], axis=1)



#Let's convert to recarray
#toedrains_in_rec=toedrains_in.to_records

#We sample elevations
new_dt=np.dtype(ml.drn.stress_period_data[0].dtype.descr + [('top', '<f4')]+
                [('PT_bot', '<f4')]+
                [('TM_bot', '<f4')]+
                [('h_PT', '<f4')]+
                [('h_TM', '<f4')]+
                [('h_SP', '<f4')]+
                [('Grad', '<f4')])
drns=np.zeros(ml.drn.stress_period_data[0].shape, dtype=new_dt)
drns[['k','i','j','elev','cond','iface']]=ml.drn.stress_period_data[0][['k','i',
                                                                        'j','elev','cond','iface']]
drns['top']=ml.dis.top[drns['i'],drns['j']]
drns['PT_bot']=ml.dis.botm[0][drns['i'],drns['j']]
drns['TM_bot']=ml.dis.botm[1][drns['i'],drns['j']]



#Indices of constant head cells
CH=np.where(ibound<0)

#Let's create drains shapefile template
vertices = []
for row, col in zip(drns['i'], drns['j']):
    vertices.append(grid.get_cell_vertices(row, col))
polygons = [flopy.utils.geometry.Polygon(vrt) for vrt in vertices]

#Let's create toedrains shapefile template
vertices = []
for row, col in zip(toedrains_in['i'], toedrains_in['j']):
    vertices.append(grid.get_cell_vertices(row, col))
polygons_toedrn = [flopy.utils.geometry.Polygon(vrt) for vrt in vertices]

#Let's create constant head cells shapefile
CH_df = pd.DataFrame({"row":np.where(bas.strt[0][:]==np.unique(bas.strt[0][:])[-1])[0],
                   "col":np.where(bas.strt[0][:]==np.unique(bas.strt[0][:])[-1])[1]})

if False:


    pyhf.flopyf.df_to_shp(grid,CH_df,shp_dir,"CH.shp")

#SLR time series
SLR=pd.read_csv("SLR.csv",index_col=0)

#Let's calculate rates
for scenario in SLR.columns[SLR.columns!="Year"].values:
    #Let's add new column
    SLR[scenario+"_rate"]=0
    SLR.loc[ SLR.index.values[0], scenario + "_rate"]=0
    for year in SLR.index.values[1:]:
        SLR.loc[year, scenario + "_rate"]=SLR.loc[year, scenario]-SLR.loc[year-1, scenario]







#2. Import elevations

#Dataframe with average subsidence
avg_sub_df=pd.DataFrame(columns = ['Year', 'Sub_r_ft_yr'])


##Let's export shapefile of levees
Levees_df = pd.DataFrame({"row":levees[0],
                       "col":levees[1]})
pyhf.flopyf.df_to_shp(grid,Levees_df,shp_dir,"Levees.shp")

    ##Let's export shapefile of active cells
Active_cells_df = pd.DataFrame({"row":Active_cells[0],
                       "col":Active_cells[1]})
if False:
    pyhf.flopyf.df_to_shp(grid,Active_cells_df,shp_dir,"Active_Cells.shp")

#Active cells but no levees
#Let's join both dataframes
AC_Lv_left=pd.merge(Active_cells_df,Levees_df,how="left",on=["row","col"],indicator=True)
#Let's remove the ones that are found in both

ACNL_df=AC_Lv_left[~(AC_Lv_left._merge=='both')].reset_index()
if False:
    pyhf.flopyf.df_to_shp(grid,ACNL_df,shp_dir,"Active_Cells_no_levees.shp")

toedrains_df=toedrains_in.copy()
toedrains_df.columns=["row","col"]

ACNL_df=ACNL_df.drop(columns=['_merge'])

#Let's also remove toedrains
ACNLNT_df=pd.merge(ACNL_df,toedrains_in,how="left",left_on=["row","col"],right_on=["i","j"],indicator=True)
ACNLNT_df=ACNLNT_df[ACNLNT_df._merge=="left_only"].reset_index()
if False:
    pyhf.flopyf.df_to_shp(grid,ACNLNT_df,shp_dir,"Active_Cells_no_levees_no_toedrains.shp")

#Import indexes and elevations of transects
Transects=pd.read_csv("Transects.csv")

#Let´s get mask of rice
rice_df=pd.read_csv(os.path.join("Base","WLR","Rice.csv"))
rice_mask=list(zip(rice_df.row,rice_df.column_lef))
rice_mask=tuple(np.array(rice_mask).T)
if False:
    pyhf.flopyf.df_to_shp(grid,rice_df,shp_dir,"Rice.shp",col_key="column_lef")

#Let´s get create mask of wetlands
wetland_df=pd.merge(ACNL_df,rice_df,how="left",left_on=["row","col"],right_on=["row","column_lef"],indicator=True)
wetland_df=wetland_df[wetland_df._merge=="left_only"].reset_index()
if False:
    pyhf.flopyf.df_to_shp(grid,wetland_df,shp_dir,"Wetlands.shp")

wetland_mask=list(zip(wetland_df.row,wetland_df.col))
wetland_mask=tuple(np.array(wetland_mask).T)

# SUBCALC input
SC_Input = {'fom': np.load(os.path.join("Base", 'fom_0.npy')),
            'fomund': np.load(os.path.join("Base", 'fomund_0.npy')),
            'st': np.load(os.path.join("Base", 'st_0.npy')),
            'km': np.load(os.path.join("Base", 'km_0.npy')),
            'firstkm': np.load(os.path.join("Base", 'firstkm_0.npy')),
            'vmax': np.load(os.path.join("Base", 'vmax_0.npy')),
            'bd': np.load(os.path.join("Base", 'bd_0.npy')),
            'bdund': np.load(os.path.join("Base", 'bdund_0.npy')),
            'massmin': np.load(os.path.join("Base", 'massmin_0.npy'))}

# Let's set nas to 0
SC_Input['fom'][np.isnan(SC_Input['fom'])] = 0

#Uncertainty terms (cm/yr)
uncert_SUBCALC={"LB":-0.1, "Base":0,"UB":0.1}
uncert_SEDCALC={"LB":-1.4, "Base":0,"UB":1.4}

if False:
    ml.modelgrid.write_shapefile(os.path.join(shp_dir, "Grid.shp"))


#Now, let's loop through sensitivity scenarios
for sens in ["LB","Base","UB"]:


    ##BAU

    #Location of rasters for BAU
    BAU_ras_dir=os.path.join(ras_dir,sens,"BAU")
    if not os.path.exists(BAU_ras_dir):
        os.makedirs(BAU_ras_dir)

    RPF_BAU_out =os.path.join(wdir,sens, "BAU", "Output", "RPF")
    if not os.path.exists(RPF_BAU_out):
        os.makedirs(RPF_BAU_out)

    for year in range(Start_Year, 2021):
    #3. Let's loop through years now
    #for year in range(Start_Year,End_Year+1):

        # 3.1 Initiate MODFLOW for year 1
        if year==Start_Year:
            ml = flopy.modflow.Modflow.load('Base/BAU/MF_inputs/Bacon.nam')
            ml.change_model_ws(new_pth=os.path.join(sens, "BAU",str(year)))
            ml.pcg.rclose = 864
            ml.exe_name = mf_exe_path
            ml.modelgrid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)
            ml.write_input()
            ml.run_model()
            # Peat thickness
            PT_thck = ml.dis.gettop()[0]-ml.dis.getbotm()[0]
            PT_thck[Inactive_cells] = 0
            PT_thck[levees] = 0
            #Let's calculate starting RPFs
            pyhf.RPF.rpf(Transects,ml,year-1,SLR,RPF_BAU_out)



        #Depth to groundwater
        h = flopy.utils.HeadFile(os.path.join(ml.model_ws,"Bacon.hds"), model=ml)
        heads=h.get_data()

        wt=flopy.utils.postprocessing.get_water_table(heads,masked_values=[999]).data

        DTW = np.maximum(0, (ml.dis.top[:] - wt))


        dtwm = DTW / 3.28084
        peatdepth = PT_thck / 3.28084 * 100



        # Let's export DTW raster
        flopy.export.utils.export_array(grid, os.path.join(BAU_ras_dir, "DTW_m_" + str(year) + ".tif"), dtwm)

        # Let's export Peat thickness cm
        flopy.export.utils.export_array(grid, os.path.join(BAU_ras_dir, "PT_thck_cm_" + str(year) + ".tif"), peatdepth)

        SC_Input=pyhf.subcalc_2021_npy(SC_Input['fom'],
                                  SC_Input['fomund'],
                                  dtwm,
                                  peatdepth,
                                  year,
                                  st=SC_Input['st'],
                                  km=SC_Input['km'],
                                  firstkm=SC_Input['firstkm'],
                                  vmax=SC_Input['vmax'],
                                  bd=SC_Input['bd'],
                                  bdund=SC_Input['bdund'],
                                  massmin=SC_Input['massmin']
                                  )
        #Let's apply uncertainty
        SC_Input['tdelev']=np.maximum(0,SC_Input['tdelev']+uncert_SUBCALC[sens])
        subs=SC_Input['tdelev']*3.28084/100
        #For levees, subsidence is zer0
        subs[levees]=0
        #where subsidence is not na, we update top elevation by subsidence
        ml.dis.top[~np.isnan(subs)]=np.maximum(ml.dis.top[~np.isnan(subs)]-subs[~np.isnan(subs)],ml.dis.botm[0][~np.isnan(subs)])
        PT_thck=np.maximum(PT_thck-subs,0)

        #average subsidence
        subs_avg=np.average(subs[subs>0])
        avg_sub_df=avg_sub_df.append({"Year":year,'Sub_r_ft_yr':subs_avg},ignore_index = True)

        #Let's excavate drains
        #We sample elevations
        drns['top']=ml.dis.top[drns['i'],drns['j']]

        #Let's convert recarray to pandas
        drns_pd=pd.DataFrame(drns)

        #Drains that will remain in layer 1
        cond=(drns_pd.k==0)&(drns_pd.elev-subs_avg>=drns_pd.PT_bot)
        drns_pd.loc[cond,'elev']=drns_pd[cond]['elev']-subs_avg
        #Drains that will switch from layer 1 to layer 2
        cond=(drns_pd.k==0)&(drns_pd.elev-subs_avg<drns_pd.PT_bot)
        drns_pd.loc[cond,'k']=1
        drns_pd.loc[cond,'elev']=drns_pd[cond]['elev']-subs_avg
        #Drains that will remain in layer 2
        cond=(drns_pd.k==1)&(drns_pd.elev-subs_avg>=drns_pd.TM_bot)
        drns_pd.loc[cond,'elev']=drns[cond]['elev']-subs_avg


        drns=drns_pd.to_records(index=False)

        # LEt's update drains
        ml.drn.stress_period_data[0][['k', 'elev']] = drns[['k', 'elev']]

        # Let's update constant heads

        for layer in range(3):
            ml.bas6.strt[layer] = SLR.loc[year, "2_ft"]
        #Change model working space
        ml.change_model_ws(new_pth=os.path.join(sens, "BAU",str(year)))
        ml.write_input()

        # Let's run MODFLOW
        ml.run_model()

        #Update heads
        h = flopy.utils.HeadFile(os.path.join(ml.model_ws,"Bacon.hds"), model=ml)
        heads=h.get_data()
        #Path for output
        out_BAU_path=os.path.join(sens,"BAU","Output", str(year))
        if not os.path.exists(out_BAU_path):
            os.makedirs(out_BAU_path)
        np.save(os.path.join(out_BAU_path, "heads.npy"), heads)

        drns['h_PT']=heads[0][drns['i'],drns['j']]
        drns['h_TM']=heads[1][drns['i'],drns['j']]
        drns['h_SP']=heads[2][drns['i'],drns['j']]

        #Let's calculate gradients
        #If water table in peat
        cond=np.where(drns['h_PT']>ml.hdry)
        #Gradient is the head difference between the peat and the sand
        #divided by the remaining PT + TM
        drns['Grad'][cond]=(-drns['h_PT'][cond]+drns['h_SP'][cond])/(drns['elev'][cond]-drns['TM_bot'][cond])
        #If water table is in the tidal muc
        cond=np.where(drns['h_PT']==ml.hdry)
        #Gradient is the head difference between the tidal mud and the sand
        #divided by the remaining TM
        drns['Grad'][cond]=(-drns['h_TM'][cond]+drns['h_SP'][cond])/(drns['elev'][cond]-drns['TM_bot'][cond])

        #Directory for shapefiles
        shp_BAU_path=os.path.join(shp_dir,sens,"BAU")
        if not os.path.exists(shp_BAU_path):
            os.makedirs(shp_BAU_path)

        flopy.export.shapefile_utils.recarray2shp(drns,
                                                  geoms=polygons,
                                                  shpname=os.path.join(shp_BAU_path,"DRN_"+str(year)+".shp"),
                                                  epsg=grid.epsg)

       #Let's export subsidence
        flopy.export.utils.export_array(grid, os.path.join(BAU_ras_dir,"Subs_ft_"+str(year)+".tif"), subs)
        #Let's export heads
        flopy.export.utils.export_array(grid, os.path.join(BAU_ras_dir, "H_Pt_ft_" + str(year) + ".tif"), heads[0])
        flopy.export.utils.export_array(grid, os.path.join(BAU_ras_dir, "H_TM_ft_" + str(year) + ".tif"), heads[1])
        flopy.export.utils.export_array(grid, os.path.join(BAU_ras_dir, "H_SP_ft_" + str(year) + ".tif"), heads[2])


        ti=datetime.datetime.now()
        drns_pd=pd.DataFrame(drns)


        #Let's subset to toedrains
        toedrains_dum=pd.merge(drns_pd, toedrains_in, how="right", on=["i","j"])
        #Let's convert to recarray
        toedrains_dum_rec=toedrains_dum.to_records()
        flopy.export.shapefile_utils.recarray2shp(toedrains_dum_rec,
                                                  geoms=polygons_toedrn,
                                                  shpname=os.path.join(shp_BAU_path,"TOEDRNS_"+str(year)+".shp"),
                                                  epsg=grid.epsg)

        drns_pd["Year"]=year
        toedrains_dum["Year"]=year
        drns_pd.to_csv(os.path.join(csv_dir,sens+"_BAU_DRNS"+str(year)+".csv"),index=False)
        toedrains_dum.to_csv(os.path.join(csv_dir,sens+"_BAU_TOEDRNS"+str(year)+".csv"),index=False)


        #Let's export top elevation
        flopy.export.utils.export_array(grid, os.path.join(BAU_ras_dir, "Top_ft_" + str(year) + ".tif"), ml.dis.top[:])


        #Let's export constant heads
        CH_df["CH"]=bas.strt[0][CH]
        CH_rec=CH_df.to_records(index=False)

        pyhf.RPF.rpf(Transects, ml, year, SLR, RPF_BAU_out)

    # Let's create band plots
    pyhf.utils.band_plot(Start_Year-1, End_Year, RPF_BAU_out)

    np_dir_BAU=os.path.join(np_dir,sens,"BAU")
    if not os.path.exists(np_dir_BAU):
        os.makedirs(np_dir_BAU)
    #Let's export average subsidence dataframe
    avg_sub_df.to_csv(os.path.join(np_dir_BAU,"Avg_Sub.csv"),index=False)

    ##Now, let's start with the WLR scenario
    sim_years = End_Year-Start_Year+2
    sedcalc_dum = pyhf.SEDCALC.sedcalc(h2oin=pd.read_csv(os.path.join("Base","WLR","h2oin.csv")),
                                       minin=pd.read_csv(os.path.join("Base", "WLR", "minin.csv")),
                                       orgin=pd.read_csv(os.path.join("Base", "WLR", "orgin.csv")),
                                       porelim=pd.read_csv(os.path.join("Base", "WLR", "porelim.csv")),
                                       endtim=sim_years)
    relelv_dum = sedcalc_dum["relelv"].values
    depth_dum = sedcalc_dum["depth"].values
    porg_dum = sedcalc_dum["porg"].values
    bulkd_dum = sedcalc_dum["bulkd"].values
    acc_rate_dum = pd.DataFrame(columns=["Year", "Yearly Accretion (cm)", "Yearly Accretion (ft)"])
    totC_dum = pd.DataFrame(columns=["Year", "Yearly Accretion (cm)", "gC/cm3", "gC/cm2"])
    year = Start_Year

    for i in range(sim_years - 1):
        # Index for the backwards sum
        i_bw = sim_years - i - 1

        #        else:
        acc_rate_dum_dum = pd.DataFrame({"Year": [year],
                                         "Yearly Accretion (cm)": [(depth_dum[i_bw] - depth_dum[i_bw - 1])],
                                         "Yearly Accretion (ft)": [
                                             (depth_dum[i_bw] - depth_dum[i_bw - 1]) * 0.0328084]})
        totC_dum_dum = pd.DataFrame({"Year": [year],
                                     "Yearly Accretion (cm)": [(depth_dum[i_bw] - depth_dum[i_bw - 1])],
                                     "gC/cm3": [porg_dum[i_bw] * bulkd_dum[i_bw] / 2],
                                     "gC/cm2": [
                                         (depth_dum[i_bw] - depth_dum[i_bw - 1]) * porg_dum[i_bw] * bulkd_dum[
                                             i_bw] / 2]})

        acc_rate_dum = acc_rate_dum.append(acc_rate_dum_dum, ignore_index=True)
        totC_dum = totC_dum.append(totC_dum_dum, ignore_index=True)
        year = year + 1
    #Let's add uncertainty
    totC_dum["Yearly Accretion (cm)"]=np.maximum(0,totC_dum["Yearly Accretion (cm)"]+uncert_SEDCALC[sens])
    acc_rate_dum["Yearly Accretion (cm)"]=np.maximum(0,acc_rate_dum["Yearly Accretion (cm)"]+uncert_SEDCALC[sens])

    # Let's export the total carbon table as a csv file
    totC_dum.to_csv(os.path.join(sens,"WLR","SEDCALC_totC_output.csv"), index=False)
    # Let's export the SEDCALC accretion rates
    acc_rate_dum.to_csv(os.path.join(sens,"SEDCALC_Accretion.csv"), index=False)
    sedcalc_ts = acc_rate_dum

    WLR_ras_dir=os.path.join(ras_dir,sens,"WLR")
    if not os.path.exists(WLR_ras_dir):
        os.makedirs(WLR_ras_dir)

    # Directory for shapefiles
    shp_WLR_path = os.path.join(shp_dir, sens, "WLR")
    if not os.path.exists(shp_WLR_path):
        os.makedirs(shp_WLR_path)

    RPF_WLR_out = os.path.join(wdir, sens, "WLR", "Output", "RPF")
    if not os.path.exists(RPF_WLR_out):
        os.makedirs(RPF_WLR_out)

    # 3. Let's loop through years now
    for year in range(Start_Year, End_Year + 1):
        # 3.1 Initiate MODFLOW for year 1
        if year == Start_Year:
            ml = flopy.modflow.Modflow.load('Base/WLR/MF_inputs/Bacon.nam')
            ml.change_model_ws(new_pth=os.path.join(sens, "WLR", str(year)))
            ml.pcg.rclose = 864
            ml.exe_name = mf_exe_path
            ml.modelgrid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)
            #Let's set rice and wetland to constant head
            ml.bas6.ibound[0][rice_mask] = -1
            ml.bas6.ibound[0][wetland_mask] = -1
            # Let's add constant heads for rice
            ml.bas6.strt[0][rice_mask] = ml.dis.top[rice_mask]


            ml.write_input()
            ml.run_model()
            pyhf.RPF.rpf(Transects, ml, year-1, SLR, RPF_WLR_out)


        # Let´s add wetlands accretion to land surface
        ml.dis.top[wetland_mask] = np.minimum(ml.dis.top[wetland_mask] + float(sedcalc_ts.loc[sedcalc_ts.Year == year,
                                                                                   "Yearly Accretion (ft)"]),
                                          SLR.loc[SLR.Year==year, "2_ft"].values[0])
        # Let´s add constant head cells for wetlands
        ml.bas6.strt[0][wetland_mask] = ml.dis.top[wetland_mask]

        # Peat thickness
        PT_thck = ml.dis.gettop()[0]-ml.dis.getbotm()[0]

        # Depth to groundwater
        h = flopy.utils.HeadFile(os.path.join(ml.model_ws, "Bacon.hds"), model=ml)
        heads = h.get_data()
        wt = flopy.utils.postprocessing.get_water_table(heads, masked_values=[999]).data
        DTW = np.maximum(0, (ml.dis.top[:] - wt))

        # We sample elevations
        drns['top'] = ml.dis.top[drns['i'], drns['j']]

        # Let's update constant heads

        for layer in range(3):
            ml.bas6.strt[layer][CH] = SLR.loc[SLR.Year==year, "2_ft"]

        ml.change_model_ws(new_pth=os.path.join(os.path.join(sens, "WLR",str(year))))
        ml.write_input()
        # Let's run MODFLOW
        ml.run_model()

        #Update heads
        h = flopy.utils.HeadFile(os.path.join(ml.model_ws, "Bacon.hds"), model=ml)
        heads=h.get_data()
        drns['h_PT']=heads[0][drns['i'],drns['j']]
        drns['h_TM']=heads[1][drns['i'],drns['j']]
        drns['h_SP']=heads[2][drns['i'],drns['j']]

        #Let's calculate gradients
        #If water table in peat
        cond=np.where(drns['h_PT']>ml.hdry)
        #Gradient is the head difference between the peat and the sand
        #divided by the remaining PT + TM
        drns['Grad'][cond]=(-drns['h_PT'][cond]+drns['h_SP'][cond])/(drns['elev'][cond]-drns['TM_bot'][cond])
        #If water table is in the tidal muc
        cond=np.where(drns['h_PT']==ml.hdry)
        #Gradient is the head difference between the tidal mud and the sand
        #divided by the remaining TM
        drns['Grad'][cond]=(-drns['h_TM'][cond]+drns['h_SP'][cond])/(drns['elev'][cond]-drns['TM_bot'][cond])



        flopy.export.shapefile_utils.recarray2shp(drns,
                                                  geoms=polygons,
                                                  shpname=os.path.join(shp_WLR_path,"WLR_DRN_"+str(year)+".shp"),
                                                  epsg=grid.epsg)
        #Let's export heads
        flopy.export.utils.export_array(grid, os.path.join(WLR_ras_dir, "H_Pt_ft_" + str(year) + ".tif"), heads[0])
        flopy.export.utils.export_array(grid, os.path.join(WLR_ras_dir, "H_TM_ft_" + str(year) + ".tif"), heads[1])
        flopy.export.utils.export_array(grid, os.path.join(WLR_ras_dir, "H_SP_ft_" + str(year) + ".tif"), heads[2])

        drns_pd = pd.DataFrame(drns)

        # Let's subset to toedrains
        toedrains_dum = pd.merge(drns_pd, toedrains_in, how="right", on=["i", "j"])

        # Let's convert to recarray
        toedrains_dum_rec = toedrains_dum.to_records()

        flopy.export.shapefile_utils.recarray2shp(toedrains_dum_rec,
                                                  geoms=polygons_toedrn,
                                                  shpname=os.path.join(shp_WLR_path,"WLR_TOEDRNS_"+str(year)+".shp"),
                                                  epsg=grid.epsg)
        drns_pd["Year"] = year
        toedrains_dum["Year"] = year
        drns_pd.to_csv(os.path.join(csv_dir, "WLR_DRNS" + str(year) + ".csv"), index=False)
        toedrains_dum.to_csv(os.path.join(csv_dir,"WLR_TOEDRNS"+str(year)+".csv"),index=False)

        #Let's export top elevation
        flopy.export.utils.export_array(grid, os.path.join(WLR_ras_dir, "Top_ft_" + str(year) + ".tif"), ml.dis.top[:])

        pyhf.RPF.rpf(Transects, ml, year, SLR, RPF_WLR_out)

    # Let's create band plots
    pyhf.utils.band_plot(Start_Year-1, End_Year, RPF_WLR_out)


    #Let's export shapefile of cross sections for first and last time step
    pyhf.utils.shp_df_join(RPF_BAU_out,["Transects_"+str(Start_Year-1)+".csv",
                                        "Transects_" + str(End_Year) + ".csv"]
                           ,os.path.join(shp_dir,"RPF","Transect_Cells_v3.shp"),
                           os.path.join(shp_dir, "RPF"),
                           sens + "_BAU_"
                           )

    pyhf.utils.shp_df_join(RPF_WLR_out,["Transects_"+str(Start_Year-1)+".csv",
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