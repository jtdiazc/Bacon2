import sys
#path to the hydrofocus functions module
sys.path.insert(0, r'\\hydro-nas\Team\Projects\5630_DSC\Codes')
import pyhf
#Path to flopy module
sys.path.insert(0, r'\\hydro-nas\Team\Projects\5630_DSC\Codes\flopy')
import flopy
import os
import pandas as pd
import numpy as np
import rasterio
import subprocess
import datetime

#We set directory
wdir=r"\\hydro-nas\Team\Projects\5630_DSC\Bacon Island Model\Model\20221111"
os.chdir(wdir)

#Numpy arrays
np_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\Numpy"

#Shapefiles
shp_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\Final"

#Rasters
ras_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\Final"

#CSVs
csv_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV"

#MODFLOW output files
mf_out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Bacon Island Model\Model\20221111\MF_files"

##Sea Level Rise timeseries
SLR=pd.read_csv("SLR.csv",index_col=0)

#First year of simulation
Start_Year=2018

#Last year of simulation
End_Year=2070

#Let's load model
ml = flopy.modflow.Modflow.load('MF_inputs/Bacon.nam')
ml.modelgrid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)
nrg=ml.nrow
ncg=ml.ncol

#Let's get active cells
bas = flopy.modflow.ModflowBas.load('MF_inputs/Bacon.bas', ml)

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

pyhf.flopyf.df_to_shp(grid,CH_df,shp_dir,"CH.shp")

#SLR time series
SLR=pd.read_csv("SLR.csv",index_col=0)

#SUBCALC input
SC_Input={'fom':np.load('fom_0.npy'),
          'fomund':np.load('fomund_0.npy'),
          'st':np.load('st_0.npy'),
          'km':np.load('km_0.npy'),
          'firstkm':np.load('firstkm_0.npy'),
          'vmax':np.load('vmax_0.npy'),
          'bd':np.load('bd_0.npy'),
          'bdund':np.load('bdund_0.npy'),
          'massmin':np.load('massmin_0.npy')}

#Let's set nas to 0
SC_Input['fom'][np.isnan(SC_Input['fom'])]=0

#2. Import elevations

#Dataframe with average subsidence
avg_sub_df=pd.DataFrame(columns = ['Year', 'Sub_r_ft_yr'])

##Let's export shapefile of levees
Levees_df = pd.DataFrame({"row":levees[0],
                   "col":levees[1]})
pyhf.flopyf.df_to_shp(grid,Levees_df,shp_dir,"Levees.shp")


##BAU

#3. Let's loop through years now
for year in range(Start_Year,End_Year+1):
    # 3.1 Initiate MODFLOW for year 1
    if year==Start_Year:
        ml = flopy.modflow.Modflow.load('MF_inputs/Bacon.nam')
        ml.write_input()
        t0 = datetime.datetime.now()
        subprocess.check_output(["mf2005", "Bacon_fix.nam"])
        tf = datetime.datetime.now()
        print(tf-t0)
        # Peat thickness
        PT_thck = ml.dis.gettop()[0]-ml.dis.getbotm()[0]
        PT_thck[Inactive_cells] = 0
        PT_thck[levees] = 0

    #Depth to groundwater
    h = flopy.utils.HeadFile("Bacon.hds", model=ml)
    heads=h.get_data()
    #ml.BAS6.strt[0][:]=heads[0]
    #ml.BAS6.strt[1][:] = heads[1]
    #ml.BAS6.strt[2][:] = heads[2]

    #ml.pcg.damp=0.95
    #ml.pcg.relax = 0.99


    ml.modelgrid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)

    wt=flopy.utils.postprocessing.get_water_table(heads,masked_values=[999]).data

    DTW = np.maximum(0, (ml.dis.top[:] - wt))

    t0 = datetime.datetime.now()

    dtwm = DTW / 3.28084
    peatdepth = PT_thck / 3.28084 * 100

    # Let's export DTW raster
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "DTW_m_" + str(year) + ".tif"), dtwm)

    # Let's export Peat thickness cm
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "PT_thck_cm_" + str(year) + ".tif"), peatdepth)

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

    tf = datetime.datetime.now()
    print(tf-t0)
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


    ml.change_model_ws(new_pth=os.path.join(mf_out_path,str(year)))
    ml.write_input()
    ml.change_model_ws(new_pth=wdir)
    ml.write_input()
    print("Started MODFLOW at: "+str(datetime.datetime.now()))
    t0 = datetime.datetime.now()
    #Let's run MODFLOW
    subprocess.check_output(["mf2005", "Bacon_fix.nam"])

    tf = datetime.datetime.now()

    print(tf-t0)

    #Update heads
    h = flopy.utils.HeadFile("Bacon.hds", model=ml)
    heads=h.get_data()
    np.save(os.path.join(mf_out_path, str(year), "heads.npy"), heads)

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
                                              shpname=os.path.join(shp_dir,"DRN_"+str(year)+".shp"),
                                              epsg=grid.epsg)
   #Let's export subsidence
    flopy.export.utils.export_array(grid, os.path.join(ras_dir,"Subs_ft_"+str(year)+".tif"), subs)
    #Let's export heads
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "H_Pt_ft_" + str(year) + ".tif"), heads[0])
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "H_TM_ft_" + str(year) + ".tif"), heads[1])
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "H_SP_ft_" + str(year) + ".tif"), heads[2])

    #Let's export head contours

    ti=datetime.datetime.now()
    drns_pd=pd.DataFrame(drns)


    #Let's subset to toedrains
    toedrains_dum=pd.merge(drns_pd, toedrains_in, how="right", on=["i","j"])
    #Let's convert to recarray
    toedrains_dum_rec=toedrains_dum.to_records()
    flopy.export.shapefile_utils.recarray2shp(toedrains_dum_rec,
                                              geoms=polygons_toedrn,
                                              shpname=os.path.join(shp_dir,"TOEDRNS_"+str(year)+".shp"),
                                              epsg=grid.epsg)

    drns_pd["Year"]=year
    toedrains_dum["Year"]=year
    drns_pd.to_csv(os.path.join(csv_dir,"BAU_DRNS"+str(year)+".csv"),index=False)
    toedrains_dum.to_csv(os.path.join(csv_dir,"BAU_TOEDRNS"+str(year)+".csv"),index=False)
    #drns.tofile(os.path.join(np_dir,"DRNS"+str(year)+".csv"),sep=",")

    #Let's export top elevation
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "Top_ft_" + str(year) + ".tif"), ml.dis.top[:])


    #Let's export constant heads
    CH_df["CH"]=bas.strt[0][CH]
    CH_rec=CH_df.to_records(index=False)

#    flopy.export.shapefile_utils.recarray2shp(CH_rec,
#                                              geoms=polygons_CH,
#                                              shpname=os.path.join(shp_dir,"CH_"+str(year)+".shp"),
#                                              epsg=grid.epsg)
    print(year,ti-t0)

#Let's export average subsidence dataframe
avg_sub_df.to_csv(os.path.join(np_dir,"Avg_Sub.csv"),index=False)