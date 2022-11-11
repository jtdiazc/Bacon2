import sys
#path to the hydrofocus functions module
sys.path.insert(0, r'\\hydro-nas\Team\Projects\5630_DSC\Codes\pyhf')
import pyhf
#Path to flopy module
sys.path.insert(0, r'\\hydro-nas\Team\Projects\5630_DSC\Codes\flopy')
import flopy
import os
import pandas as pd
import numpy as np
import rasterio
import subprocess

#We set directory
os.chdir(r"\\hydro-nas\Team\Projects\5630_DSC\Bacon Island Model\Model\20221111")

#Numpy arrays
np_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\Numpy"

#Shapefiles
shp_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\Final"

#Rasters
ras_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\Final"

#CSVs
csv_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV"

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
toedrains_in_rec=toedrains_in.to_records

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

vertices = []
for row, col in zip(CH_df["row"], CH_df["col"]):
    vertices.append(grid.get_cell_vertices(row, col))
polygons_CH = [flopy.utils.geometry.Polygon(vrt) for vrt in vertices]

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

##BAU

#3. Let's loop through years now
for year in range(Start_Year,End_Year+1):
    # 3.1 Initiate MODFLOW for year 1
    if year==Start_Year:
        ml = flopy.modflow.Modflow.load('MF_inputs/Bacon.nam')
        ml.write_input()
        subprocess.check_output(["mf2005", "Bacon_fix.nam"])
        # Peat thickness
        PT_thck = ml.dis.gettop()[0]-ml.dis.getbotm()[0]
        PT_thck[Inactive_cells] = 0
        PT_thck[levees] = 0