# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 16:20:43 2021

@author: JDiaz
"""

import sys
import os
import datetime
import pandas as pd
import numpy as np
import pickle
import operator
import subprocess
import flopy
import matplotlib as mpl
import time
from math import exp, log, log10, isnan
import matplotlib.pyplot as plt
import csv
import platform

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 



from multiprocessing import Pool,cpu_count


def subcalc_2021_npy(fom,
                 fomund,
                 dtwm,
                 peatdepth,
                 startyr,
                 
                 sim_id='0',
                 st=None,
                 km=None,
                 firstkm=None,
                 vmax=None,
                 bd=None,
                 bdund=None,
                 bdund2=1.3,
                 fomund2=0,
                 massmin=None,
                 
                 tempeff=False,
                 tempc="model",
                 rcp="4.5",
                 initelev=0,
                 output="table",
                 cols="all",
                 silent=True):
    
    '''
    
    This is the main SUBCALC function for processing emisisons and subsidence
    within the Sacramento-San Joaquin delta. This is the updated model as of
    2021 building on the most recent published 2016 model. The original SUBCALC model 
    was developed by Deverel and Leighton (2010).
    
    This model makes the following revisions:
        
    - km equation revised based on collected eddy data and substrate relationship
    - substrate used in Michaelis-Menten in place of FOC
    - water fraction no longer scales emissions; assumes whole unsaturated zone is emitting
    - bulk density euqation for unsaturated zone adjusted to Staten-based relationship with OM
    - bulk density below the water table changed to Drexler 2009 equation for >120cm
    - temperature effect added from Deverel & Rojstaczer 1996
    
    sim_id = unique ID for set of inputs
    fom = fraction of organic matter in unsaturated zone (weighted average); decimal form
    fomund = representative OM for beneath water table
    dtwm = averge DTW in meters
    tempc = initial average soil temperature in celsius to calculate vmax; "model" option uses the Bradford et al. 2019 prediction
    startyr = first year of the simulation
    endyr = last year of the simulation
    st = total substrate in unsaturated zone in gC
    km = Michelis-Menten constant
    vmax = maximum rate; Michelis-Menten parameter
    bd = bulk density in unsaturated zone
    bdund = representative Db for beneath water table
    bdund2 = Db of material beneath the peat
    fomund2 = fraction of organic matter beneath peat
    massmin = mass of the mineral fraction
    massom = mass of the organic fraction
    peatdepth = maximum depth of the peat in cm
    tempeff = flag for applying temperature relationship
    rcp = which climate scenario to use
    initelev = inital elevation
    output = 'table', 'sum', 'average'; 
              table return the full output table of selected columns
              sum returns the sum of selected columns for the period
              average returns the average of selected columns for the period
    cols = select specific columns to return by list; default all columns
    

    '''
    

    
    
    
     

    
    endyr=startyr
    ### Simulation depth and inital elevation ###
    simdpth = dtwm * 100
#    elev = initelev
    
    
    
    ##### Calcualte BD for unsat and sat zones based on OM if not provided #####
    if bd is None:
        bd = -0.354 * np.log(fom*100) + 1.7561
    
    if bdund is None:
        bdund = -0.21* np.log(fomund*100) + 1.01
        
        
    ##### Calculate total substrate from OM if not provided #####
    if st is None:
        st = (fom * 0.5) * bd * simdpth
        
        
    ##### Calculate mineral mass if not provided #####
    if massmin is None:
        massmin = (1-fom) * bd * simdpth
    
    
    
    
    
    
    ### Inital check if water table is below the peat ###
    
#    if peatdepth is not None:
        
        # if DTW exceeds peat depth, change fomund & bdund and recalc substrate
#    if simdpth >= peatdepth:
#        st = (fom * 0.5) * bd * peatdepth
    st[np.where(simdpth >= peatdepth)]=fom[np.where(simdpth >= peatdepth)]*0.5*bd[np.where(simdpth >= peatdepth)]*peatdepth[np.where(simdpth >= peatdepth)]
#        fomund = fomund2
    fomund[np.where(simdpth >= peatdepth)]= fomund2
#        bdund = bdund2
    bdund[np.where(simdpth >= peatdepth)]=bdund2
#        massmin = ((1-fom)*bd)*peatdepth
    massmin[np.where(simdpth >= peatdepth)]=((1-fom[np.where(simdpth >= peatdepth)])*bd[np.where(simdpth >= peatdepth)])*peatdepth[np.where(simdpth >= peatdepth)]

#    else:
#        peatdepth = np.nan
    peatdepth[np.where(simdpth < peatdepth)]=np.nan
    
    
    
    
    
    ##### Get vector of years #####
    # endyr should be the last year of the simulation
#    years = list(range(int(startyr),int(endyr)+1))



    ## Calculate km value
    if km is None:
        km = (-4.2732 * st) + 51.709
#        firstkm = False # flag for resetting km in loop
        firstkm =np.full(km.shape, False)
        
#    if km > 0:
#        firstkm = True
    firstkm[np.where(km > 0)]=True
    
        
    
    
    
    
    
    ### Constants ###
    preexp = 3038800000000
    activ = 73.30
    conconst = 0.00325
    condepth = 150
    

    
    
    
    ### Calculate vmax ###
    if vmax is None:
        
        if tempc == "model":
            
            if startyr < 2050:
                
                if rcp == "4.5":
                    tempc = 0.027961395 * startyr - 39.64990971
                if rcp == "8.5":
                    tempc = 0.035192952 * startyr - 54.18533977
                
            if startyr >= 2050:
                
                if rcp == "4.5":
                    tempc = 0.021424862 * startyr - 26.25001712
                if rcp == "8.5":
                    tempc = 0.049837842 * startyr - 84.20736384
        
        tempk = tempc + 273.15
        vmax = preexp * exp(-activ / (0.0083144 * tempk))
    
    
    
    
    
   
    
    
    
    
    
    
    
    
    
    
    
    ##### Begin SUBCALC loop #####
#    for year in years:
        
    ## output table row ##
    idx = 0
    
    year=startyr
    ### Calculate temperature effect ###
#    if tempeff and year >=2010:
#        
#        if years[0] < 2050:
#            init45 = 0.027961395 * years[0] - 39.64990971
#            init85 = 0.035192952 * years[0] - 54.18533977
            
#        if years[0] >= 2050:
#            init45 = 0.021424862 * years[0] - 26.25001712
#            init85 = 0.049837842 * years[0] - 84.20736384
            
        # calculate temp given year based on Bradford et al., 2019
        # linear interpolations between averages within delta
#        if year < 2050:
            
#            tempyr45 = 0.027961395 * year - 39.64990971
#            tempyr85 = 0.035192952 * year - 54.18533977
            
            # determine change in temp from slopes
            # deltc45 = 0.027961395
            # deltc85 = 0.035192952
            
            # determine cumulative change in temp
#            deltc45 = tempyr45 - init45
#            deltc85 = tempyr85 - init85
            
#        if year >= 2050:
            
#            tempyr45 = 0.021424862 * year - 26.25001712
#            tempyr85 = 0.049837842 * year - 84.20736384
            
            # determine change in temp from slopes
            # deltc45 = 0.02142862
            # deltc85 = 0.049837842
            
#            deltc45 = tempyr45 - init45
#            deltc85 = tempyr85 - init85
        
        
        # choose right temperature based on desired scenario
#        if rcp == "4.5":
#            tempyr = tempyr45
#            deltc = deltc45
#        if rcp == "8.5":
#            tempyr = tempyr85
#            deltc = deltc85
        
        
#        # Assume initial temperature for first year with no increase
#        if idx == 0:
#            delcflux = 0
#        else:
            # Calculate theoretical flux from increase temperature
            # based on Van't Hoff
#            init_flux = output_table.at[0,'cflux']
#            vant_hoff = init_flux * exp( (log(1.94) / 10) * deltc)
            
            # Added cflux from temperature increase #
#            delcflux = vant_hoff - init_flux

    
#        output_table.at[idx,'deltc'] = deltc
#        output_table.at[idx,'tempyr'] = tempyr
    
#    else:
        
    delcflux = 0
        
    
#    output_table.at[idx,'delcflux'] = delcflux
    
    
    
    cflux=np.zeros(km.shape)
    ### Calculate oxidation ###
    ## Set cflux to vmax when km out of range
#    if km <= 0:
        # recalc vmax for increasing temp when at max
#        if tempeff:
#            # vmax = preexp * exp(-activ / (0.0083144 * (tempyr+273.15)))
#            vmax = vmax
#        cflux = vmax + delcflux
    ind_dum=np.where(km <= 0)
    cflux[ind_dum]=vmax + delcflux
    
#    else:

#        cflux = ((vmax * st) / (km + st)) + delcflux
    ind_dum=np.where(km > 0)
    cflux[ind_dum]=((vmax * st[ind_dum]) / (km[ind_dum] + st[ind_dum])) + delcflux
    
    
    omflux = cflux * 2
    delelvc = omflux / (bd * fom)
    
    
    
    
    
    ### Calculate wind erosion (past years) ###
    
    
    
    ### Calculate burning (paste years) ###
    
    
    
    ### Calculate consolidation ###
    delp = delelvc
    delelvcons = conconst * condepth * delp
    
    ### combine sources of subsidence ###
    tdelev = delelvc + delelvcons
#    elev = elev - tdelev
    
    
    
    if peatdepth is not None:
        
        ## Update remaining peat ##
#        peatleft = peatdepth + elev # elev is negative
        peatleft = peatdepth- tdelev
#        output_table.at[idx,'peatdepth'] = peatleft
        
        # if total subsidence plus the DTW exceeds peat depth
#        if -elev+simdpth >= peatdepth:       
#            fomund = fomund2
        ind_dum=np.where(simdpth+tdelev>=peatdepth)
        fomund[ind_dum]=fomund2
#            bdund = bdund2
        bdund[ind_dum]= bdund2
        
    else:
        peatleft = simdpth
    
    
    ### Update variables ###
    volloss = tdelev
    massomlft = (st*2) - omflux 
    massom = massomlft + (fomund * bdund * volloss)
    massmin = massmin + ((1 - fomund) * bdund * volloss)
    
    fom = massom / (massom + massmin)
    
#    if simdpth <= peatleft: 
#        bd = (massom + massmin) / simdpth
    ind_dum=np.where(simdpth <= peatleft)
    bd[ind_dum]=(massom[ind_dum] + massmin[ind_dum]) / simdpth[ind_dum]
    
#    else:
#        bd = (massom + massmin) / peatleft
    ind_dum=np.where(simdpth > peatleft)
    bd[ind_dum] = (massom[ind_dum] + massmin[ind_dum]) / peatleft[ind_dum]
    
    st = massom/2
    kmnew = (-4.2732 * st) + 51.709 # used to "turn off" clfux = vmax
    
 #   if kmnew > 0 and firstkm==False:

#        km = kmnew
    ind_dum=np.where((kmnew > 0) & (firstkm==False))
    km[ind_dum]=kmnew[ind_dum]
#        firstkm = True
    firstkm[ind_dum]=True

    Output={'km':km,
            'vmax':vmax,
            'fomund':fomund,
            'bdund':bdund,
            'cflux':cflux,
            'delelvc':delelvc,
            'delelvcons':delelvcons,
            'tdelev':tdelev,
            'volloss':volloss,
            'massomlft':massomlft,
            'massom':massom,
            'massmin':massmin,
            'fom':fom,
            'bd':bd,
            'st':st,
            'firstkm':firstkm}
    return Output



t0=datetime.datetime.now()

#We set directory
os.chdir(r"\\hydro-nas\Team\Projects\5630_DSC\Bacon Island Model\Model\Final")

#Numpy arrays
np_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\Numpy"

#Shapefiles
shp_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\Final"

#Rasters
ras_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\Final"

###Number of processors
n_cpu=cpu_count()

##Sea Level Rise timeseries
SLR=pd.read_csv("SLR.csv",index_col=0)

#First year of simulation
Start_Year=2018

#Last year of simulation
End_Year=2070

#Let's load model
ml = flopy.modflow.Modflow.load('MF_inputs/Bacon.nam')

grid = ml.modelgrid
grid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)
nrg=ml.nrow
ncg=ml.ncol

#Let's export grid shapefile
ml.dis.export(os.path.join(shp_dir,"Grid.shp"))

#Let's export hydraulic conductivities raster
lpf = flopy.modflow.ModflowLpf.load('MF_inputs/Bacon.lpf', ml)
flopy.export.utils.export_array(grid, os.path.join(ras_dir, "HK_Lay1.tif"), lpf.hk[0][:])

#Let's plot cross section


fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(1, 1, 1)
line = flopy.plot.plotutil.shapefile_get_vertices(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\GW Model\25ft\CrossSection\X_Sec2.shp")
xsect = flopy.plot.PlotCrossSection(model=ml, line={"line": line[0]})

patches = xsect.plot_ibound()
patches = xsect.plot_bc("DRN", color="pink")
#patches = xsect.plot_bc("CHD", color="red")
linecollection = xsect.plot_grid()
#cb = plt.colorbar(csa, shrink=0.75)
plt.savefig(r"\\hydro-nas\Team\Projects\5630_DSC\Report\Longer_Version\CrossSection\XSec.svg")


subsidence=np.zeros(shape=(nrg,ncg))

class_length=int(nrg/n_cpu)

#Let's get active cells
bas = flopy.modflow.ModflowBas.load('MF_inputs/Bacon.bas', ml)

ibound=np.array(bas.ibound[0][:])
Active_cells=np.where(ibound!=0)

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


#SLR time series
SLR=pd.read_csv("SLR.csv",index_col=0)

# 1. Import initial soil properties

SC_Input={'fom':np.load('fom_0.npy'),
          'fomund':np.load('fomund_0.npy'),
          'st':np.load('st_0.npy'),
          'km':np.load('km_0.npy'),
          'firstkm':np.load('firstkm_0.npy'),
          'vmax':np.load('vmax_0.npy'),
          'bd':np.load('bd_0.npy'),
          'bdund':np.load('bdund_0.npy'),
          'massmin':np.load('massmin_0.npy')}

#2. Import elevations

#Dataframe with average subsidence
avg_sub_df=pd.DataFrame(columns = ['Year', 'Sub_r_ft_yr'])

#3. Let's loop through years now


for year in range(Start_Year,End_Year+1):
    
    # 3.1 Initiate MODFLOW for year 1
    if year==Start_Year:
        ml = flopy.modflow.Modflow.load('MF_inputs/Bacon.nam')
        ml.write_input()
        subprocess.check_output(["mf2005", "Bacon_fix.nam"])
        
    #3.2 Let's run SUBCALC
    
    #Peat thickness
    PT_thck=ml.dis.thickness[0]
    
    #Depth to groundwater
    h = flopy.utils.HeadFile("Bacon.hds", model=ml)
    heads=h.get_data()

    ml.modelgrid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)
    #levels = np.arange(int(np.min(heads[2]))+1, int(np.sort(np.unique(heads[0]))[-2]), 2)
    #levels = np.arange(-20,0, 5)
    levels =[-16,-14,-12,-8,-4,0,4]
    fig = plt.figure(figsize=(18, 5))
    ax = fig.add_subplot(1, 1, 1)
    line = flopy.plot.plotutil.shapefile_get_vertices(
        r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\GW Model\25ft\CrossSection\X_Sec2.shp")
    xsect = flopy.plot.PlotCrossSection(model=ml, line={"line": line[0]})

    patches = xsect.plot_ibound()
    patches = xsect.plot_bc("DRN", color="pink")
    # patches = xsect.plot_bc("CHD", color="red")
    linecollection = xsect.plot_grid()
    contour_set = xsect.contour_array(
        heads, masked_values=[999.0], head=heads, levels=levels, colors="k"
    )
    plt.clabel(contour_set, fmt="%.1f", colors="k", fontsize=11)
    plt.savefig(os.path.join(shp_dir,"XSec_"+str(year)+".svg"))

    wt = flopy.utils.postprocessing.get_water_table(heads=heads, nodata=np.min(heads[0]))
    DTW=np.maximum(0,(ml.dis.top[:]-wt))
    t0 = datetime.datetime.now()
    SC_Input=subcalc_2021_npy(SC_Input['fom'],
                              SC_Input['fomund'],
                              DTW/3.28084,
                              PT_thck/3.28084*100,
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
    
    #average subsidence
    subs_avg=np.average(subs[subs>0])
    avg_sub_df=avg_sub_df.append({"Year":year,'Sub_r_ft_yr':subs_avg},ignore_index = True)
    
    #Let's excavate drains
    #We sample elevations
    drns['top']=ml.dis.top[drns['i'],drns['j']]
    #Let's convert recarray to pandas
    drns_pd=pd.DataFrame(drns)
    #Drains that will remain in layer 1
    #cond=(drns['k']==0)&(drns['elev']-subs_avg>=drns['PT_bot'])
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

    #Convert back to recarray
    
    #LEt's update drains
    ml.drn.stress_period_data[0][['k','elev']]=drns[['k','elev']]
    ml.write_input()
    #Let's update constant heads
    
    for layer in range(3):
        bas.strt[layer][CH]=SLR.loc[year,"2_ft"]
    bas.write_file()
    #Let's run MODFLOW

    subprocess.check_output(["mf2005", "Bacon_fix.nam"])
    
    #Update heads
    h = flopy.utils.HeadFile("Bacon.hds", model=ml)
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
                                              shpname=os.path.join(shp_dir,"DRN_"+str(year)+".shp"), 
                                              epsg=grid.epsg)
    #Let's export subsidence
    flopy.export.utils.export_array(grid, os.path.join(ras_dir,"Subs_ft_"+str(year)+".tif"), subs)
    #Let's export heads
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "H_Pt_ft_" + str(year) + ".tif"), heads[0])
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "H_TM_ft_" + str(year) + ".tif"), heads[1])
    flopy.export.utils.export_array(grid, os.path.join(ras_dir, "H_SP_ft_" + str(year) + ".tif"), heads[2])
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
    drns_pd.to_csv(os.path.join(np_dir,"DRNS"+str(year)+".csv"),index=False)
    toedrains_dum.to_csv(os.path.join(np_dir,"TOEDRNS"+str(year)+".csv"),index=False)
    #drns.tofile(os.path.join(np_dir,"DRNS"+str(year)+".csv"),sep=",")
    print(year,ti-t0)

#Let's export average subsidence dataframe
avg_sub_df.to_csv(os.path.join(np_dir,"Avg_Sub.csv"),index=False)
         
        
        
