# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 16:20:43 2021

@author: JDiaz
"""

import os
import datetime
import pandas as pd
import numpy as np
import pickle
import operator
import subprocess
import flopy
import time
from math import exp, log, log10, isnan
import matplotlib.pyplot as plt
import csv

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 



from multiprocessing import Pool,cpu_count


def rtemin(y,
           tdrnge=300,
           mw=0):
    # *************************************************************
    # *                       RTEMIN                              *
    # * Subroutine for determining the rate of mineral sed input  *
    # *************************************************************
    #
    #
    # this section calculates the amount of mineral sediment
    # input each year - based on the relative elevation of the
    # marsh surface.  Y is the relative elevation of
    # the core at at given time.
    # tdrnge is the tidal range in meters
    # mhw is the relative elevation of mhw in meters (0 in this case)
    #
    tdhght = (y - mw) / (tdrnge / 2.0)
    if tdhght <= 0:
        rtemin = 1.0
    else:
        rtemin = 1 - (min(tdhght, 1.0))
    return rtemin


def poresp(c, k1):
    # *****************************************************
    # *            PORESP                                 *
    # * Subroutine for determining changes in pore space  *
    # *****************************************************
    #
    # the next section calculates the pore space for each section.
    # This is where changes due to compaction occur. I am assuming that
    # all of the pore spaces are filled with water, and that any compaction
    # is due to the loss of water and decrease in pore space volume.
    # Pore space is assumed to be a function of the amount of material
    # (both organic and mineral) that is in a given section, as well
    # as the mass above that particular section.
    #
    # a = totorg(t2)
    # b = min(t2)
    # c = densabv(t2)
    # d = oldorg(t2)
    # e = min(t2)
    #
    # k1 is a constant that affects the curve for compaction
    # k2 affects the relative importance of organic versus mineral
    # matter in determining pore space.
    # K2 > 1 - organic matter more important
    # k2 < 1 - organic matter less important
    # k2 = 1 - organic and mineral matter the same
    #
    #
    # p1 & p2 are just temporary variable to make calculations easier.
    #
    poresp = 1 - (c / (k1 + c))
    #	k2 = 0.1
    #	poresp = (1/(1+(k1*c)))
    #
    # everything below here has been commented out in order to
    # "simplify" the calculation of pore space.
    #
    #  	write(6,2005) a,b,c,d,e
    #  	if (((k2*d)+e).le.0) then
    #		p2 = 1
    #	else
    #		p2 = sqrt(((k2*a)+b)/((k2*d)+e))
    #	end if
    #	poresp = p1*p2

    #     Last change:  SJD  23 Aug 2007    4:35 pm
    return poresp


def rtprod(d):
    # c **********************************************
    # c *  RTPROD                                    *
    # c *  Root production subroutine                *
    # c **********************************************
    # c
    # c
    # c  What follows is a subroutine for determining organic
    # c  production at various time / depths
    # c  1/11/94 - I am changing this so root production decreases \
    # c  exponentially with depth.
    # c  sed2.for has the old version of root production.
    # c
    # c
    # c parameters:
    # c depth(t2) is the only parameter - it is passed
    # c as the single variable "d"
    # c
    # c variables:
    # c	undpro - total underground production (g/cm^3)
    # c	kdist - controls the decay of the root production function
    # c
    undpro = 0.06
    kdist = 0.40
    rtprod = np.exp(-kdist * d) * kdist * undpro

    # c
    # c
    # c this section calculates root production at a particular depth
    # c based on the parameters/curve that are designated above.
    # c
    # c

    # this is sed5.for  -  sediment accretion program for fortran
    # this version uses an exponentially decreasing underground production function
    # this change is in the subroutine - rtprod
    #
    # additionally - this version uses a new decomposition model - still 3
    # different rates - but the rates for each year class are from an
    # exponential decay curve.
    #
    # this version (4/20/2007) reads organic and mineral inputs from a file
    # declaring variables
    #
    #      real*8 org, min, orgden, minden, h2oden, pore, minin, orgin
    #      real*8 orgbd, minbd, bulkd, porg, depth, massabv, intelv, relelv
    #      real*8 slr, subsid, totorg, totvol, orgvol, minvol, densabv
    #      real*8 dt, mindev, porelim, refrac,h2oin
    #      real*8 pctpore, tpore
    #      real*8 acc1000, org1000, min1000, acc4000, org4000, min4000
    #      real*8 acc4900, org4900, min4900, finelv
    #      integer time, t2, endtim
    #      dimension org(0:7000,4), min(0:7000), pore(0:7000), totvol(7000)
    #      dimension orgbd(7000), minbd(7000), porg(7000), minvol(7000)
    #      dimension depth(0:7000),  massabv(0:7000), densabv(0:7000)
    #      dimension bulkd(7000),relelv(0:7000),totorg(0:7000),orgvol(7000)
    #      dimension pctpore(0:7000)
    #      dimension tpore(0:7000)
    #      dimension minin(7000)
    #      dimension orgin(7000)
    #      dimension h2oin(7000)
    #      dimension porelim(7000)
    # c
    return rtprod


def decmp1(d, kdec1):
    # c
    # c *************************************************
    # c *           DECMP1                              *
    # c * Decomposition of "youngest" organic material  *
    # c *************************************************
    # c
    # c
    # c the next section is the decomposition subroutine FOR 1st year org matter
    # c it is a function that gives a decomposition rate (from 0 to 1)
    # c based on the depth of each section.
    # c
    # c decomp is a RATE (units g lost/g present) so it has
    # c to be multiplied by the organic mass (org) of each section
    # c
    # c as with the production function the only thing that determines
    # c this is the depth.
    # c
    # c Variables:
    # c 	mx1 - maximum rate of decay for this age class
    # c	kdec1 - k for exponential decay curve for this
    # c 		age class decompostion curve
    mx1 = 0.92

    #	decmp1 = (exp(-kdec1*d))*mx1
    decmp1 = np.exp(-kdec1 * d) * mx1
    return decmp1


def decmp2(d):
    # c
    # c
    # c
    # c *************************************************
    # c *           DECMP2                              *
    # c * Decomposition of "medium" organic material    *
    # c *************************************************
    # c
    # c
    # c
    # c the next section is the decomposition subroutine FOR 2nd year org matter
    # c it is a function that gives a decomposition rate (from 0 to 1)
    # c based on the depth of each section.
    # c
    # c decomp is a RATE (units g lost/g present) so it has
    # c to be multiplied by the organic mass (org) of each section
    # c
    #	real*8 function decmp2(d)
    #	real*8 mx2, kdec2, d
    # c
    # c as with the production function the only thing that determines
    # c this is the depth.
    # c
    # c Variables:
    # c 	mx2 - maximum rate of decay for this age class
    # c	kdec2 - k for exponential decay curve for this
    # c 		age class decompostion curve
    mx2 = 0.37
    kdec2 = 0.57
    decmp2 = (np.exp(-kdec2 * d)) * mx2
    return decmp2


def decmp3(d):
    # c
    # c
    # c
    # c *************************************************
    # c *           DECMP3                              *
    # c * Decomposition of "oldest" organic material    *
    # c *************************************************
    # c
    # c
    # c
    # c the next section is the decomposition subroutine FOR old org matter
    # c it is a function that gives a decomposition rate (from 0 to 1)
    # c based on the depth of each section.
    # c
    # c The rates are LOWEST for this group of organic material
    # c
    # c decomp is a RATE (units g lost/g present) so it has
    # c to be multiplied by the organic mass (org) of each section
    # c
    #	real*8 mx3, kdec3, d
    # c
    # c as with the production function the only thing that determines
    # c this is the depth.
    # c
    # c Variables:
    # c 	mx3 - maximum rate of decay for this age class
    # c	kdec3 - k for exponential decay curve for this
    # c 		age class decompostion curve
    mx3 = 0.16
    kdec3 = 0.1
    # c
    decmp3 = (np.exp(-kdec3 * d)) * mx3
    #	end
    return decmp3


def sedcalc(endtim=40,
            h2oden=1.00,
            h2oin=pd.read_csv("h2oin.csv"),  # Initial pore space (fraction)
            #################################This file needs to be in the working
            #################################directory
            intelv=-26,
            k1=2.5,  # Consolidation constant
            mindev=0.0,
            minden=2.61,  # Mineral particle density (g cm-2)
            minin=pd.read_csv("minin.csv"),  ##Surface mineral matter deposition (g/cm2)
            #################################This file needs to be in the working
            #################################directory
            orgden=1.14,  # Organic particle density (g cm-2)
            orgin=pd.read_csv("orgin.csv"),  # Surface organic matter deposition (g/cm2)
            #################################This file needs to be in the working
            #################################directory
            porelim=pd.read_csv("porelim.csv"),  # Final pore space (fraction)
            #################################This file needs to be in the working
            #################################directory
            refrac=0.4,
            slr=0.0,  # Sea level rise (cm yr -1)
            strint=5.0,
            strdev=0.3,
            subsid=0.0,  # Subsidence (cm yr -1) from gas withdrawal
            kdec1=0.41
            ):
    # h20in is a percent.  It has to be converted to a volume
    # to be useful for calculations.  The conversion from % to volume is:
    # porespace volume = ((%)/(1-%))*(minvol + orgvol)
    # Let's initialize bulkd
    bulkd = np.zeros(endtim + 1)
    # Let's initialize densabv
    densabv = np.zeros(endtim + 1)

    # Let's initialize depth
    depth = np.zeros(endtim + 1)

    # Let's initialize massabv
    massabv = np.zeros(endtim + 1)

    # Let's initialize minbd
    minbd = np.zeros(endtim + 1)

    # Let's initialize mini (min is a reserved name in Python)
    mini = np.zeros(endtim + 1)

    # Let's initialize minvol
    minvol = np.zeros(endtim + 1)

    # Let's initialize org
    org = np.zeros((endtim + 1, 5))

    # Let's initialize orgbd
    orgbd = np.zeros(endtim + 1)

    # Let's initialize orgvol
    orgvol = np.zeros(endtim + 1)

    # Let's initialize pctpore
    pctpore = np.zeros(endtim + 1)

    # Let's initialize pore
    pore = np.zeros(endtim + 1)

    # Let's initialize porg
    porg = np.zeros(endtim + 1)

    # Let's initialize relelv
    relelv = np.zeros(endtim + 1)

    # Let's initialize totorg
    totorg = np.zeros(endtim + 1)

    # Let's initialize totvol
    totvol = np.zeros(endtim + 1)

    # Let's initialize tpore
    tpore = np.zeros(endtim + 1)

    # ************************************************
    # this is the beginning of the main control loop *
    # ************************************************

    for time in range(1, endtim + 1):

        # c this section moves all values down one section
        # c before the next round of growth, new input and decomposition
        # c

        for t2 in range(time - 1, -1, -1):
            org[t2 + 1, 4] = org[t2, 4]
            org[t2 + 1, 3] = org[t2, 3] + org[t2, 2]
            org[t2 + 1, 2] = org[t2, 1]
            org[t2 + 1, 1] = 0
            mini[t2 + 1] = mini[t2]
            pctpore[t2 + 1] = pctpore[t2]
            tpore[t2 + 1] = tpore[t2]

        # c ************************************************************
        # c * these are the new inputs of material onto the surface of *
        # c * the marsh (into the first position in the array).        *
        # c ************************************************************

        org[1, 1] = float(orgin.loc[time - 1]) * (1 - refrac)
        org[1, 4] = float(orgin.loc[time - 1]) * (refrac)
        mini[1] = float(minin.loc[time - 1]) * rtemin(relelv[time - 1])
        tpore[1] = 1
        pctpore[1] = float(h2oin.loc[0])
        pore[1] = ((float(h2oin.loc[time - 1]) / (1 - float(h2oin.loc[time - 1]))) * (
                    (org[1, 1] / orgden) + (mini[1] / minden)))

        # ******************************************************************
        # * the following section is where the "yearly" calculations       *
        # * take place.  It combines all of the other calculation sections *
        # * from earlier versions of the model (7/10/93).                  *
        # ******************************************************************
        #
        # this section calculates the volume of each section
        # based on the mass of organic matter, mineral matter, and
        # water.  It will also use a compaction subfunction in
        # the future. Compaction will be a function of the
        # mass that is on top of the current section.
        #
        #
        # this is where the new roots and rhizomes are put into
        # the sediment.  rtprod is a subroutine/function
        # that will determine root production based on depth/time
        #
        # this is also the decomposition section.  Again decomp is
        # a subroutine based on depth/time.
        #

        for t2 in range(1, time + 1):

            istep = 10
            dt = 1 / istep
            for ie in range(1, istep + 1):
                totorg[t2] = org[t2, 1] + org[t2, 2] + org[t2, 3] + org[t2, 4]
                massabv[t2] = massabv[t2 - 1] + totorg[t2 - 1] + mini[t2 - 1] + pore[t2 - 1]
                if depth[t2 - 1] == 0:
                    densabv[t2] = 0
                else:
                    densabv[t2] = massabv[t2] / depth[t2 - 1]
                orgvol[t2] = totorg[t2] / orgden
                minvol[t2] = mini[t2] / minden
                if t2 <= 1:
                    pctpore[t2] = pctpore[t2]
                else:
                    dum1 = float(h2oin.loc[time - 1])
                    dum2 = float(porelim.loc[time - 1])
                    tpore[t2] = tpore[t2] - (tpore[t2] - tpore[t2] * poresp(densabv[t2], k1)) * dt
                    pctpore[t2] = dum2 + (dum1 - dum2) * tpore[t2]

                pore[t2] = ((pctpore[t2] / (1 - pctpore[t2])) * (orgvol[t2] + minvol[t2]))

                # c the line above is for running the model without compaction
                # c pore space is constant for all sections.
                # c
                totvol[t2] = orgvol[t2] + minvol[t2] + pore[t2]
                depth[t2] = depth[t2 - 1] + totvol[t2]
                porg[t2] = totorg[t2] / (totorg[t2] + mini[t2])
                bulkd[t2] = (totorg[t2] + mini[t2]) / totvol[t2]
                orgbd[t2] = totorg[t2] / totvol[t2]
                minbd[t2] = mini[t2] / totvol[t2]
                org[t2, 1] = org[t2, 1] + ((dt * (rtprod(depth[t2]) * totvol[t2])) * (1 - refrac)) \
                             - ((dt * (((decmp1(depth[t2], kdec1)) * org[t2, 1]))))
                org[t2, 2] = org[t2, 2] - \
                             (dt * (((decmp2(depth[t2])) * org[t2, 2])))
                org[t2, 3] = org[t2, 3] \
                             - (dt * ((decmp3(depth[t2])) * org[t2, 3]))
                org[t2, 4] = org[t2, 4] + ((dt * (rtprod(depth[t2]) * totvol[t2])) * refrac)

        # c
        # c commenting out the 3 "&" lines above, cuts out decompostion
        # c
        # c           write(29,*) decmp3(depth(t2)), depth(t2)
        # c

        # c
        # c **********************************************************************
        # c * this section calculates the relative elevation of the marsh at the *
        # c * end of the year                                                    *
        # c **********************************************************************
        # c

        relelv[time] = intelv + depth[time] - (slr * time) - (subsid * time)

    sksxx = pd.DataFrame({"totorg": totorg[1:],
                          "min": mini[1:],
                          "pore": pore[1:],
                          "totvol": totvol[1:],
                          "orgvol": orgvol[1:],
                          "minvol": minvol[1:],
                          "porg": porg[1:],
                          "bulkd": bulkd[1:],
                          "depth": depth[1:],
                          "massabv": massabv[1:],
                          "densabv": densabv[1:],
                          "relelv": np.flip(relelv, 0)[1:],
                          "time": np.array(range(1, endtim + 1))})
    sksxx = sksxx.reindex(
        columns=['totorg', "min", "pore", "totvol", "orgvol", "minvol", "porg", "bulkd", "depth", "massabv", "densabv",
                 "relelv", "time"])

    sksxx.to_csv("sksxx.csv", index=False)

    return sksxx


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
os.chdir(r"C:\Projects\5630\12072021")

#Numpy arrays
np_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\Numpy"

#Shapefiles
shp_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\GW Model\25ft\20211207\Output"

#Rasters
ras_dir=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\raster\20211207\Output"

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

#Let's export shapefile of grid

ml.modelgrid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)
#ml.dis.export(os.path.join(shp_dir, "Grid.shp"))
grid = ml.modelgrid
#grid.set_coord_info(xoff=6249820, yoff=2165621, epsg=2227)
nrg=ml.nrow
ncg=ml.ncol

subsidence=np.zeros(shape=(nrg,ncg))

class_length=int(nrg/n_cpu)

#Let's get active cells
bas = flopy.modflow.ModflowBas.load('MF_inputs/Bacon.bas', ml)

ibound=np.array(bas.ibound[0][:])
Active_cells=np.where(ibound!=0)

#Let's get masks of wetland and rice scenario
#bas_wlr = flopy.modflow.ModflowBas.load('MF_inputs/BaconWLR_year0.bas', ml)

#Let's get mask of levees
levees=np.where(ml.lpf.hk[0][:]==min(np.unique(ml.lpf.hk[0][:])))

#Let's get indexes of toedrains
toedrains_in=pd.read_csv(r"\\hydro-nas\Team\Projects\5630_DSC\GIS\vector\GW Model\25ft\Drains\ToeDrains_Index_Flopy.csv")
toedrains_in=toedrains_in.drop(['k','elev', 'cond', 'iface', 'top', 'PT_bot', 'TM_bot','h_PT', 'h_TM', 'h_SP', 'Grad'], axis=1)

#Let's convert to recarray
toedrains_in_rec=toedrains_in.to_records

#Let´s get mask of rice
rice_df=pd.read_csv(r"C:\Projects\5630\12072021\Rice.csv")
rice_mask=list(zip(rice_df.row,rice_df.column_lef))

#Let´s get mask of wetlands
wetland_df=pd.read_csv(r"C:\Projects\5630\12072021\Wetlands.csv")
wetland_mask=list(zip(wetland_df.row,wetland_df.column_lef))


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

#Let´s run SEDCALC
sedcalc_ts=sedcalc(endtim=End_Year+1,)

#3. Let's loop through years now


for year in range(Start_Year,End_Year+1):
    
    # 3.1 Initiate MODFLOW for year 1
    if year==Start_Year:
        ml = flopy.modflow.Modflow.load('MF_inputs/Bacon.nam')
        ml.write_input()
        subprocess.check_output(["mf2005", "Bacon_fix.nam"])
        #Let's set rice and wetland to constant head
        bas.ibound[0][rice_mask] = -1
        bas.ibound[0][wetland_mask] = -1
        bas.write_file()

    #Let´s add constant head cells for wetlands and rice

    #3.2 Let's run SUBCALC
    
    #Peat thickness
    PT_thck=ml.dis.thickness[0]
    
    #Depth to groundwater
    h = flopy.utils.HeadFile("Bacon.hds", model=ml)
    heads=h.get_data()
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
    drns_pd.to_csv(os.path.join(np_dir,"DRNS"+str(year)+".csv"),index=False)
    #drns.tofile(os.path.join(np_dir,"DRNS"+str(year)+".csv"),sep=",")
    print(year,ti-t0)

#Let's export average subsidence dataframe
avg_sub_df.to_csv(os.path.join(np_dir,"Avg_Sub.csv"),index=False)
         
        
        
