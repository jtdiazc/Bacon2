import rasterio
import pandas as pd
import os
import numpy as np


fom=np.load(r"C:\Projects\5630\20221122\Base\fom_0.npy")
fomund=np.load(r"C:\Projects\5630\20221122\Base\fomund_0.npy")
dtwm=np.load(r"C:\Projects\5630\20221122\Base\BAU\DTW_0.npy")
peatdepth=np.load(r"C:\Projects\5630\20221122\Base\BAU\PeatDepth_0.npy")
startyr=2018
sim_id='0'
st=None
km=None
firstkm=None
vmax=None
bd=None
bdund=None
bdund2=1.3
fomund2=0
massmin=None

tempeff=False
tempc="model"
rcp="4.5"
initelev=0
output="table"
cols="all"
silent=True

endyr = startyr
### Simulation depth and inital elevation ###
simdpth = dtwm * 100
#    elev = initelev

##### Calcualte BD for unsat and sat zones based on OM if not provided #####
if bd is None:
    bd = -0.354 * np.log(fom * 100) + 1.7561

if bdund is None:
    bdund = -0.21 * np.log(fomund * 100) + 1.01

##### Calculate total substrate from OM if not provided #####
if st is None:
    st = (fom * 0.5) * bd * simdpth

##### Calculate mineral mass if not provided #####
if massmin is None:
    massmin = (1 - fom) * bd * simdpth

### Inital check if water table is below the peat ###


st[np.where(simdpth >= peatdepth)] = fom[np.where(simdpth >= peatdepth)] * 0.5 * bd[
    np.where(simdpth >= peatdepth)] * peatdepth[np.where(simdpth >= peatdepth)]
#        fomund = fomund2
fomund[np.where(simdpth >= peatdepth)] = fomund2
#        bdund = bdund2
bdund[np.where(simdpth >= peatdepth)] = bdund2
#        massmin = ((1-fom)*bd)*peatdepth
massmin[np.where(simdpth >= peatdepth)] = ((1 - fom[np.where(simdpth >= peatdepth)]) * bd[
    np.where(simdpth >= peatdepth)]) * peatdepth[np.where(simdpth >= peatdepth)]

#    else:
#        peatdepth = np.nan
peatdepth[np.where(simdpth < peatdepth)] = np.nan



## Calculate km value
if km is None:
    km = (-4.2732 * st) + 51.709
    #        firstkm = False # flag for resetting km in loop
    firstkm = np.full(km.shape, False)

#    if km > 0:
#        firstkm = True
firstkm[np.where(km > 0)] = True

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
    vmax = preexp * np.exp(-activ / (0.0083144 * tempk))

##### Begin SUBCALC loop #####
#    for year in years:

## output table row ##
idx = 0

year = startyr


delcflux = 0

#    output_table.at[idx,'delcflux'] = delcflux

cflux = np.zeros(km.shape)
### Calculate oxidation ###

ind_dum = np.where(km <= 0)
cflux[ind_dum] = vmax + delcflux


ind_dum = np.where(km > 0)
cflux[ind_dum] = ((vmax * st[ind_dum]) / (km[ind_dum] + st[ind_dum])) + delcflux

omflux = cflux * 2
delelvc = omflux / (bd * fom)


### Calculate consolidation ###
delp = delelvc
delelvcons = conconst * condepth * delp

### combine sources of subsidence ###
tdelev = delelvc + delelvcons
#    elev = elev - tdelev

if peatdepth is not None:


    peatleft = peatdepth - tdelev

    ind_dum = np.where(simdpth + tdelev >= peatdepth)
    fomund[ind_dum] = fomund2

    bdund[ind_dum] = bdund2


else:
    peatleft = simdpth

### Update variables ###
volloss = tdelev
massomlft = (st * 2) - omflux
massom = massomlft + (fomund * bdund * volloss)
massmin = massmin + ((1 - fomund) * bdund * volloss)

#Let's deal with nans
cond=np.isnan(massom / (massom + massmin))
fom[~cond] = np.maximum(0, massom[~cond] / (massom[~cond] + massmin[~cond]))

#    if simdpth <= peatleft:
#        bd = (massom + massmin) / simdpth
ind_dum = np.where(simdpth <= peatleft)
bd[ind_dum] = (massom[ind_dum] + massmin[ind_dum]) / simdpth[ind_dum]

#    else:
#        bd = (massom + massmin) / peatleft
ind_dum = np.where(simdpth > peatleft)
bd[ind_dum] = (massom[ind_dum] + massmin[ind_dum]) / peatleft[ind_dum]

st = massom / 2
kmnew = (-4.2732 * st) + 51.709  # used to "turn off" clfux = vmax

#   if kmnew > 0 and firstkm==False:

#        km = kmnew
ind_dum = np.where((kmnew > 0) & (firstkm == False))
km[ind_dum] = kmnew[ind_dum]
#        firstkm = True
firstkm[ind_dum] = True

Output = {'km': km,
          'vmax': vmax,
          'fomund': fomund,
          'bdund': bdund,
          'cflux': cflux,
          'delelvc': delelvc,
          'delelvcons': delelvcons,
          'tdelev': tdelev,
          'volloss': volloss,
          'massomlft': massomlft,
          'massom': massom,
          'massmin': massmin,
          'fom': fom,
          'bd': bd,
          'st': st,
          'firstkm': firstkm}
