import pyhf
import os
import pandas as pd
import flopy
import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
#from moviepy.editor import VideoClipPM
from moviepy.video.io.bindings import mplfig_to_npimage


#Path to csv files
csv_path=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV\20221122"

#Names of series
names=["BAU","WLR"]

out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_12\Gradients"



#key of column with values
valkey="Grad"
lb_prefix="LB"
ub_prefix="UB"
exp_prefix="Base"
suffix='TOEDRNS'

serie=names[0]
nbins=10000

leg_x = 1
leg_y = 0.8

x=np.linspace(-1,5,nbins)

year=2070

fig, ax = plt.subplots()

kernel='linear'
bandwidth=0.01
line_width=0.1
#BAU
serie="BAU"
exp_df = pd.read_csv(os.path.join(csv_path,exp_prefix+"_"+serie+"_"+suffix+str(2070)+".csv"))
kde_exp = KernelDensity(kernel=kernel, bandwidth=bandwidth).fit(exp_df[valkey].values.reshape(-1,1))
log_dens_exp = kde_exp.score_samples(x.reshape(-1,1))
#ax.plot(x, np.exp(log_dens_exp),label=serie,linewidth=line_width)
lb_df = pd.read_csv(os.path.join(csv_path,lb_prefix+"_"+serie+"_"+suffix+str(year)+".csv"))
ub_df = pd.read_csv(os.path.join(csv_path,ub_prefix+"_"+serie+"_"+suffix+str(year)+".csv"))
# Let's calculate kernel density of lb
kde_lb = KernelDensity(kernel=kernel, bandwidth=bandwidth).fit(lb_df[valkey].values.reshape(-1, 1))
exp_dens=np.exp(log_dens_exp)

# Log densities
log_dens_lb = kde_lb.score_samples(x.reshape(-1, 1))
lb_dens=np.exp(log_dens_lb)
# Let's calculate kernel density of ub
kde_ub = KernelDensity(kernel=kernel, bandwidth=bandwidth).fit(ub_df[valkey].values.reshape(-1, 1))
# Log densities
log_dens_ub = kde_ub.score_samples(x.reshape(-1, 1))
ub_dens=np.exp(log_dens_ub)
ax.fill_between(x, lb_dens, ub_dens)

band_width=log_dens_ub-log_dens_lb

band=pd.DataFrame({'x':x,'lb_dens':lb_dens,'ub_dens':ub_dens})
band['width']=band.ub_dens-band.lb_dens
band['lb_plot']=np.minimum(band.ub_dens,band.lb_dens)
band['ub_plot']=np.maximum(band.ub_dens,band.lb_dens)
band=band[band.width!=0]

ax.fill_between(band.x, band.lb_plot, band.ub_plot)

#WLR
serie="WLR"
exp_df = pd.read_csv(os.path.join(csv_path,exp_prefix+"_"+serie+"_"+suffix+str(year)+".csv"))
kde_exp = KernelDensity(kernel=kernel, bandwidth=bandwidth).fit(exp_df[valkey].values.reshape(-1,1))
log_dens_exp = kde_exp.score_samples(x.reshape(-1,1))
ax.plot(x, np.exp(log_dens_exp),label=serie,linewidth=line_width)
ax.legend(bbox_to_anchor=(leg_x, leg_y))

plt.savefig(os.path.join(out_path, "density_" + str(year) + ".png"))
plt.close()