import pandas as pd
import os
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
import pyhf



#Path to csv files
csv_path=r"\\hydro-nas\Team\Projects\5630_DSC\GIS\CSV\20221122"

#Names of series
names=["BAU","WLR"]

out_path=r"\\hydro-nas\Team\Projects\5630_DSC\Paper\2022_12\Gradients"

min, max=pyhf.utils.ecdf(csv_path,names,out_path,video=True,calc_x=True,End_Year=2034)

#key of column with values
valkey="Grad"
lb_prefix="LB"
ub_prefix="UB"
exp_prefix="Base"
suffix='TOEDRNS'
year=2034
nbins=1000
bands=True
calc_x=True
Start_Year=2018
End_Year=2070
min_x=-1
max_x=5

if calc_x:
    all_grads = np.empty(0)
    for year in range(Start_Year, End_Year + 1):
        for serie in names:
            exp_df = pd.read_csv(os.path.join(csv_path, exp_prefix + "_" + serie + "_" + suffix + str(year) + ".csv"))
            if bands:
                all_grads = np.concatenate((all_grads, lb_df.Grad.values, ub_df.Grad.values, exp_df.Grad.values))
            else:
                all_grads = np.concatenate((all_grads, exp_df.Grad.values))
    x = np.sort(np.unique(all_grads))
else:
    min = min_x
    max = max_x
    # Let's create plot points
    x = np.linspace(min, max, nbins)