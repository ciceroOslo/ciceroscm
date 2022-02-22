#Compare CSCM results to
#von Schuckmann et al. Heat stored in the Earth system:
#Where does the energy #go? Earth Syst. Sci. Data, 12,2013-2041,
#10.5194/essd-12-2013-2020, 2020.

#The python environment for cicero scm did not work for nc-file.
#Read file prepared for UTRICS.

#Can add uncertainties and other timeseries as well.

#For comparing with observations, important to add volcanoes as well.

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr


#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_ohc=pd.read_csv(outdir+'/output_ohc.txt', sep='\t', index_col=0)
print(df_ohc)

#Did not work with the CICERO SCM environment
#filename =  'data_compare/GCOS_all_heat_content_1960-2018_ZJ_v22062020.nc'
#gcos = xr.open_dataset(filename)
#print(gcos)

baseyears = [1961,1990]
baseyears_str = str(baseyears[0]) +'-' + str(baseyears[1])


df_gcos = pd.read_csv('/div/qbo/utrics/Observations/OHC/GCOS20_vonSchuckmann/GCOS_ohc.csv',index_col=0)
df_gcos = df_gcos*0.1 #ZJ -> 10^22 J
df_gcos = df_gcos - df_gcos.loc[baseyears[0]:baseyears[1]].mean()
print(df_gcos)


fig, axs = plt.subplots(nrows=1, ncols=2,sharex=True,figsize=(12,6))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Ocean heat content')


comp = 'OHC700'
df_ohc['anomaly 0-700'] = df_ohc[comp] - df_ohc[comp].loc[baseyears[0]:baseyears[1]].mean()
df_ohc['anomaly 0-700'].plot(ylabel='OHC [10^22 J]',ax=axs[0],label=scen)
df_gcos['ohc_0_700m'].plot(ax=axs[0],color='black',label='GCOS')

axs[0].set_title(comp)


comp = 'OHCTOT'
df_ohc['anomaly TOT'] = df_ohc[comp] - df_ohc[comp].loc[baseyears[0]:baseyears[1]].mean()
df_ohc['anomaly TOT'].plot(ylabel='OHC [10^22 J]',ax=axs[1],label=scen)
df_gcos['tot'] = df_gcos['ohc_0_700m'] + df_gcos['ohc_700_2000m']+ df_gcos['ohc_below_2000m']
df_gcos['tot'].plot(ax=axs[1],color='black',label='GCOS')

axs[1].set_title(comp)

axs[0].legend()
axs[1].legend()

plt.show()



exit()
