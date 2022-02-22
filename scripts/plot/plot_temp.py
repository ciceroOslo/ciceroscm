#To be fixed in scm-code: dT_SHsea

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams['font.size'] = 4

#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_temp=pd.read_csv(outdir+'/output_temp.txt', sep='\t', index_col=0)
print(df_temp)

antcomp = len(df_temp.columns)
print(antcomp)
#print(len(df_comp.index))


fig, axs = plt.subplots(nrows=1, ncols=3,sharex=True,figsize=(12,6))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Temperature')


i=0
comp = 'dT_glob'
df_temp[comp].plot(ylabel='Temp. [C$^\circ$ ]',ax=axs[i],label=scen)

i=1
df_temp['dT_glob'].plot(ylabel='Temp. [C$^\circ$ ]',ax=axs[i],label='dT_glob')
df_temp['dT_glob_air'].plot(ax=axs[i],label='dT_glob_air')
df_temp['dT_glob_sea'].plot(ax=axs[i],label='dT_glob_sea')

i=2
df_temp['dT_glob'].plot(ylabel='Temp. [C$^\circ$ ]',ax=axs[i],label='dT_glob')
df_temp['dT_glob_air'].plot(ax=axs[i],label='dT_glob_air')
df_temp['dT_glob_sea'].plot(ax=axs[i],label='dT_glob_sea')

df_temp['dT_NH'].plot(linestyle=':',ax=axs[i],label='dT_NH')
df_temp['dT_NH_air'].plot(ax=axs[i],linestyle=':',label='dT_NH_air')
df_temp['dT_NH_sea'].plot(ax=axs[i],linestyle=':',label='dT_NH_sea')

df_temp['dT_SH'].plot(linestyle='--',ax=axs[i],label='dT_SH')
df_temp['dT_SH_air'].plot(ax=axs[i],linestyle='--',label='dT_SH_air')
df_temp['dT_SHsea'].plot(ax=axs[i],linestyle='--',label='dT_SHsea')

axs[0].legend()
axs[1].legend()
axs[2].legend()

axs[0].axhline(y=0,color='k',linestyle=':',linewidth=0.5)


fig, axs = plt.subplots(nrows=1, ncols=2,sharex=True,figsize=(12,6))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, SLR (!!!!!need to be revised-checked!!!!)')
df_temp['dSL(m)'].plot(ylabel='Sea level rise [m]',ax=axs[0],label=scen)

df_temp['dSL(m)'].plot(ylabel='Sea level rise [m]',ax=axs[1],label='dSL(m)')
df_temp['dSL_thermal(m)'].plot(ylabel='Sea level rise [m]',
                               ax=axs[1],label='dSL_thermal(m)')
df_temp['dSL_ice(m)'].plot(ylabel='Sea level rise [m]',
                           ax=axs[1],label='dSL_ice(m)')

axs[0].legend()
axs[1].legend()

plt.show()
