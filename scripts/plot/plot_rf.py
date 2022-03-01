#Plot radiative forcing all components in the CICERO SCM
#Not all components have radiative forcing output.

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams['font.size'] = 4

#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_rf=pd.read_csv(outdir+'/output_forc.txt', sep='\t', index_col=0)
#print(df_rf)

#Read components, to get the units:
inputdir = '/div/no-backup/users/ragnhibs/ciceroscm/tests/test-data/'
df_comp =pd.read_csv(inputdir + 'gases_v1RCMIP.txt', sep='\t', index_col=0)
print(df_comp.loc['CO2'])

antcomp = len(df_rf.columns)
print(antcomp)
#print(len(df_comp.index))


#Plot first 16 components:
complist = np.arange(0,16)
fig, axs = plt.subplots(nrows=4, ncols=4,sharex=True,figsize=(12,20))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Radiative Forcing 1 of 3')
for i,c in enumerate(complist):
    print(i)
    comp = df_rf.columns[c]
    print(comp)
    df_rf[comp].plot(ylabel='RF [Wm$^{-2}$ ]',ax=axs[i],label=scen)
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].axhline(y=0,color='k',linestyle=':',linewidth=0.5)
    
#Plot next 16 components:
complist = np.arange(16,32)
fig, axs = plt.subplots(nrows=4, ncols=4,sharex=True,figsize=(12,20))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Radiative Forcing 2 of 3')
for i,c in enumerate(complist):
    print(i)
    comp = df_rf.columns[c]
    print(comp)
    df_rf[comp].plot(ylabel='RF [Wm$^{-2}$]',ax=axs[i],label=scen)
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].axhline(y=0,color='k',linestyle=':',linewidth=0.5)
#Plot rest of components:
complist = np.arange(32,antcomp)
fig, axs = plt.subplots(nrows=4, ncols=4,sharex=True,figsize=(12,20))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Radiative Forcing 3 of 3')
for i,c in enumerate(complist):
    print(i)
    comp = df_rf.columns[c]
    print(comp)
    df_rf[comp].plot(ylabel='RF [Wm$^{-2}$]',ax=axs[i],label=scen)
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].axhline(y=0,color='k',linestyle=':',linewidth=0.5)

plt.show()
exit()
