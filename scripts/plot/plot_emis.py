#Plot emissions all components in the CICERO SCM
#If unit is 'X', no emissions is provided, but component in output file.

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams['font.size'] = 4

#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_em=pd.read_csv(outdir+'/output_em.txt', sep='\t', index_col=0)
#print(df_em)

#Read components, to get the units:
inputdir = '/div/no-backup/users/ragnhibs/ciceroscm/tests/test-data/'
df_comp =pd.read_csv(inputdir + 'gases_v1RCMIP.txt', sep='\t', index_col=0)
#print(df_comp.loc['CO2']['EM_UNIT'])

antcomp = len(df_em.columns)
print(antcomp)
#print(len(df_comp.index))


#Plot first 16 components:
complist = np.arange(0,16)
fig, axs = plt.subplots(nrows=4, ncols=4,sharex=True,figsize=(12,20))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Emissions 1 of 3')
for i,c in enumerate(complist):
    print(i)
    comp = df_em.columns[c]
    print(comp)
    df_em[comp].plot(ylabel='Emissions ['+df_comp.loc[comp]['EM_UNIT']+']',ax=axs[i],label=scen)
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].set_ylim(bottom=0)

    
#Plot next 16 components:
complist = np.arange(16,32)
fig, axs = plt.subplots(nrows=4, ncols=4,sharex=True,figsize=(12,20))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Emissions 2 of 3')
for i,c in enumerate(complist):
    print(i)
    comp = df_em.columns[c]
    print(comp)
    df_em[comp].plot(ylabel='Emissions ['+df_comp.loc[comp]['EM_UNIT']+']',ax=axs[i],label=scen)
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].set_ylim(bottom=0)
    
#Plot rest of components:
complist = np.arange(32,antcomp)
fig, axs = plt.subplots(nrows=4, ncols=4,sharex=True,figsize=(12,20))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Emissions 3 of 3')
for i,c in enumerate(complist):
    print(i)
    comp = df_em.columns[c]
    print(comp)
    df_em[comp].plot(ylabel='Emissions ['+df_comp.loc[comp]['EM_UNIT']+']',ax=axs[i],label=scen)
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].set_ylim(bottom=0)

plt.show()
exit()
