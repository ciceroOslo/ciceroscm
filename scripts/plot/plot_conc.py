#Plot concentration all components in the CICERO SCM
#If unit is '-', no concentration is provided, but component in output file.
#To do: Trop O3, concentration output given, unit '-'. This is DobsenUnit. 

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams['font.size'] = 4

#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_conc=pd.read_csv(outdir+'/output_conc.txt', sep='\t', index_col=0)
#print(df_conc)

#Read components, to get the units:
inputdir = '/div/no-backup/users/ragnhibs/ciceroscm/tests/test-data/'
df_comp =pd.read_csv(inputdir + 'gases_v1RCMIP.txt', sep='\t', index_col=0)
#print(df_comp.loc['CO2']['EM_UNIT'])

antcomp = len(df_conc.columns)
print(antcomp)
#print(len(df_comp.index))


#Plot first 16 components:
complist = np.arange(0,16)
fig, axs = plt.subplots(nrows=4, ncols=4,sharex=True,figsize=(12,20))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, concentrations 1 of 3')
for i,c in enumerate(complist):
    print(i)
    comp = df_conc.columns[c]
    print(comp)
    df_conc[comp].plot(ylabel='Conc. ['+df_comp.loc[comp]['CONC_UNIT']+']',
                       ax=axs[i],label=scen)
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
    comp = df_conc.columns[c]
    print(comp)
    df_conc[comp].plot(ylabel='Conc. ['+df_comp.loc[comp]['CONC_UNIT']+']',ax=axs[i],label=scen)
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
    comp = df_conc.columns[c]
    print(comp)
    df_conc[comp].plot(ylabel='Conc. ['+df_comp.loc[comp]['CONC_UNIT']+']',
                       ax=axs[i],label=scen)
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].set_ylim(bottom=0)

plt.show()
exit()
