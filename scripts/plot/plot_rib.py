#Plot radiative imbalance
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_rib=pd.read_csv(outdir+'/output_rib.txt', sep='\t', index_col=0)
print(df_rib)


fig, axs = plt.subplots(nrows=1, ncols=2,sharex=True,figsize=(12,6))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, RIB')


comp = 'RIB_glob'
df_rib[comp].plot(ylabel='check',ax=axs[0],label=scen)
axs[0].set_title(comp)



df_rib['RIB_glob'].plot(ylabel='check ',ax=axs[1],label=scen)
df_rib['RIB_N'].plot(ax=axs[1],linestyle=':',label='RIB_N')
df_rib['RIB_S'].plot(ax=axs[1],linestyle='--',label='RIB_S')
axs[1].set_title(comp)

axs[0].legend()
axs[1].legend()

plt.show()



exit()
