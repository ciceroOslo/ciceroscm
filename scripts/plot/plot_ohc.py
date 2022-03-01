import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_ohc=pd.read_csv(outdir+'/output_ohc.txt', sep='\t', index_col=0)
print(df_ohc)


fig, axs = plt.subplots(nrows=1, ncols=2,sharex=True,figsize=(12,6))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Ocean heat content')


comp = 'OHC700'
df_ohc[comp].plot(ylabel='OHC [$10^{22}$ J]',ax=axs[0],label=scen)
axs[0].set_title(comp)


comp = 'OHCTOT'
df_ohc[comp].plot(ylabel='OHC [$10^{22}$ J]',ax=axs[1],label=scen)
axs[1].set_title(comp)

axs[0].legend()
axs[1].legend()

plt.show()



exit()
