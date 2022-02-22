#Plot concentration and compare to CMIP6 input 
#Meinshausen, M., Vogel, E., Nauels, A., Lorbacher, K., Meinshausen,
#N., Etheridge, D. M., Fraser, P. J., Montzka, S. A., Rayner,
#P. J., Trudinger, C. M., Krummel, P. B., Beyerle, U., Canadell,
#J. G., Daniel, J. S., Enting, I. G., Law, R. M., Lunder, C. R.,
#O'Doherty, S., Prinn, R. G., Reimann, S., Rubino, M., Velders,
#G. J. M., Vollmer, M. K., Wang, R. H. J., and Weiss,
#R.: Historical greenhouse gas concentrations for climate
#modelling (CMIP6), Geosci. Model Dev., 10,2057-2116,
#10.5194/gmd-10-2057-2017, 2017.

#Several components here not included in the CSCM
#Only compared the components with identical names.

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

print(df_comp.index)

#Selected components are those with identical component names as in CMIP6 file.
complist=['CFC-11', 'CFC-12', 'CFC-113', 'CFC-114',
          'CFC-115', 'CH3Br', 'CCl4', 'CH3CCl3', 'HCFC-22', 'HCFC-141b',
          'HCFC-142b', 
          'C2F6', 'C6F14', 'CF4', 'SF6']

compskiped = ['HCFC-123',          
              'H-1211', 'H-1301', 'H-2402','HFC125','HFC134a',
              'HFC143a','HFC227ea', 'HFC23', 'HFC245fa','HFC32',
              'HFC4310mee']

cmip6_conc = pd.read_excel('data_compare/Supplementary_Table_UoM_GHGConcentrations-1-1-0_annualmeans_v23March2017.xls',
                           sheet_name='historical-annualmean-Global',
                           header=21,index_col=0)
cmip6_conc.index.name='Years'
print(cmip6_conc)


fig, axs = plt.subplots(nrows=4, ncols=4,sharex=True,figsize=(12,16))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, concentrations')
for i,comp in enumerate(complist):
    df_conc[comp].plot(ylabel= comp +' ['+df_comp.loc[comp]['CONC_UNIT']+']',
                       ax=axs[i],label=scen)
    
    cmip6_conc[comp].plot(ax=axs[i],color='black',label='CMIP6input')
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].set_ylim(bottom=0)

df_conc[compskiped].plot(ylabel='ppt',ax=axs[15])    
axs[15].set_title('Rest of CSCM components')
axs[15].legend(loc='lower left')
plt.show()
