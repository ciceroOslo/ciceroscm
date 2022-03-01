#Plot radiative forcing  in the CICERO SCM
#and compare to IPCC AR6 forcing timeseries.
#https://github.com/chrisroadmap/ar6/blob/main/data_output/AR6_ERF_1750-2019.csv
#https://raw.githubusercontent.com/chrisroadmap/ar6/main/data_output/AR6_ERF_1750-2019.csv
#Can add 90% CI as well.


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams['font.size'] = 4

#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_rf=pd.read_csv(outdir+'/output_forc.txt', sep='\t', index_col=0)
#print(df_rf)


#Read ipcc AR6 forcing timeseries
df_rf_ar6 = pd.read_csv('data_compare/AR6_ERF_1750-2019.csv',sep=',',header=0,index_col=0)

#print(df_rf_ar6.columns)
#Components given in IPCC
#'co2', 'ch4', 'n2o', 'other_wmghg', 'o3', 'h2o_stratospheric',
#       'contrails', 'aerosol-radiation_interactions',
#       'aerosol-cloud_interactions', 'bc_on_snow', 'land_use', 'volcanic',
#       'solar', 'nonco2_wmghg', 'aerosol', 'chapter2_other_anthro',
#       'total_anthropogenic', 'total_natural', 'total'

complist_ar6 = ['co2', 'ch4', 'n2o', 'other_wmghg',
                'o3','h2o_stratospheric','land_use','aerosol',
                'total_anthropogenic']


fig, axs = plt.subplots(nrows=3, ncols=3,sharex=True,figsize=(12,20))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation, Radiative Forcing')
for i,comp in enumerate(complist_ar6):
    print(i)
    print(comp)
    df_rf_ar6[comp].plot(ylabel='RF [Wm$^{-2}$ ]',
                         ax=axs[i],color='black',label='IPCC AR6')
    if comp == 'co2':
        df_rf['CO2'].plot(ax=axs[i],label=scen)
    elif comp == 'ch4':
        df_rf['CH4'].plot(ax=axs[i],label=scen)
    elif comp == 'n2o':
        df_rf['N2O'].plot(ax=axs[i],label=scen)
    elif comp == 'other_wmghg':
        complist_other=['CFC-11', 'CFC-12', 'CFC-113', 'CFC-114',
                        'CFC-115', 'CH3Br', 'CCl4', 'CH3CCl3', 'HCFC-22', 'HCFC-141b',
                        'HCFC-142b', 
                        'C2F6', 'C6F14', 'CF4', 'SF6','HCFC-123',          
                        'H-1211', 'H-1301', 'H-2402','HFC125','HFC134a',
                        'HFC143a','HFC227ea', 'HFC23', 'HFC245fa','HFC32',
                        'HFC4310mee']
        df_rf['other_wmghg'] = df_rf[complist_other].sum(axis=1, skipna=True)
        #print(df_rf['other_wmghg'])
        df_rf['other_wmghg'].plot(ax=axs[i],label=scen)
    elif comp == 'o3':
        print(df_rf.columns)
        df_rf['o3'] = df_rf['TROP_O3']+df_rf['STRAT_O3'] 
        df_rf['o3'].plot(ax=axs[i],label=scen)
    elif comp == 'h2o_stratospheric':
        df_rf['STRAT_H2O'].plot(ax=axs[i],label=scen)
    elif comp == 'land_use':
        df_rf['LANDUSE'].plot(ax=axs[i],label=scen)
    elif comp == 'aerosol':
        complist_aerosols = ['SO2','SO4_IND','BMB_AEROS_BC',
                             'BMB_AEROS_OC','BMB_AEROS', 'BC', 'OC']
        df_rf['aerosol'] = df_rf[complist_aerosols].sum(axis=1, skipna=True)
        df_rf['aerosol'].plot(ax=axs[i],label=scen)
    elif comp == 'total_anthropogenic':
        df_rf['Total_forcing'].plot(ax=axs[i],label=scen)
    axs[i].set_title(comp)
    axs[i].legend()
    axs[i].axhline(y=0,color='k',linestyle=':',linewidth=0.5)
    

    
fig, axs = plt.subplots(nrows=1, ncols=1,sharex=True,figsize=(6,6))
fig.suptitle('CICERO SCM simulation, Aerosol Radiative Forcing')


comp = 'aerosol-radiation_interactions'
df_rf_ar6[comp].plot(ylabel='RF [Wm$^{-2}$ ]',
                     ax=axs,color='green',label='IPCC AR6 '+comp)

comp = 'aerosol-cloud_interactions'
df_rf_ar6[comp].plot(ylabel='RF [Wm$^{-2}$ ]',
                     ax=axs,color='blue',label='IPCC AR6 '+comp)

complist_aerosols_direct = ['SO2','BMB_AEROS_BC',
                            'BMB_AEROS_OC','BMB_AEROS', 'BC', 'OC']
df_rf['aerosol_direct'] = df_rf[complist_aerosols_direct].sum(axis=1, skipna=True)
df_rf['aerosol_direct'].plot(ax=axs,color='purple',label='SCM ' + scen + ' sum aerosols direct')
df_rf['SO4_IND'].plot(ax=axs,color='orange',label='SCM ' + scen + ' SO4 indirect')



axs.legend()

plt.show()
exit()
