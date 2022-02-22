#Can add other temperature series than HadCRUT.

from io import StringIO
import pandas as pd
import requests
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.size'] = 8

#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test/'
scen = 'test'
df_temp=pd.read_csv(outdir+'/output_temp.txt', sep='\t', index_col=0)
print(df_temp)

antcomp = len(df_temp.columns)
print(antcomp)
#print(len(df_comp.index))

#Read HADCRUT temp
headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/56.0.2924.76 Safari/537.36'}

url_hadcrut = "https://www.metoffice.gov.uk/hadobs/hadcrut5/data/current/analysis/diagnostics/HadCRUT.5.0.1.0.analysis.summary_series.global.annual.csv"
s=requests.get(url_hadcrut, headers= headers).text
hadcrut=pd.read_csv(StringIO(s),index_col=0)

#Possible to add several other dataseries as well.

#Same base period as HadCRUT
baseyears = [1961,1990]
baseyears_str = str(baseyears[0]) +'-' + str(baseyears[1])



fig, axs = plt.subplots(nrows=1, ncols=1,sharex=True,figsize=(8,6))

fig.suptitle('CICERO SCM simulation, Temperature')

df_temp['anomaly'] = df_temp['dT_glob']- df_temp['dT_glob'].loc[baseyears[0]:baseyears[1]].mean()

df_temp['anomaly'].plot(ax=axs,label=scen)
hadcrut['Anomaly (deg C)'].plot(ax=axs,color='black',label='HadCRUT')

axs.set_ylabel('GMST relative to ' + baseyears_str + '[$^\circ$C]')

axs.legend()

plt.show()
