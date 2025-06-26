
import sys
import re
import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
import pandas as pd
import pandas.testing as pdt
import warnings
import logging
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import concurrent
from ciceroscm.parallel._configdistro import _ConfigDistro
from ciceroscm.parallel.calibrator import Calibrator
from ciceroscm.parallel.distributionrun import DistributionRun


sys.path.insert(0,os.path.join(os.getcwd(), '../', 'src'))
from ciceroscm import CICEROSCM

# %% [markdown]
# Define some input handling functions to give us example inputs

# %%
from ciceroscm.input_handler import read_inputfile,read_components,read_natural_emissions

# %% [markdown]
# Define a function to convert model output to a dataframe

# %% [markdown]
# set up input directories

# %%
test_data_dir = os.path.join(os.getcwd(), '../../', 'tests', 'test-data')

# %% [markdown]
# # Read in datafiles into dataframes

# %%
# NBVAL_IGNORE_OUTPUT
#Read gas parameters
gaspam =read_components(test_data_dir + '/gases_v1RCMIP.txt')
gaspam.head()

# %% [markdown]
# Read natural emissions

# %%
# NBVAL_IGNORE_OUTPUT
df_nat_ch4 =read_natural_emissions(test_data_dir + '/natemis_ch4.txt','CH4')
df_nat_n2o =read_natural_emissions(test_data_dir + '/natemis_n2o.txt','N2O')
df_nat_ch4.head()


# %% [markdown]
# Read forcing

# %%
df_ssp2_conc =read_inputfile(test_data_dir + '/ssp245_conc_RCMIP.txt')
df_ssp2_conc.head()

# %%
emi_input =read_inputfile(test_data_dir + '/ssp245_em_RCMIP.txt')
emi_input.rename(columns={"CO2": "CO2_FF", "CO2.1": "CO2_AFOLU"}, inplace=True)
emi_input.head()

# %% [markdown]
# # Set up model run with defined input variables

# %%
scendata={
            "gaspam_data": gaspam,
            "emstart": 1750,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2100,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
            "udir": test_data_dir,
            "scenname": "ssp245",
        }

# %%
    calibdata = pd.DataFrame(
        data={
            "Variable Name": [
                "Heat Content|Ocean",
                "Surface Air Ocean Blended Temperature Change",
            ],
            "Yearstart_norm": [1971, 1961],
            "Yearend_norm": [1971, 1990],
            "Yearstart_change": [2018, 2000],
            "Yearend_change": [2018, 2019],
            "Central Value": [320.69251537323, 0.5372],
            "sigma": [17.020342912051203, 0.039028311931729676],
        })

# %%
ordering=[
    "rlamdo",
    "akapa",
    "cpi",
    "W",
    "beto",
    "lambda",
    "mixed",
    "qo3",
    "qdirso2",
    "qindso2",
    "qbc",
    "qoc",
    "beta_f",
    "mixed_carbon",
    "qbmb",
    "qh2o_ch4",
]

# %%
len(ordering)

# %%
prior_flat_array = np.array(
    [
        [5, 25],
        [0.06, 0.8],
        [0.161, 0.569],
        [0.55, 2.55],
        [0, 7],
        [2 / 3.71, 5 / 3.71],
        [25, 125],
        [0.4, 0.6],
        [-0.55, -0.2],
        [-1.5, -0.5],
        [0.1, 0.2],
        [-0.1, -0.06],
        [0.110, 0.465],
        [25, 125],
        [0, 2],
        [0.08, 0.1],
    ]
)

# %%
prior_flat_array.shape

# %%
    testconfig = _ConfigDistro(
        distro_array=prior_flat_array,
        setvalues={
            "threstemp": 7.0,
            "lm": 40,
            "ldtime": 12,
        },
        ordering=ordering,
    )

# %%
len(testconfig.ordering)

# %%


# %%
# NBVAL_IGNORE_OUTPUT
scen = 'test'
cscm_dir=CICEROSCM({
            "gaspam_data": gaspam,
            "emstart": 1750,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2100,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
        })

# %%
distrorun1 = DistributionRun(testconfig, numvalues=50000)
output_vars = ["Heat Uptake", "Surface Air Temperature Change"]


# %%
def get_results(cfg):
    try:
        cscm_dir._run({"results_as_dict": True},pamset_udm=cfg['pamset_udm'],pamset_emiconc=cfg['pamset_emiconc'])
        res=cscm_dir.results

    except:
        res=None
    return [cfg,res]

# %%
def run_parallel(cfgs,nworkers=4):
    results=len(cfgs)*[None]
    with ProcessPoolExecutor(nworkers) as exe:
            # execute tasks concurrently and process results in order
            pres=list(tqdm(exe.map(get_results, cfgs)))
            for result in pres:
                # get the corresponding index of the config
                ind=int(result[0]['Index'])
                # put it in the right element of the results vector
                results[ind]=result[1]
    return results

# %%
results=run_parallel(distrorun1.cfgs,nworkers=200)

# %%
flds=['dT_glob','OHC700','concentrations-CO2','RIB_glob']

# %%
isgd=np.where([r!=None for r in results])[0].astype(int)

# %%
fresults = [results[i] for i in isgd]
fcfgs=[distrorun1.cfgs[i] for i in isgd]

# %%
def to_df(rs):

    out=pd.concat([pd.DataFrame(v).reset_index(drop=True) for k, v in rs.items()], keys = rs.keys(),axis=1) 
    out.index=rs['forcing'].index
    return out

# %%
def make_ensdf(results,flds):
    resdf=[]
    for i,res in enumerate(results):
        df=to_df(res)
        test_list=df.columns.map('{0[0]}-{0[1]}'.format).tolist()
        df.columns=[sub.replace('-0', '') for sub in test_list]
        resdf.append(df[flds].unstack())
    ensdf=pd.concat(resdf,axis=1)
    return ensdf

# %%
df=make_ensdf(fresults,flds)

# %%
issane=np.where(df.max()<1e8)[0]

# %%
df

# %%
df1 = df.iloc[:,issane]
fcfgs1=[fcfgs[i] for i in issane]

# %%
df1.index.names = ['variable', 'year']
df1.columns.names=['run_id']

# %%
def plot_range(df, var, ax,col='k'):
    Tdf=df.xs(var).T
    lower = Tdf.quantile(0.10)
    upper = Tdf.quantile(0.90)
    ax.fill_between(Tdf.columns, lower, upper, color=col, alpha=0.2,edgecolor=None)
    ax.plot(Tdf.columns, Tdf.mean(), color=col)



# %%
fig, ax = plt.subplots( 1,len(flds) ,figsize=(15, 4))

ax=ax.flatten()
for i,f in enumerate(flds):
    plot_range(df1, f, ax[i])
    ax[i].set_title(f)
ax[0].set_ylim([0,5])
ax[1].set_ylim([0,200])
ax[2].set_ylim([200,700])
ax[3].set_ylim([0,5])



# %%
dates=[1955,1975,1995,2008,2018]
targ=df1.loc[(slice(None), dates), :].T
nflds=targ.shape[1]


# %%
targ['OHC700']=(targ['OHC700'].T-targ['OHC700'][1995]).T

# %%
targ['OHC700']

# %%
def merge_dicts(dc):
    x=dc['pamset_udm']
    y=dc['pamset_emiconc']
    z = x.copy()
    z.update(y)
    return z

# %%
pdict=fcfgs1
mdict=[ merge_dicts(d) for d in pdict ]
pmat=pd.DataFrame(mdict)

# %%
parammat=pmat.loc[:, (pmat != pmat.iloc[0]).any()]
parammat

# %%
store = pd.HDFStore('data/data.h5')
store['targ'] = targ
store['parammat'] = parammat
store.close()


