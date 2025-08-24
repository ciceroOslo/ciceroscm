#!/usr/bin/env python
# coding: utf-8

# # CICERO SCM notebook - parallel application

# Import some stuff

# In[1]:


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

try:
    from pandas.core.common import SettingWithCopyWarning
except:
    from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.filterwarnings("ignore", message=".*Parameter.*")


# In[2]:


from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm


# In[3]:


import concurrent


# In[4]:


from ciceroscm.parallel._configdistro import _ConfigDistro
from ciceroscm.parallel.calibrator import Calibrator
from ciceroscm.parallel.distributionrun import DistributionRun


# Import the model

# In[5]:


sys.path.insert(0,os.path.join(os.getcwd(), '../', 'src'))
from ciceroscm import CICEROSCM


# Define some input handling functions to give us example inputs

# In[6]:


from ciceroscm import input_handler 


# Define a function to convert model output to a dataframe

# set up input directories

# In[7]:


test_data_dir = os.path.join(os.getcwd(), '../../', 'tests', 'test-data')


# # Read in datafiles into dataframes

# In[8]:


# NBVAL_IGNORE_OUTPUT
#Read gas parameters
gaspam = input_handler.read_components(test_data_dir + '/gases_vupdate_2022_AR6.txt')
print(gaspam.head())

#sys.exit(4)
# Read natural emissions

# In[9]:


# NBVAL_IGNORE_OUTPUT
df_nat_ch4 = input_handler.read_natural_emissions(test_data_dir + '/natemis_ch4.txt','CH4')
df_nat_n2o = input_handler.read_natural_emissions(test_data_dir + '/natemis_n2o.txt','N2O')
print(df_nat_ch4.head())


# Read forcing

# In[10]:


df_ssp2_conc =input_handler.read_inputfile(test_data_dir + '/ssp245_conc_RCMIP.txt', True)
print(df_ssp2_conc.head())
print(df_ssp2_conc.index)


# In[11]:

ih = input_handler.InputHandler({"nyend": 2025, "nystart": 1750, "emstart": 1850})
emi_input =ih.read_emissions(test_data_dir + '/ssp245_em_RCMIP.txt')
emi_input.rename(columns={"CO2": "CO2_FF", "CO2.1": "CO2_AFOLU"}, inplace=True)
print(emi_input.head())


# # Set up model run with defined input variables

# In[12]:


scendata={
            "gaspam_data": gaspam,
            "emstart": 1850,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2025,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
            "udir": test_data_dir,
            "scenname": "ssp245",
        }

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
    "qh2o_ch4",
]


# In[15]:


len(ordering)


# In[16]:


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
        [0.08, 0.1],
    ]
)


# In[17]:


print(prior_flat_array.shape)


# In[18]:


testconfig = _ConfigDistro(
    distro_array=prior_flat_array,
    setvalues={
        "threstemp": 7.0,
        "lm": 40,
        "ldtime": 12,
        "qbmb": 0
    },
    ordering=ordering,
)


print(len(testconfig.ordering))

scen = 'test'
cscm_dir=CICEROSCM({
            "gaspam_data": gaspam,
            "emstart": 1850,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2025,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24
        })

scenariodata = [{
            "gaspam_data": gaspam,
            "emstart": 1850,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2025,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
            "scenname" : "ssp245-short"
        }]


# In[21]:
calibdata = pd.DataFrame(
    data={
        "Variable Name": [
            "Heat Content|Ocean",
            "Surface Air Ocean Blended Temperature Change",
            "Effective Radiative Forcing|Aerosols",
            "Atmospheric Concentrations|CO2",
            "Ocean carbon flux", 
            "Biosphere carbon flux" 
        ],
        "Yearstart_norm": [1971, 1850, 1750, 1750, 1750, 1750],
        "Yearend_norm": [1971, 1900, 1750, 1750, 1750, 1750],
        "Yearstart_change": [2023, 2015, 2024, 2024, 2023, 2023],
        "Yearend_change": [2023, 2024, 2024, 2024, 2023, 2023],
        "Central Value": [484.82157000000063, 1.24, -1.09, 422.5, 2.880720693, 2.302042382],
        "sigma": [36.8551891022091 , 0.073, 0.6, 3.0, 0.4, 0.5],

    })

distrorun1 = DistributionRun(testconfig, numvalues=50000)
output_vars = calibdata["Variable Name"]
results = distrorun1.run_over_distribution(scenariodata, output_vars, max_workers=200)

# In[22]:

print(results.head())
print(results.shape)
#sys.exit(4)

# In[33]:
results_for_fit_dict = {}
for idx, data in calibdata.iterrows():
    print(data)
    results_sub = results.loc[results["variable"] == data["Variable Name"]]
    if data["Yearstart_norm"] == data["Yearend_norm"] and data["Yearstart_norm"]== 1750:
        results_for_fit_dict[data["Variable Name"]] = results_sub.iloc[:, data["Yearstart_change"]-1750+7].values
    elif data["Yearstart_norm"] == data["Yearend_norm"]:
        results_for_fit_dict[data["Variable Name"]] = (results_sub.iloc[:, data["Yearstart_change"]-1750+7] - results_sub.iloc[:,data["Yearstart_norm"]-1750+7]).values
    else:
        results_for_fit_dict[data["Variable Name"]] = ((
            results_sub.iloc[:,data["Yearstart_change"]-1750+7: data["Yearend_change"]-1750+8]).mean(axis = 1) - (
            results_sub.iloc[:,data["Yearstart_norm"]-1750+7: data["Yearend_norm"]-1750+8]).mean(axis = 1)
        ).values
print(results_for_fit_dict)
targ = pd.DataFrame(
    data = results_for_fit_dict
)
targ.index.set_names("run_id", inplace=True)
print(targ)

# TODO: Possibly put back sanity check?
pdict=distrorun1.cfgs


# In[36]:


def merge_dicts(dc):
    x=dc['pamset_udm']
    y=dc['pamset_emiconc']
    z = x.copy()
    z.update(y)
    return z


# In[41]:

mdict=[ merge_dicts(d) for d in pdict ]
pmat=pd.DataFrame(mdict)


# In[42]:


parammat=pmat.loc[:, (pmat != pmat.iloc[0]).any()]
parammat


# In[43]:


store = pd.HDFStore('data/data.h5')
store['targ'] = targ
store['parammat'] = parammat
store.close()

