# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # CICERO SCM notebook example - interactive input

# %% [markdown]
# Import some stuff

# %%
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import warnings
try:
    from pandas.core.common import SettingWithCopyWarning
except:
    from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

# %% [markdown]
# Import the model

# %%
sys.path.insert(0,os.path.join(os.getcwd(), '../', 'src'))
from ciceroscm import CICEROSCM

# %% [markdown]
# Define some input handling functions to give us example inputs

# %%
from ciceroscm.input_handler import read_inputfile,read_components,read_natural_emissions


# %% [markdown]
# Define a function to convert model output to a dataframe

# %%
def to_df(cscm):
    out=pd.concat([pd.DataFrame(v) for k, v in cscm.results.items()], axis = 1, keys = list(cscm.results.keys()))
    return out


# %% [markdown]
# set up input directories

# %%
test_data_dir = os.path.join(os.getcwd(), '../', 'tests', 'test-data')

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
# Here we will set up 4 different model runs to illustrate structural switches:
# 1. Default CICERO SCM setup (diffusion ocean, standard carbon cycle)
# 2. CICERO SCM with two-layer ocean model, but standard carbon cycle
# 3. CICERO SCM with box carbon cycle model, but standard diffusion ocean
# 4. CICERO SCM with both two-layer ocean model and box carbon cycle model

# %%
# NBVAL_IGNORE_OUTPUT
scen = 'test'
cscm_dir=CICEROSCM({
            "gaspam_data": gaspam,
            "emstart": 1751,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2100,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
        })

cscm_two_layer_ocean=CICEROSCM({
            "gaspam_data": gaspam,
            "emstart": 1751,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2100,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
            "thermal_model":"twolayer"
        })
cscm_box_carbon=CICEROSCM({
            "gaspam_data": gaspam,
            "emstart": 1751,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2100,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
            "carbon_cycle_model":"box"
        })
cscm_box_carbon_and_two_layer=CICEROSCM({
            "gaspam_data": gaspam,
            "emstart": 1751,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2100,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
            "carbon_cycle_model":"box",
            "thermal_model":"twolayer"
        })

# %% [markdown]
# and run it!

# %% [markdown]
# ## Parameter Organization
#
# The model parameters are organized into different categories:
#
# - **pamset_udm**: Physical oceanographic parameters (upwelling-diffusion model)
# - **pamset_emiconc**: Emissions and concentrations parameters, including forcing coefficients
# - **pamset_carbon**: Carbon cycle specific parameters
#
# In this example we use default parameters for all models, just to display the various structural switch options.

# %%
# NBVAL_IGNORE_OUTPUT
cscm_dir._run({
            "results_as_dict":True,
            "carbon_cycle_outputs":True
        })
cscm_two_layer_ocean._run({
            "results_as_dict":True,
            "carbon_cycle_outputs":True
        })   
cscm_box_carbon._run({
            "results_as_dict":True,
            "carbon_cycle_outputs":True
        })   
cscm_box_carbon_and_two_layer._run({
            "results_as_dict":True,
            "carbon_cycle_outputs":True
        })   

# %% [markdown]
# Convert the output to a dataframe for easy handling

# %%

df_temp = to_df(cscm_dir)
df_temp_two_layer = to_df(cscm_two_layer_ocean)
df_temp_box_carbon = to_df(cscm_box_carbon)
df_temp_box_carbon_and_two_layer = to_df(cscm_box_carbon_and_two_layer)


# %% [markdown]
# # Plot output

# %%
# NBVAL_IGNORE_OUTPUT
fig, axs = plt.subplots(nrows=2, ncols=2,figsize=(15,10))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation')

df_temp['dT_glob'].plot(ylabel='(K)',ax=axs[0])
df_temp_box_carbon['dT_glob'].plot(ylabel='(K)',ax=axs[0],linestyle=':')
df_temp_two_layer['dT_glob'].plot(ylabel='(K)',ax=axs[0],linestyle=':')
df_temp_box_carbon_and_two_layer['dT_glob'].plot(ylabel='(K)',ax=axs[0],linestyle=':')
axs[0].legend(['GMST default','GMST box carbon','GMST two layer ocean','GMST box carbon + two layer ocean'])

df_temp['concentrations']['CO2'].plot(ylabel='(ppm)',ax=axs[1])
df_temp_box_carbon['concentrations']['CO2'].plot(ylabel='(ppm)',ax=axs[1], linestyle=':')
df_temp_two_layer['concentrations']['CO2'].plot(ylabel='(ppm)',ax=axs[1], linestyle=':')
df_temp_box_carbon_and_two_layer['concentrations']['CO2'].plot(ylabel='(ppm)',ax=axs[1], linestyle=':')
axs[1].legend([r'CO$_2$ concentration','Box carbon','Two layer ocean','box carbon + two layer ocean'])

# Plot also Carbon cycle outputs and OHC
df_temp['OHCTOT'].plot(ylabel='(K)',ax=axs[2])
df_temp_box_carbon['OHCTOT'].plot(ylabel='(K)',ax=axs[2],linestyle=':')
df_temp_two_layer['OHCTOT'].plot(ylabel='(K)',ax=axs[2],linestyle=':')
df_temp_box_carbon_and_two_layer['OHCTOT'].plot(ylabel='(K)',ax=axs[2],linestyle=':')
axs[2].legend(['OHC default','box carbon','two layer ocean','box carbon + two layer ocean'])


df_temp["carbon cycle"]["Ocean carbon flux"].plot(ylabel='(W m-2)',ax=axs[3])
df_temp_box_carbon["carbon cycle"]["Ocean carbon flux"].plot(ylabel='(W m-2)',ax=axs[3], linestyle=':')
df_temp_two_layer["carbon cycle"]["Ocean carbon flux"].plot(ylabel='(W m-2)',ax=axs[3], linestyle=':')
df_temp_box_carbon_and_two_layer["carbon cycle"]["Ocean carbon flux"].plot(ylabel='(W m-2)',ax=axs[3], linestyle=':')
axs[3].legend(['Ocean carbon flux default','box carbon','two layer ocean','box carbon + two layer ocean'])
