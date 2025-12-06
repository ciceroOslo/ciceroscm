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
#     display_name: ciceroscm (3.8.10)
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
from ciceroscm.pub_utils import make_regional_aerosol_gaspamdata


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
#Read gas parameters and make a regional aerosol gas pamdata version
df_gas =read_components(test_data_dir + '/gases_vupdate_2022_AR6.txt')
reg_aerosol_RF_file = os.path.join(test_data_dir, "HTAP_reg_aerosol_RF.txt")
reg_aerosol_df = pd.read_csv(reg_aerosol_RF_file, sep="\t", index_col=0)
reg_aerosol_df.rename(columns={"sulfate": "SO2"}, inplace=True)
df_gas_updated = make_regional_aerosol_gaspamdata(df_gas, reg_aerosol_df)

df_gas_updated.tail()

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
emi_input =read_inputfile(test_data_dir + '/ssp245_with_regional_aerosols_em_RCMIP.txt')
emi_input.rename(columns={"CO2": "CO2_FF", "CO2.1": "CO2_AFOLU"}, inplace=True)
emi_input.head()

# %% [markdown]
# # Set up model run with defined input variables

# %%
# NBVAL_IGNORE_OUTPUT
scen = 'test'
cscm_dir=CICEROSCM({
            "gaspam_data": df_gas_updated,
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
# ### Carbon Cycle Parameters
#
# Carbon cycle parameters should be passed in `pamset_carbon`. This includes:
#
# - Physical parameters: `beta_f`, `mixed_carbon`, `npp0`, `ml_w_sigmoid`, `ml_fracmax`, `ml_t_half`, etc.
# - Function parameters: Can use either dictionary or flat format
#
# ### Using Flat Carbon Cycle Parameters
#
# The carbon cycle function parameters can now be specified using "flat" parameter names instead of requiring dictionary structures. For example, instead of:
#
# ```python
# rb_function = {"coeffs": [0.5, 0.25, 0.25], "timescales": [2.5, 10.0, 60.0]}
# ```
#
# You can now use:
#
# ```python
# 'rb_coef0': 0.5, 'rb_coef1': 0.25, 'rb_coef2': 0.25,
# 'rb_tim0': 2.5, 'rb_tim1': 10.0, 'rb_tim2': 60.0
# ```

# %%
# NBVAL_IGNORE_OUTPUT
cscm_dir._run({
            "results_as_dict":True,
            "carbon_cycle_outputs":True
        },
    pamset_udm={"threstemp": 7.0, #scales vertical velocity as a function of mixed layer temperature
                    "rlamdo":16.0,#air-sea heat exchange coefficient (wm^-2K^-1)
                    "akapa":0.634, #vertical heat diffusivity
                    "cpi":0.4, #temperature change ratio: polar to nonpolar region
                    "W":4, #vertical velocity (m/yr)
                    "beto":3.5, #ocean interhemispheric heat exchange coeff (Wm^-2K^-1)
                    "lambda":0.54,
                    "mixed":60.0,  #mixed layer depth
                    "foan":0.61, #fraction of northern hemisphere covered by ocean
                    "foas":0.81, #fraction of northern hemisphere covered by ocean
                    "ebbeta":0.0,#atmospheric interhemispheric heat exchange 
                    "fnso":0.7531, #ocean area ratio, northern to southern hemisphere
                    "lm":40, #number of vertical layers
                    "ldtime":12,
                   },
    pamset_emiconc={"qbmb": 0.0,
                    "qo3": 0.5,
                    "qdirso2": -0.00308,
                    "qindso2": -0.97 / 57.052577209999995,
                    "qbc": 0.0279,
                    "qoc": -0.00433,
                    "qh2o_ch4": 0.091915,
                    "ref_yr": 2010
                    },
    pamset_carbon={
                    "beta_f": 0.,
                    "mixed_carbon": 75.0,
                    "qnmvoc": 0.0,
                    "qnh3": 0.0,
                    "qnox": 0.0,
                    "npp0": 60.0,
                    "npp_t_half": 0.5,
                    "npp_w_sigmoid": 4,
                    "npp_t_threshold": 6,
                    "npp_w_threshold": 4,
                    "ml_w_sigmoid": 3.0,
                    "ml_fracmax": 0.5,
                    "ml_t_half": 0.5,
                    'rb_coef0': 0.5,
                    'rb_coef1': 0.25,
                    'rb_coef2': 0.25,
                    'rb_tim0': 2.5,
                    'rb_tim1': 10.0,
                    'rb_tim2': 60.0,
                    'rs_coef0': 0.1,
                    'rs_coef1': 0.6,
                    'rs_coef2': 0.15,
                    'rs_coef3': 0.15,
                    'rs_tim0': .8,
                    'rs_tim1': 7,
                    'rs_tim2': 80,
                    'beta_f': 1.0
                }
            )   



# %% [markdown]
# Convert the output to a dataframe for easy handling

# %%

df_temp = to_df(cscm_dir)

# %% [markdown]
# # Plot output

# %%
# NBVAL_IGNORE_OUTPUT
fig, axs = plt.subplots(nrows=2, ncols=4,figsize=(30,10))
axs=axs.flatten()
fig.suptitle('CICERO SCM simulation')

df_temp['dT_glob'].plot(ylabel='(K)',ax=axs[0])
df_temp['dT_NH'].plot(ylabel='(K)',ax=axs[0],linestyle=':')
df_temp['dT_SH'].plot(ylabel='(K)',ax=axs[0],linestyle=':')
axs[0].legend(['Global mean temperature','NH temperature','SH temperature'])

df_temp['forcing']['Total_forcing'].plot(ylabel='(ppm)',ax=axs[1])
axs[1].legend([r'Total Forcing'])

aerosols = ["SO2", "BC", "OC"]
regions = ["ASIA", "LAM", "MAF", "REF", "OECD"]
for aer_num, aerosol in enumerate(aerosols):
    for region in regions:
        df_temp["emissions"][f"{aerosol}_{region}"].plot(
            ylabel=f'(Tg {aerosol})'
            ,ax=axs[2 + aer_num*2]
            )
        df_temp["forcing"][f"{aerosol}_{region}"].plot(
            ylabel='(W m-2)'
            ,ax=axs[3 + aer_num*2]
            )
    axs[2 + aer_num*2].legend([f'{aerosol}_{region}' for region in regions])
    axs[3 + aer_num*2].legend([f'{aerosol}_{region}' for region in regions])
plt.show()
