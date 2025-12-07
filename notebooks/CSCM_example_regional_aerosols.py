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
# Change sign of all values in the dataframe
reg_aerosol_df *= -1
#reg_aerosol_df.values = - reg_aerosol_df.values  # Convert from W/m2 to W/m2 per Tg
#sys.exit(4)
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
emi_input_vanilla =read_inputfile(test_data_dir + '/ssp245_em_RCMIP.txt')
emi_input_vanilla.rename(columns={"CO2": "CO2_FF", "CO2.1": "CO2_AFOLU"}, inplace=True)
emi_input_vanilla.head()

# %% [markdown]
# # Set up model run with defined input variables

# %%
# NBVAL_IGNORE_OUTPUT
scen = 'test'
cscm_dir=CICEROSCM({
            "gaspam_data": df_gas_updated,
            "emstart": 1850,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2100,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input,
            "nat_ch4_data": df_nat_ch4,
            "nat_n2o_data": df_nat_n2o,
            "idtm":24,
        })

scen = 'test'
cscm_dir_vanilla=CICEROSCM({
            "gaspam_data": df_gas,
            "emstart": 1850,  
            "conc_run":False,
            "nystart": 1750,
            "nyend": 2100,
            "concentrations_data": df_ssp2_conc,
            "emissions_data": emi_input_vanilla,
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
cscm_dir._run({"results_as_dict":True})   

cscm_dir_vanilla._run({"results_as_dict":True,})   


# %% [markdown]
# Convert the output to a dataframe for easy handling

# %%

df_temp = to_df(cscm_dir)
df_temp_vanilla = to_df(cscm_dir_vanilla)

# %% [markdown]
# # Plot output

# %%
# NBVAL_IGNORE_OUTPUT
fig, axs = plt.subplots(nrows=2, ncols=4,figsize=(30,10))
axs=axs.flatten()
fig.suptitle('CICERO SCM regional aerosols example run - SSP2-4.5 scenario', fontsize=16)

df_temp['dT_glob'].plot(ylabel='(K)',ax=axs[0])
df_temp['dT_NH'].plot(ylabel='(K)',ax=axs[0],linestyle=':')
df_temp['dT_SH'].plot(ylabel='(K)',ax=axs[0],linestyle=':')
df_temp_vanilla['dT_glob'].plot(ylabel='(K)',ax=axs[0],linestyle='--')
axs[0].legend(['GMST','NH temperature','SH temperature', 'GMST_non_regional'])
axs[0].set_title('Temperature')

df_temp['forcing']['Total_forcing'].plot(ylabel='(ppm)',ax=axs[1])
df_temp_vanilla['forcing']['Total_forcing'].plot(ylabel='(ppm)',ax=axs[1], linestyle='--')
axs[1].legend([r'Total Forcing', r'Total Forcing_non_regional'])
axs[1].set_title('Total Forcing')

aerosols = ["SO2", "BC", "OC"]
regions = ["ASIA", "LAM", "MAF", "REF", "OECD"]
for aer_num, aerosol in enumerate(aerosols):
    em_sum = pd.Series(0, index=df_temp.index)
    forc_sum = pd.Series(0, index=df_temp.index)
    for region in regions:
        df_temp["emissions"][f"{aerosol}_{region}"].plot(
            ylabel=f'(Tg {aerosol})',
            ax=axs[2 + aer_num*2]
            )
        df_temp["forcing"][f"{aerosol}_{region}"].plot(
            ylabel='(W m-2)',
            ax=axs[3 + aer_num*2]
            )
        em_sum += df_temp["emissions"][f"{aerosol}_{region}"]
        forc_sum += df_temp["forcing"][f"{aerosol}_{region}"]
    # Also plot the non-regional aerosol emissions and forcing
    em_sum.plot(
        ylabel=f'(Tg {aerosol})',
        ax=axs[2 + aer_num*2],
        linestyle='--'
        )
    forc_sum.plot(
        ylabel=f'ERF {aerosol} (W m-2)',
        ax=axs[3 + aer_num*2],
        linestyle='--'
        )
    df_temp_vanilla["emissions"][aerosol].plot(
        ylabel=f'(Tg {aerosol})',
        ax=axs[2 + aer_num*2],
        linestyle='--'
        )
    print(df_temp_vanilla["forcing"].keys())
    if aerosol == "SO2":
        df_temp_vanilla["forcing"]["SO4_DIR"].plot(
            ylabel=f'ERF {aerosol} (W m-2)',
            ax=axs[3 + aer_num*2],
            linestyle='--'
            )
    else:
        df_temp_vanilla["forcing"][aerosol].plot(
            ylabel=f'ERF {aerosol}(W m-2)',
            ax=axs[3 + aer_num*2],
            linestyle='--'
            )  
    legend_em =  [f'{aerosol}_{region}' for region in regions]
    legend_em.append(f"{aerosol}_total")
    legend_em.append(f'{aerosol}_non_regional')            
    axs[2 + aer_num*2].legend(legend_em)
    axs[2 + aer_num*2].set_title(f'Emissions of {aerosol}')
    axs[3 + aer_num*2].set_title(f'Forcing of {aerosol}')
    axs[3 + aer_num*2].legend(legend_em)
plt.show()
