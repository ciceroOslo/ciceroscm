import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from ciceroscm import concentrations_emissions_handler, input_handler

# Change data_dir to get inputs from somewhere else
data_dir = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test-data")

# Change these years if you wish to adjust the estimated time range:
nyend = 2500
nystart = 1750

# Choose the lifetime mode assumption for methane 
# lifetime_modes = ["TAR", "CONSTANT_12", "CONSTANT_from_file", "WIGLEY"]
# Change this line to your preferred lifetime:
lf_mode = "TAR"

pamset = {"nyend": nyend, "nystart": nystart, "lifetime_mode": lf_mode}
ih_temp = input_handler.InputHandler(pamset)

# If you want to fit from different data input files, then change these paths
em_data = ih_temp.read_emissions(os.path.join(data_dir, "ssp245_em_RCMIP.txt"))
conc_data =  input_handler.read_inputfile(os.path.join(data_dir, "ssp245_conc_RCMIP.txt"))
gaspam_data = input_handler.read_components(os.path.join(data_dir, "gases_vupdate_2022_AR6.txt"))

def get_lifetime(tracer, yr, ce_handler):
    """
    Get lifetime for tracer as calculated in the code
    """
    q = 1/ce_handler.df_gas["TAU1"][tracer]
    if tracer != "CH4":
        return q
    q = ce_handler.methane_lifetime(q, ce_handler.conc_in[tracer][yr-1], yr)
    return q

def old_get_nat_em_timeseries(tau, beta, em_nat, em_series, conc_series, sp="N2O"):
    """
    Get natural emissions using the previously used old natural emissions
    method. 
    
    We recommend using the ode method instead for more efficient
    and accurate calculations
    """
    em_nat_hist = np.full(len(conc_series), em_nat,np.double)
    CONC_NEW_last = conc_series[0]
    for i in range(len(conc_series)):
        not_done = True
        times_around = 0
        print(i)
        q = get_q(sp, i+nystart)
        
        while not_done: #and times_around < 100:
            pc = (em_series[i] + em_nat)/(beta)
            CONC_NEW = (pc/q+(CONC_NEW_last-pc/q)*np.exp(-q*1.0))
            print(CONC_NEW)
            print(conc_series[i])
            if (CONC_NEW/conc_series[i] >  1.005): 
                em_nat = em_nat*0.995
            elif (CONC_NEW/conc_series[i] < 0.995): 
                em_nat = em_nat*1.005
            else:
                em_nat_hist[i] = em_nat
                not_done = False
            times_around = times_around +1
        CONC_NEW_last = CONC_NEW
    return em_nat_hist

def ode_get_nat_em_timeseries(tau, beta, em_series, conc_series, ce_handler, sp = "N2O"):
    """
    Use ode exact solution backwards to get natural emissions
    needed to match the concentrations from anthropogenic emissions
    series
    """
    conc_extended = np.concatenate(([conc_series[0]], conc_series), axis = 0)
    if sp != "CH4":
        em_nat_hist = beta/tau/(1-np.exp(-1/tau))*(conc_series -conc_extended[:-1]*np.exp(-1/tau)) - em_series[:len(conc_series)]
    else:
        em_nat_hist = np.zeros(len(conc_series))
        for i in range(len(conc_series)):
            q = get_lifetime(sp, i+nystart, ce_handler)
            em_nat_hist[i] = beta*q/(1-np.exp(-q))*(conc_series[i] - conc_extended[i]*np.exp(-q)) - em_series[i]
    return em_nat_hist



nat_ch4_data = pd.DataFrame(
                data={"CH4": np.ones(nyend - nystart + 1) * 242.09},
                index=np.arange(nystart, nyend + 1),
            )
nat_n2o_data = pd.DataFrame(
                data={"CH4": np.ones(nyend - nystart + 1) * 242.09},
                index=np.arange(nystart, nyend + 1),
            )
ce_handler = concentrations_emissions_handler.ConcentrationsEmissionsHandler(input_handler.InputHandler({"gaspam_data":gaspam_data, "emissions_data":em_data, "concentrations_data": conc_data, "nat_ch4_data": nat_ch4_data, "nat_n2o_data": nat_n2o_data, "nystart": nystart, "nyend": nyend}), pamset)

ce_handler.reset_with_new_pams(ce_handler.pamset)
conc_y0 = ce_handler.conc_in.index.values[0]

# In principle you could use this for other tracers as well or do just a single tracer at a time
tracers = ["N2O", "CH4"]
year_out = range(nystart, nyend + 1)
for j, tracer in enumerate(tracers):
    tau = ce_handler.df_gas["TAU1"][tracer]
    beta = ce_handler.df_gas["BETA"][tracer]
    em_nat_const = ce_handler.df_gas["NAT_EM"][tracer]
    em_series = ce_handler.emis[tracer].values
    conc_series = ce_handler.conc_in[tracer][nystart-conc_y0:nyend+1-conc_y0].values
    print(print(conc_series))
    print(ce_handler.conc_in[tracer])

    # Option for old version calculation (comment in next three lines for old version):
    # em_nat_hist_old = old_get_nat_em_timeseries(tau, beta, em_nat_const, em_series, conc_series, sp=tracer)
    # em_nat_out_old = np.concatenate((em_nat_hist_old,  np.full((len(year_out)-len(em_nat_hist_old)), np.mean(em_nat_hist_old[-11:]),np.double)), axis=0)
    # np.savetxt(f"natemis_{tracer}_old_method_from_rcmip.txt", em_nat_out_old,fmt='%1.4f')

    # Option for new (recommended) version calculation (comment out next three lines for old version):
    em_nat_hist_ode = ode_get_nat_em_timeseries(tau, beta, em_series, conc_series, ce_handler, sp=tracer)
    em_nat_out_ode = np.concatenate((em_nat_hist_ode,  np.full((len(year_out)-len(em_nat_hist_ode)), np.mean(em_nat_hist_ode[-11:]),np.double)), axis=0)
    np.savetxt(f"natemis_{tracer}_ode_method_from_rcmip.txt", em_nat_out_ode,fmt='%1.4f')


