import os
import sys

# Adding location of source code to system path
# os.path.dirname(__file__) gives the directory of
# current file. Put in updated path if running script from elsewhere
# os.path joins all the folders of a path together in a
# system independent way (i.e. will work equally well on Windows, linux etc)
sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))

from ciceroscm import CICEROSCM

data_dir = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")
pert_dir = os.path.join(os.path.dirname(__file__), "pre_script")
outdir = os.path.join(os.getcwd(), "output_test")

# Emission perturbation test:
cscm = CICEROSCM(
    {
        "gaspam_file": os.path.join(data_dir, "gases_v1RCMIP.txt"),
        "nystart": 1900,
        "emstart": 1950,
        "nyend": 2050,
        "concentrations_file": os.path.join(data_dir, "ssp245_conc_RCMIP.txt"),
        "emissions_file": os.path.join(data_dir, "ssp245_em_RCMIP.txt"),
        "nat_ch4_file": os.path.join(data_dir, "natemis_ch4.txt"),
        "nat_n2o_file": os.path.join(data_dir, "natemis_n2o.txt"),
        "perturb_em_file": os.path.join(pert_dir, "pertem_test.txt"),        
    },
)

cscm._run({"output_folder": outdir, "output_prefix": "emissions_perturbation"})

# Forcing perturbation test:
cscm = CICEROSCM(
    {
        "gaspam_file": os.path.join(data_dir, "gases_v1RCMIP.txt"),
        "nystart": 1900,
        "emstart": 1950,
        "nyend": 2050,
        "concentrations_file": os.path.join(data_dir, "ssp245_conc_RCMIP.txt"),
        "emissions_file": os.path.join(data_dir, "ssp245_em_RCMIP.txt"),
        "nat_ch4_file": os.path.join(data_dir, "natemis_ch4.txt"),
        "nat_n2o_file": os.path.join(data_dir, "natemis_n2o.txt"),
        "perturb_forc_file": os.path.join(pert_dir, "pertforc_test.txt"),        
    },
)

cscm._run({"output_folder": outdir, "output_prefix": "forcing_perturbation"})
