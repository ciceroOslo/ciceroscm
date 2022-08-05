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
#os.getcwd() gets the path of where you are running from
outdir = os.path.join(os.getcwd(), "output_test")
cscm = CICEROSCM(
    {
        "gaspam_file": os.path.join(data_dir, "gases_v1RCMIP.txt"),
        "nyend": 2100,
        "concentrations_file": os.path.join(data_dir, "ssp245_conc_RCMIP.txt"),
        "emissions_file": os.path.join(data_dir, "ssp245_em_RCMIP.txt"),
        "nat_ch4_file": os.path.join(data_dir, "natemis_ch4.txt"),
        "nat_n2o_file": os.path.join(data_dir, "natemis_n2o.txt"),
        "idtm": 24
    },
)

cscm._run(
    {"output_folder": outdir, "output_prefix": "pams_other"},
    pamset_udm={
        "rlamdo": 7.935264,
        "akapa": 0.33679135961526513,
        "cpi": 0.5513349,
        "W": 1.081163,
        "beto": 0.8246083,
        "threstemp": 6.0,
        "lambda": 0.6086934,
        "mixed": 78.46448,
        "foan": 0.62,
        "foas": 0.82,
        "ebbeta": 0.1,
        "fnso": 0.7532,
        "lm": 39,
        "ldtime": 10,
    },
    pamset_emiconc={
        "qbmb": 0.0,
        "qo3": 0.3,
        "qdirso2": -0.5203740910977426,
        "qindso2": -0.5852790880577513,
        "qbc": 0.2287316893023989,
        "qoc": -0.11778341014690517,
        "ref_yr": 2015,
    },
)
