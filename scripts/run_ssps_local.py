import os
import sys

# Adding location of source code to system path
# os.path.dirname(__file__) gives the directory of
# current file. Put in updated path if running script from elsewhere
# os.path joins all the folders of a path together in a
# system independent way (i.e. will work equally well on Windows, linux etc)
sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))

from ciceroscm import CICEROSCM

# Not system independent, can only be run on amoc or qbo
# If running on different system, this must be changed
# to where you have ssp input files
data_dir = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")
input_dir = "/div/amoc/CSCM/SCM_Linux_v2019/RCMIP/input/"

outdir = os.path.join(os.getcwd(), "./output_test")

ssps = ["119", "126", "245", "245methane", "370-lowNTCF", "370", "434", "460", "534-over", "585"]

for s in ssps:
    concf = f"{input_dir}ssp{s}_conc_RCMIP.txt"
    emif = f"{input_dir}ssp{s}_em_RCMIP.txt"

    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(data_dir, "gases_v1RCMIP.txt"),
            "concentrations_file": concf,
            "emissions_file": emif,
            "nat_ch4_file": os.path.join(data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(data_dir, "natemis_n2o.txt"),
        },
    )
    cscm._run({"output_folder": outdir, "output_prefix": f"ssp{s}"})
