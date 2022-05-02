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

# Setup for forcing run from forcing file, with solar and volcanic forcing
cscm = CICEROSCM(
    {
        "sunvolc": 1,
        "nyend": 2100,
        "forc_file": os.path.join(data_dir, "CO2_1pros.txt"),
    },
)

cscm._run({"output_folder": outdir})
