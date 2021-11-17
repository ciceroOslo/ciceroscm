"""
CICEROSCM
"""
import logging
import os
import sys
import numpy as np
import pandas as pd

from .upwelling_diffusion_model import UpwellingDiffusionModel
from .concentrations_emissions_handler import ConcentrationsEmissionsHandler


LOGGER = logging.getLogger(__name__)

default_data_dir = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "default_data"
)



def check_pamset(pamset):
    """
    Check that parameterset has necessary values for run
    Otherwise set to default values
    """
    required = {
        "rlamdo": 16.0,
        "akapa": 0.634,
        "cpi": 0.4,
        "W": 4.0,
        "beto": 3.5,
        "threstemp": 7.0,
        "lambda": 0.540,
        "mixed": 60.0,
    }
    for pam in required:
        if pam not in pamset:
            LOGGER.warning(
                "Parameter %s not in pamset. Using default value %f"%
                (pam, required[pam]),
            )
            pamset[pam] = required[pam]
        elif not isinstance(pamset[pam], int) and not isinstance(pamset[pam], float):
            LOGGER.warning(
                "Parameter %s must be a number. Using default value %f"%
                (pam, required[pam]),
            )
            pamset[pam] = required[pam]
    return pamset


def check_pamset_nonforc_run(pamset):
    """
    Check that parameterset has necessary values for run
    Otherwise set to default values
    """
    required = {
        "lamb": 0.8,
        "qbmb": 0.03,
        "qo3": 0.34,
        "qdirso2":-0.4,
        "qindso2": -0.8,
        "qbc": 0.22,
        "qoc":-0.05
    }
    for pam in required:
        if pam not in pamset:
            LOGGER.warning(
                "Parameter %s not in pamset. Using default value %f"%
                (pam, required[pam]),
            )
            pamset[pam] = required[pam]
        elif not isinstance(pamset[pam], int) and not isinstance(pamset[pam], float):
            LOGGER.warning(
                "Parameter %s must be a number. Using default value %f"%
                (pam, required[pam]),
            )
            pamset[pam] = required[pam]
    return pamset


def read_forc(forc_file):
    """
    Read in forcing from forc
    """
    components = False
    with open(forc_file, "r") as fread:
        first_line = fread.readline()
        if first_line[:4].lower() == "year":
            components = True
    if not components:
        df_forc = np.loadtxt(forc_file)
    else:
        # Decide on formatting for this
        df_forc = pd.read_csv(forc_file, index_col=0)
    return df_forc


class CICEROSCM:
    """
    Main ciceroscm class
    """

    def __init__(self):
        """
        Intialise CICEROSCM
        """
        self.nystart = 1750
        self.nyend = 2100
        """
        self.idtm = 24
        self.scenstart = 1991
        self.scenend = self.nyend
        self.emend = self.nyend
        self.nyears = self.nyend-self.nystart +1
        self.lamb = 0.8
        self.qbmb = 0.03
        self.qo3 = 0.34
        self.qdirso2 = -0.4
        self.qindso2 = -0.8
        self.qbc = 0.22
        self.qoc = -0.05
        """
        # Dont need the above?
        self.rf = None

    def initialise_output_arrays(self, pamset):
        """
        Initialise dataframe to hold data for run
        """
        # Add test to check that nystart and nyend are numbers
        if "nystart" in pamset:
            self.nystart = int(pamset["nystart"])
        if "nyend" in pamset:
            self.nyend = int(pamset["nyend"])
        self.ohc_700 = np.zeros(self.nyend - self.nystart + 1)
        self.ohc_tot = np.zeros(self.nyend - self.nystart + 1)
        self.rib = np.zeros(self.nyend - self.nystart + 1)
        self.rib_n = np.zeros(self.nyend - self.nystart + 1)
        self.rib_s = np.zeros(self.nyend - self.nystart + 1)
        self.dT_glob = np.zeros(self.nyend - self.nystart + 1)
        self.dT_NH = np.zeros(self.nyend - self.nystart + 1)
        self.dT_SH = np.zeros(self.nyend - self.nystart + 1)
        self.dT_glob_air = np.zeros(self.nyend - self.nystart + 1)
        self.dT_glob_NH_air = np.zeros(self.nyend - self.nystart + 1)
        self.dT_glob_SH_air = np.zeros(self.nyend - self.nystart + 1)
        self.dT_glob_sea = np.zeros(self.nyend - self.nystart + 1)
        self.dT_glob_NH_sea = np.zeros(self.nyend - self.nystart + 1)
        self.dT_glob_SH_sea = np.zeros(self.nyend - self.nystart + 1)
        self.dSL = np.zeros(self.nyend - self.nystart + 1)
        self.dSL_ice = np.zeros(self.nyend - self.nystart + 1)
        self.dSL_thermal = np.zeros(self.nyend - self.nystart + 1)
        self.forcing = np.zeros(self.nyend - self.nystart + 1)

    def read_data_on_year_row(self, volc_datafile):
        """
        Read in data from file with no headers
        and each year being a row. Typically the format for
        volcano and solar data
        """
        indices = np.arange(self.nystart, self.nyend + 1)
        nrows = len(indices)
        if self.nystart > 1750:
            skiprows = self.nystart - 1750
            df = pd.read_csv(
                volc_datafile,
                header=None,
                skiprows=skiprows,
                nrows=nrows,
                delim_whitespace=True,
            )
        else:
            df = pd.read_csv(
                volc_datafile, header=None, nrows=nrows, delim_whitespace=True
            )

        df.set_axis(labels=indices, inplace=True)
        return df

    def forc_set(self, yr):
        """
        Read the forcing for this year
        """
        row_index = yr - self.nystart
        # Add support for other forcing formats
        if isinstance(self.rf, np.ndarray):
            # Add luc albedo later
            forc = self.rf[row_index]  # + self.rf_luc.iloc[row_index][0]
        else:
            forc = self.rf["total"][row_index]
        forc = forc + self.rf_sun.iloc[row_index, 0]
        return forc

        
    def add_year_data_to_output(self, values, forc, index):
        """
        Add single year output to output arrays
        """
        self.ohc_700[index] = values["OHC700"]
        self.ohc_tot[index] = values["OHCTOT"]
        self.rib[index] = values["RIB"]
        self.rib_n[index] = values["RIBN"]
        self.rib_s[index] = values["RIBS"]
        self.dT_glob[index] = values["dtemp"]
        self.dT_NH[index] = values["dtempnh"]
        self.dT_SH[index] = values["dtempsh"]
        self.dT_glob_air[index] = values["dtemp_air"]
        self.dT_glob_NH_air[index] = values["dtempnh_air"]
        self.dT_glob_SH_air[index] = values["dtempsh_air"]
        self.dT_glob_sea[index] = values["dtemp_sea"]
        self.dT_glob_NH_sea[index] = values["dtempnh_sea"]
        self.dT_glob_SH_sea[index] = values["dtempsh_sea"]
        self.dSL[index] = values["deltsl"][0] + values["deltsl"][1]
        self.dSL_ice[index] = values["deltsl"][1]
        self.dSL_thermal[index] = values["deltsl"][0]
        self.forcing[index] = forc

    def _run(self, pamset, cfg):
        """
        Run CICEROSCM
        """
        # Parameters as dict, possible make pamfilereading an option
        pamset = check_pamset(pamset)
        # Add something to adjust start and end of simulation
        
        self.initialise_output_arrays(pamset)
        # Setting up UDM
        udm = UpwellingDiffusionModel(pamset)
                
        # Reading in solar and volcanic forcing
        if "sunvolc" in pamset and pamset["sunvolc"] == 1:
            # Possibly change to allow for other files
            # And for SH to differ from NH
            rf_volc_n = self.read_data_on_year_row(
                os.path.join(default_data_dir, "meanVOLCmnd_ipcc_NH.txt")
            )
            rf_volc_s = rf_volc_n
            # Test, adjust for not including spin up. See Gregory et al.
            # Se regneark.
            rf_volc_n =rf_volc_n + 0.371457071
            rf_volc_s = rf_volc_s + 0.353195076
            self.rf_sun = self.read_data_on_year_row(
                os.path.join(default_data_dir, "solar_IPCC.txt")
            )
        # Add support for sending filename in pamset
        else:
            indices = np.arange(self.nystart, self.nyend + 1)
            rf_volc_n = pd.DataFrame(
                data=np.zeros((self.nyend - self.nystart + 1, 12)),
                index=indices,
                columns=range(12),
            )
            rf_volc_s = rf_volc_n
            self.rf_sun = pd.DataFrame(
                data={0: np.zeros(self.nyend - self.nystart + 1)}
            )

        # Add support for sending filename in pamset
        self.rf_luc = self.read_data_on_year_row(
            os.path.join(default_data_dir, "IPCC_LUCalbedo.txt")
        )

        rf_run = False
        conc_run = False
        if "forc_file" in cfg:
            rf_run = True
            if not os.path.exists(cfg["forc_file"]):
                raise FileNotFoundError(
                    "Forcing input file {} not found".format(cfg["forc_file"])
                )
            self.rf = read_forc(cfg["forc_file"])
        elif "concentrations_file" in cfg:
            conc_run = True
            if not os.path.exists(cfg["concentrations_file"]):
                raise FileNotFoundError(
                    "Concentration input file {} not found".format(cfg["concentrations_file"])
                )
            pamset = check_pamset_nonforc_run(pamset)
            ce_handler = ConcentrationsEmissionsHandler(pamset["gaspamfile"], cfg["concentrations_file"], cfg["emissions_file"], pamset)
            
        for yr in range(self.nystart, self.nyend + 1):
            if not rf_run:
                ce_handler.emi2conc(yr)
                forc, fn, fs = ce_handler.conc2forc(yr,self.rf_luc.loc[yr,0],self.rf_sun.loc[yr-self.nystart, 0]) 
            else:
                forc = self.forc_set(yr)
                fs = forc
                fn = forc
            values = udm.energy_budget(
                fn,
                fs,
                rf_volc_n.iloc[yr - self.nystart, :],
                rf_volc_s.iloc[yr - self.nystart, :],
            )
            self.add_year_data_to_output(values, forc, yr - self.nystart)
            
        self.write_data_to_file(pamset)

    def write_data_to_file(self, pamset):
        """
        Write results to files after run
        """
        if "output_prefix" in pamset:
            # Make os independent?
            outdir = os.path.join(os.getcwd(), pamset["output_prefix"])
        else:
            outdir = os.path.join(os.getcwd(), "output")

        indices = np.arange(self.nystart, self.nyend + 1)
        df_forc = pd.DataFrame(data={"Year": indices, "Total_forcing": self.forcing,})
        df_ohc = pd.DataFrame(
            data={"Year": indices, "OHC700": self.ohc_700, "OHCTOT": self.ohc_tot}
        )
        df_rib = pd.DataFrame(
            data={
                "Year": indices,
                "RIB_glob": self.rib,
                "RIB_N": self.rib_n,
                "RIB_S": self.rib_s,
            }
        )
        df_temp = pd.DataFrame(
            data={
                "Year": indices,
                "dT_glob": self.dT_glob,
                "dT_NH": self.dT_NH,
                "dT_SH": self.dT_SH,
                "dT_glob_air": self.dT_glob_air,
                "dT_NH_air": self.dT_glob_NH_air,
                "dT_SH_air": self.dT_glob_SH_air,
                "dT_glob_sea": self.dT_glob_sea,
                "dT_NH_sea": self.dT_glob_NH_sea,
                "dT_SHsea": self.dT_glob_SH_sea,
                "dSL(m)": self.dSL,
                "dSL_thermal(m)": self.dSL_thermal,
                "dSL_ice(m)": self.dSL_ice,
            }
        )
        df_forc.to_csv(
            os.path.join(outdir, "output_forc.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        df_ohc.to_csv(
            os.path.join(outdir, "output_ohc.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        df_rib.to_csv(
            os.path.join(outdir, "output_rib.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        df_temp.to_csv(
            os.path.join(outdir, "output_temp.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
