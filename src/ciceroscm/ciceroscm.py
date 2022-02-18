"""
CICEROSCM
"""
import logging
import os

import numpy as np
import pandas as pd

from .concentrations_emissions_handler import ConcentrationsEmissionsHandler
from .upwelling_diffusion_model import UpwellingDiffusionModel

LOGGER = logging.getLogger(__name__)

default_data_dir = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "default_data"
)


def check_inputfiles(cfg):
    """
    Check whether input files are present or not
    """
    if not os.path.exists(cfg["gaspamfile"]):
        raise FileNotFoundError(
            f"Concentration input file {cfg['concentrations_file']} not found"
        )
    if not os.path.exists(cfg["concentrations_file"]):
        raise FileNotFoundError(
            f"Concentration input file {cfg['concentrations_file']} not found"
        )
    if not os.path.exists(cfg["emissions_file"]):
        raise FileNotFoundError(
            f"Emission input file {cfg['emissions_file']} not found"
        )
    if "nat_ch4_file" not in cfg:
        LOGGER.warning(
            "Did not find prescribed nat_ch4_file or none was given. Looking in standard path",
        )
        cfg["nat_ch4_file"] = os.path.join(
            os.getcwd(), "input_OTHER", "NATEMIS", "natemis_ch4.txt"
        )
    if not os.path.exists(cfg["nat_ch4_file"]):
        raise FileNotFoundError(
            f"Natural emission input file {cfg['nat_ch4_file']} not found"
        )
    if "nat_n2o_file" not in cfg:
        LOGGER.warning(
            "Did not find prescribed nat_n2o_file or none was given. Looking in standard path",
        )
        cfg["nat_n2o_file"] = os.path.join(
            os.getcwd(), "input_OTHER", "NATEMIS", "natemis_n2o.txt"
        )
    if not os.path.exists(cfg["nat_n2o_file"]):
        raise FileNotFoundError(
            f"Natural emission input file {cfg['nat_n2o_file']} not found"
        )
    return cfg


def read_forc(forc_file):
    """
    Read in forcing from forc
    """
    components = False
    with open(forc_file, "r", encoding="utf8") as fread:
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

    # pylint: disable=too-many-instance-attributes
    # Consider whether this can be cut back later

    def __init__(self):
        """
        Intialise CICEROSCM
        """
        self.nystart = 1750
        self.nyend = 2100
        self.emstart = 1850
        self.rf = None
        self.results = {}

    def initialise_output_arrays(self, cfg):
        """
        Initialise dict with arrays to hold data for run
        """
        # Add test to check that nystart and nyend are numbers
        if "nystart" in cfg:
            self.nystart = int(cfg["nystart"])
        if "nyend" in cfg:
            self.nyend = int(cfg["nyend"])
        if "emstart" in cfg:
            self.emstart = int(cfg["emstart"])
        output_variables = [
            "OHC700",
            "OHCTOT",
            "RIB_glob",
            "RIB_N",
            "RIB_S",
            "dT_glob",
            "dT_NH",
            "dT_SH",
            "dT_glob_air",
            "dT_NH_air",
            "dT_SH_air",
            "dT_glob_sea",
            "dT_NH_sea",
            "dT_SHsea",
            "dSL(m)",
            "dSL_thermal(m)",
            "dSL_ice(m)",
            "Total_forcing",
        ]
        for output in output_variables:
            self.results[output] = np.zeros(self.nyend - self.nystart + 1)

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
            df_data = pd.read_csv(
                volc_datafile,
                header=None,
                skiprows=skiprows,
                nrows=nrows,
                delim_whitespace=True,
            )
        else:
            df_data = pd.read_csv(
                volc_datafile, header=None, nrows=nrows, delim_whitespace=True
            )

        df_data.set_axis(labels=indices, inplace=True)
        return df_data

    def read_in_volc_and_sun(self, cfg):
        """
        Read in solar and volcanic forcing and return them
        """
        if "sunvolc" in cfg and cfg["sunvolc"] == 1:
            # Possibly change to allow for other files
            # And for SH to differ from NH
            rf_volc_n = self.read_data_on_year_row(
                os.path.join(default_data_dir, "meanVOLCmnd_ipcc_NH.txt")
            )
            rf_volc_s = rf_volc_n
            # Test, adjust for not including spin up. See Gregory et al.
            # Se regneark.
            rf_volc_n = rf_volc_n + 0.371457071
            rf_volc_s = rf_volc_s + 0.353195076
            rf_sun = self.read_data_on_year_row(
                os.path.join(default_data_dir, "solar_IPCC.txt")
            )
        # Add support for sending filename in cfg
        else:
            indices = np.arange(self.nystart, self.nyend + 1)
            rf_volc_n = pd.DataFrame(
                data=np.zeros((self.nyend - self.nystart + 1, 12)),
                index=indices,
                columns=range(12),
            )
            rf_volc_s = rf_volc_n
            rf_sun = pd.DataFrame(data={0: np.zeros(self.nyend - self.nystart + 1)})
        return {"volc_n": rf_volc_n, "volc_s": rf_volc_s, "sun": rf_sun}

    def forc_set(self, yr, rf_sun):
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
        forc = forc + rf_sun.iloc[row_index, 0]
        return forc

    def add_year_data_to_output(self, values, forc, index):
        """
        Add single year output to output arrays
        """
        simple_outputs = ["OHC700", "OHCTOT"]
        for output in simple_outputs:
            self.results[output][index] = values[output]
        outputs_dict = {
            "RIB_glob": "RIB",
            "RIB_N": "RIBN",
            "RIB_S": "RIBS",
            "dT_glob": "dtemp",
            "dT_NH": "dtempnh",
            "dT_SH": "dtempsh",
            "dT_glob_air": "dtemp_air",
            "dT_NH_air": "dtempnh_air",
            "dT_SH_air": "dtempsh_air",
            "dT_glob_sea": "dtemp_sea",
            "dT_NH_sea": "dtempnh_sea",
            "dT_SHsea": "dtempsh_sea",
        }
        for output, name in outputs_dict.items():
            self.results[output][index] = values[name]
        self.results["dSL(m)"][index] = values["deltsl"][0] + values["deltsl"][1]
        self.results["dSL_ice(m)"][index] = values["deltsl"][1]
        self.results["dSL_thermal(m)"][index] = values["deltsl"][0]
        self.results["Total_forcing"][index] = forc

    def _run(
        self, cfg, pamset_udm={}, pamset_emiconc={}
    ):  # pylint: disable=dangerous-default-value
        """
        Run CICEROSCM
        """
        # Add something to adjust start and end of simulation
        self.initialise_output_arrays(cfg)
        # Setting up UDM
        udm = UpwellingDiffusionModel(pamset_udm)

        # Reading in solar and volcanic forcing
        rf_volc_sun = self.read_in_volc_and_sun(cfg)

        # Add support for sending filename in cfg
        rf_luc = self.read_data_on_year_row(
            os.path.join(default_data_dir, "IPCC_LUCalbedo.txt")
        )

        rf_run = False
        if "forc_file" in cfg:
            rf_run = True
            if not os.path.exists(cfg["forc_file"]):
                raise FileNotFoundError(
                    f"Forcing input file {cfg['forc_file']} not found"
                )
            self.rf = read_forc(cfg["forc_file"])
        else:

            cfg = check_inputfiles(cfg)
            pamset_emiconc["emstart"] = self.emstart
            pamset_emiconc["nystart"] = self.nystart
            pamset_emiconc["nyend"] = self.nyend
            ce_handler = ConcentrationsEmissionsHandler(cfg, pamset_emiconc,)

        for yr in range(self.nystart, self.nyend + 1):
            if not rf_run:
                ce_handler.emi2conc(yr)
                forc, fn, fs = ce_handler.conc2forc(
                    yr, rf_luc.loc[yr, 0], rf_volc_sun["sun"].loc[yr - self.nystart, 0]
                )

            else:
                forc = self.forc_set(yr, rf_volc_sun["sun"])
                fs = forc
                fn = forc
            values = udm.energy_budget(
                fn,
                fs,
                rf_volc_sun["volc_n"].iloc[yr - self.nystart, :],
                rf_volc_sun["volc_s"].iloc[yr - self.nystart, :],
            )
            self.add_year_data_to_output(values, forc, yr - self.nystart)

        if not rf_run:
            ce_handler.write_output_to_files(cfg)

        self.write_data_to_file(cfg, rf_run)

    def write_data_to_file(self, pamset, rf_run=True):
        """
        Write results to files after run
        """
        if "output_prefix" in pamset:
            # Make os independent?
            outdir = os.path.join(os.getcwd(), pamset["output_prefix"])
        else:
            outdir = os.path.join(os.getcwd(), "output")

        indices = np.arange(self.nystart, self.nyend + 1)
        df_ohc = pd.DataFrame(
            data={
                "Year": indices,
                "OHC700": self.results["OHC700"],
                "OHCTOT": self.results["OHCTOT"],
            }
        )
        list_rib = ["RIB_glob", "RIB_N", "RIB_S"]
        df_rib = pd.DataFrame(data={"Year": indices})
        for vari in list_rib:
            df_rib[vari] = self.results[vari]
        list_temp = [
            "dT_glob",
            "dT_NH",
            "dT_SH",
            "dT_glob_air",
            "dT_NH_air",
            "dT_SH_air",
            "dT_glob_sea",
            "dT_NH_sea",
            "dT_SHsea",
            "dSL(m)",
            "dSL_thermal(m)",
            "dSL_ice(m)",
        ]
        df_temp = pd.DataFrame(data={"Year": indices})
        for vari in list_temp:
            df_temp[vari] = self.results[vari]
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
        if rf_run:
            df_forc = pd.DataFrame(
                data={"Year": indices, "Total_forcing": self.results["Total_forcing"]}
            )
            df_forc.to_csv(
                os.path.join(outdir, "output_forc.txt"),
                sep="\t",
                index=False,
                float_format="%.5e",
            )
