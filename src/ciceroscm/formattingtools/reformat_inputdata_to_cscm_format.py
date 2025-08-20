"""
Common handler and converter of scenariodat for cicero
Class will be inherited by subversions to get specific
handling
"""

import csv

import numpy as np
import pandas as pd

cicero_comp_dict = {
    "CO2_lu": ["CO2|MAGICC AFOLU", 1],
    "CFC-113": ["CFC113", 1],
    "CFC-114": ["CFC114", 1],
    "SO2": ["Sulfur", 1],
    "NMVOC": ["VOC", 1],
    "CFC-11": ["CFC11", 1],
    "CFC-115": ["CFC115", 1],
    "CFC-12": ["CFC12", 1],
    "HCFC-141b": ["HCFC141b", 1],
    "HCFC-142b": ["HCFC142b", 1],
    "HCFC-22": ["HCFC22", 1],
    "H-1211": ["Halon1211", 1],
    "H-1301": ["Halon1301", 1],
    "H-2402": ["Halon2402", 1],
    "CO2": ["CO2|MAGICC Fossil and Industrial", 1],
    "CH4": ["CH4", 1],
    "N2O": ["N2O", 1],
    "CH3Br": ["CH3Br", 1],
    "CCl4": ["CCl4", 1],
    "CH3CCl3": ["CH3CCl3", 1],
    "HCFC-123": ["HCFC-123", 1],
    "HFC125": ["HFC125", 1],
    "HFC134a": ["HFC134a", 1],
    "HFC143a": ["HFC143a", 1],
    "HFC227ea": ["HFC227ea", 1],
    "HFC23": ["HFC23", 1],
    "HFC245fa": ["HFC245fa", 1],
    "HFC32": ["HFC32", 1],
    "HFC4310mee": ["HFC4310mee", 1],
    "C2F6": ["C2F6", 1],
    "C6F14": ["C6F14", 1],
    "CF4": ["CF4", 1],
    "SF6": ["SF6", 1],
    "NOx": ["NOx", 1],
    "CO": ["CO", 1],
    "NH3": ["NH3", 1],
    "BMB_AEROS_BC": ["BMB_AEROS_BC", 1],
    "BMB_AEROS_OC": ["BMB_AEROS_OC", 1],
    "BC": ["BC", 1],
    "OC": ["OC", 1],
}
# Halon1212, CH3Cl


class COMMONSFILEWRITER:  # pylint: disable=too-few-public-methods
    """
    Class to write scenariofiles:
    """

    def __init__(self, gaspam_file, nystart=2015, nyend=2101):
        self.components = []
        self.units = []
        self.concunits = []

        self.initialize_units_comps(gaspam_file)
        self.years = np.arange(
            nystart, nyend
        )  # Temporary default values, is updated later

    def initialize_units_comps(self, gasfile):
        """
        Get the list of gas components and units
        from the gases file:
        """
        if isinstance(gasfile, pd.DataFrame):
            self.initialize_units_comps_from_pd(gasfile)
            return
        with open(gasfile, "r", encoding="ascii") as txt_rcpfile:
            gasreader = csv.reader(txt_rcpfile, delimiter="\t")
            next(gasreader)
            for row in gasreader:
                if row[1] == "X":
                    continue

                component = row[0]
                unit = row[1]

                if component == "N2O" and unit == "Tg_N":
                    # in openscm-units, to get the mass of nitrogen, have
                    # to use the unit "Tg N2ON" (converting to "Tg N" just
                    # converts using the mass fraction of a single nitrogen
                    # atom, admittedly this isn't immediately obvious and
                    # arguably is a bug in openscm-units)
                    unit = "TgN2ON"
                elif "_" in unit:
                    unit = unit.replace("_", "")
                else:
                    comp_str = component.replace("-", "").replace("BMB_AEROS_", "")
                    unit = f"{unit}{comp_str}"

                unit = f"{unit} / yr"

                self.components.append(component)
                self.units.append(unit)
                self.concunits.append(row[2])

        self.components.insert(1, "CO2_lu")
        self.units.insert(1, "PgC / yr")
        self.concunits.insert(1, "ppm")

    def initialize_units_comps_from_pd(self, gasfile):
        """
        Intialize units and compontets if pandas of gasfile is sent deal with units
        """
        gasdata_mod = gasfile.copy()
        # Might cause problems later, hopefully not
        gasdata_mod.replace("Tg_N", "TgN2ON", inplace=True)
        gasdata_mod["CONC_UNIT"] = gasdata_mod["CONC_UNIT"].replace("_", "")
        # Not very robuste...:
        gasdata_mod.loc["BMB_AEROS_OC", "EM_UNIT"] = gasdata_mod["EM_UNIT"][
            "BMB_AEROS_OC"
        ].replace("Tg", "TgOC")
        gasdata_mod.loc["BMB_AEROS_BC", "EM_UNIT"] = gasdata_mod["EM_UNIT"][
            "BMB_AEROS_BC"
        ].replace("Tg", "TgBC")
        self.components = gasdata_mod.index.values.tolist()
        units = gasdata_mod["EM_UNIT"].values.tolist()
        self.concunits = gasdata_mod["CONC_UNIT"].values.tolist()

        self.units = [f"{unit} / yr" for unit in units]
        self.components.insert(1, "CO2_lu")
        self.units.insert(1, "PgC / yr")
        self.concunits.insert(1, "ppm")
