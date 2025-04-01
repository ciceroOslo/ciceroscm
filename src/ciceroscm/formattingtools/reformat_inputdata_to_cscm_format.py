"""
Common handler and converter of scenariodat for cicero
Class will be inherited by subversions to get specific
handling
"""

import csv
import logging
import os

import numpy as np
import openscm_units
import pandas as pd

LOGGER = logging.getLogger(__name__)


def _get_unique_index_values(idf, index_col, assert_all_same=True):
    """
    Get unique values in index column from a dataframe

    Parameters
    ----------
    idf : :obj:`pd.DataFrame`
        Dataframe to get index values from

    index_col : str
        Column in index to get the values for

    assert_all_same : bool
        Should we assert that all the values are the same before
        returning? If True, only a single value is returned. If
        False, a list is returned.

    Returns
    -------
    str, list
        Values found, either a string or a list depending on
        ``assert_all_same``.

    Raises
    ------
    AssertionError
        ``assert_all_same`` is True and there's more than one
        unique value.
    """
    out = idf.index.get_level_values(index_col).unique().tolist()
    if assert_all_same:
        if len(out) > 1:
            raise AssertionError(out)

        return out[0]

    return out


def _unit_conv_factor(unit, cicero_unit):
    """
    Convert cicero units that can't be automatically dealt with
    in openscm_units
    """
    with openscm_units.unit_registry.context("NOx_conversions"):

        if cicero_unit.startswith("GgH1"):
            conv_factor = (
                openscm_units.unit_registry(unit)
                .to(cicero_unit.replace("GgH1", "GgHalon1"))
                .magnitude
            )
        elif cicero_unit.startswith("GgH2"):
            conv_factor = (
                openscm_units.unit_registry(unit)
                .to(cicero_unit.replace("GgH2", "GgHalon2"))
                .magnitude
            )
        else:
            conv_factor = openscm_units.unit_registry(unit).to(cicero_unit).magnitude

    return conv_factor


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


class COMMONSFILEWRITER:
    """
    Class to write scenariofiles:
    """

    def __init__(self, udir, nystart=2015, nyend=2101):
        self.components = []
        self.units = []
        self.concunits = []

        self.initialize_units_comps(os.path.join(udir, "gases_v1RCMIP.txt"))
        self.years = np.arange(
            nystart, nyend
        )  # Temporary default values, is updated later
        self.udir = udir

    def initialize_units_comps(self, gasfile):
        """
        Get the list of gas components and units
        from the gases file:
        """
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

    def get_unit_convfactor(self, comp, scenarioframe):
        """
        Get unit and conversion factor for component
        """
        # Find the unit and the original unit
        cicero_unit = self.units[self.components.index(comp)]
        unit = _get_unique_index_values(
            scenarioframe[
                scenarioframe.index.get_level_values("variable")
                == f"Emissions|{cicero_comp_dict[comp][0]}"
            ],
            "unit",
        )

        return _unit_conv_factor(unit, cicero_unit)

    def transform_scenarioframe(self, scenarioframe):
        """
        Get rid of multiindex and interpolate scenarioframe
        """
        if (
            _get_unique_index_values(scenarioframe, "region") != "World"
        ):  # pragma: no cover
            raise NotImplementedError()  # emergency valve

        scenarioframe = scenarioframe.reset_index(
            ("model", "region", "scenario", "unit"), drop=True
        )
        years = scenarioframe.columns

        if not isinstance(years[0], np.int64):
            yearsint = [np.int64(d.year) for d in years]
            scenarioframe.rename(
                lambda d: np.int64(d.year), axis="columns", inplace=True
            )
        else:
            yearsint = years

        self.years = np.arange(yearsint[0], yearsint[-1] + 1)
        for year in self.years:
            if year not in scenarioframe.columns:
                scenarioframe[year] = np.nan

        scenarioframe = scenarioframe.reindex(sorted(scenarioframe.columns), axis=1)
        interpol = scenarioframe.interpolate(axis=1)

        return interpol

    def make_printoutframe(self, scenarioframe, ssp245data):
        """
        Take scenarioframe and convert to cicero format
        emissions data to write to file or pass as dataframe
        """
        logging.getLogger("pyam").setLevel(logging.ERROR)
        avail_comps = [
            c.replace("Emissions|", "")
            for c in _get_unique_index_values(
                scenarioframe, "variable", assert_all_same=False
            )
        ]
        ciceroscm_comps = [v[0] for v in cicero_comp_dict.values()]
        not_used_comps = set(avail_comps) - set(ciceroscm_comps)
        if not_used_comps:
            LOGGER.warning("%s not used by CICERO-SCM", not_used_comps)

        interpol = self.transform_scenarioframe(scenarioframe)
        printout_frame = pd.DataFrame(columns=self.components, index=self.years)

        # Setting conversion factors for components with data from scenarioframe
        for comp in self.components:
            if cicero_comp_dict[comp][0] in avail_comps:
                convfactor = self.get_unit_convfactor(comp, scenarioframe)
                if (
                    cicero_comp_dict[comp][0] in ("BC", "OC")
                    and f"BMB_AEROS_{cicero_comp_dict[comp][0]}" not in avail_comps
                ):
                    printout_frame[comp] = (
                        interpol.T[f"Emissions|{cicero_comp_dict[comp][0]}"]
                        * convfactor
                    ).to_numpy() - ssp245data[
                        f"BMB_AEROS_{cicero_comp_dict[comp][0]}"
                    ].loc[
                        str(self.years[0]) : str(self.years[-1])
                    ].to_numpy().astype(
                        np.float
                    )
                else:
                    printout_frame[comp] = (
                        interpol.T[f"Emissions|{cicero_comp_dict[comp][0]}"]
                        * convfactor
                    )
            else:
                LOGGER.warning("No %s data available, using ssp245", comp)
                printout_frame[comp] = (
                    ssp245data[comp]
                    .loc[str(self.years[0]) : str(self.years[-1])]
                    .to_numpy()
                )

        printout_frame = printout_frame.astype(float)
        return printout_frame


def _read_ssp245_em(ssp245_em_file):
    """
    Get default data from ssp245_RCMIP
    """
    ssp245df = (
        pd.read_csv(ssp245_em_file, delimiter="\t", index_col=0, skiprows=[1, 2, 3])
        .rename(columns=lambda x: x.strip())
        .rename(columns={"CO2 .1": "CO2_AFOLU"})
        .rename(columns={"CO2": "CO2_FF"})
        .astype(float)
    )

    return ssp245df


class SCENARIODATAGETTER(COMMONSFILEWRITER):
    """
    Class to write scenariofiles:
    """

    def __init__(self, udir, nystart=2015, nyend=2101):
        """
        Intialise scenario data getter
        """
        super().__init__(udir, nystart, nyend)
        self.ssp245data = _read_ssp245_em(os.path.join(udir, "ssp245_em_RCMIP.txt"))

    def get_scenario_data(self, scenarioframe):
        """
        Get printoutframe, and adding ssp245 data to
        get a frame for running
        """
        printout_frame = self.make_printoutframe(scenarioframe, self.ssp245data)
        printout_frame.rename(
            columns={"CO2_lu": "CO2_AFOLU", "CO2": "CO2_FF"}, inplace=True
        )
        final_frame = pd.concat(
            [
                self.ssp245data.iloc[: (self.years[0] - 1750)]
                .astype(float)
                .rename(columns={"Component": "Index"}),
                printout_frame,
            ]
        )
        return final_frame
