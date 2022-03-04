"""
Perturbation related methods
"""
import numpy as np
import pandas as pd


def calculate_hemispheric_forcing(tracer, q, forc_nh, forc_sh):
    """
    Calculate hemispheric forcing per tracer
    """
    if tracer in ("SO2", "SO4_IND"):
        forc_nh = forc_nh + q * 1.6
        forc_sh = forc_sh + q * 0.4
    elif tracer == "TROP_O3":
        # 1.29+0.74 != 2, update/check
        forc_nh = forc_nh + q * 1.29  # 1.3
        forc_sh = forc_sh + q * 0.74  # 0.7
    else:
        forc_nh = forc_nh + q
        forc_sh = forc_sh + q
    return forc_nh, forc_sh


def perturb_emissions(perturbation_file, emissions_df):
    """
    Add emission perturbations to emissions_df
    """
    pert_df = pd.read_csv(perturbation_file, index_col=None)
    for row in pert_df.itertuples(index=True, name="Pandas"):
        tracer = row.component
        if row.component == "CO2" and "CO2" not in emissions_df:
            tracer = "CO2_FF"
        emissions_df[tracer][row.year] = emissions_df[tracer][row.year] + row.emission
    # might not need to return, as the pandas should change itself
    return emissions_df


class ForcingPerturbation:
    """
    Class to handle  forcing perturbations
    """

    def __init__(self, perturbation_file, year0):
        self.perturb_raw = pd.read_csv(perturbation_file)
        self.years = np.unique(self.perturb_raw["year"].values)
        self.compounds = pd.unique(self.perturb_raw["component"].values)
        self.year0 = year0

    def check_if_year_in_pert(self, year):
        """
        Check if year has perturbations
        """
        return bool(year in self.years)

    def check_if_compound_in_pert(self, compound):
        """
        Check if compound is perturbed at some point
        """
        return bool(compound in self.compounds)

    def add_forcing_pert(
        self, totforc, forc_nh, forc_sh, forc, yr
    ):  # pylint: disable=too-many-arguments
        """
        Add forcing perturbations to precalculated forcing
        """
        perturb_here = self.perturb_raw[self.perturb_raw.year == yr]
        for row in perturb_here.itertuples(index=True, name="Pandas"):
            value = row.forcing
            tracer = row.component
            totforc = totforc + value
            forc_nh, forc_sh = calculate_hemispheric_forcing(
                tracer, value, forc_nh, forc_sh
            )
            forc[tracer][yr - self.year0] = forc[tracer][yr - self.year0] + value
        return totforc, forc_nh, forc_sh, forc
