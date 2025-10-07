import pandas as pd

from .carbon_cycle_abstract import AbstractCarbonCycleModel
from .carmen.src.carmen import CarbonCycle
from .._utils import update_pam_if_numeric
# from .carmen import CarbonCycle
# from .carmen.utils import load_esm_data
# from .carmen.constants import SCEN_DIR

DT_MODEL = 1/8
DT_MODEL_OCEAN = 1/8



class CarbonCycleModel(AbstractCarbonCycleModel):
    """
    CarbonCycleModel coupling the CARMEN model.
    """

    carbon_cycle_model_required_pamset = {
        "esm_to_emulate": "NorESM2-LM",
        "gpp_t_l": 0.14199847,
        "gpp_t_e": -0.11203268,
        "gpp_c_l": -0.09139857,
        "gpp_c_half": 384.32434502,
        "gpp_c_e": 0.0698166,
        "gpp_hyst": 0.01289971,
        "gpp_c_tan": 0.0,
        "gpp_fast": 0.91751023,
        "gpp_slow": 0.92258006,
        "gpp_c_tan2": 0.0,
        "lit_t_l": 1.143e-05,
        "lit_t_e": 0.01792216,
        "lit_c_l": -0.61640124,
        "lit_c_half": 256.69192703,
        "lit_c_e": -0.23731487,
        "lit_hyst": 0.0167247,
        "lit_c_tan": 0.0,
        "lit_fast": 0.92847622,
        "lit_slow": 0.91882796,
        "lit_c_tan2": 0.0,
        "vres_t_l": 8.85e-06,
        "vres_t_e": -0.07693331,
        "vres_c_l": -0.47443673,
        "vres_c_half": 808.47364366,
        "vres_c_e": -0.37496243,
        "vres_hyst": -0.04767878,
        "vres_c_tan": 0.0,
        "vres_fast": 0.90821573,
        "vres_slow": 0.84546906,
        "vres_c_tan2": 0.0,
        "sres_t_l": 8.76e-06,
        "sres_t_e": 0.06184566,
        "sres_c_l": 1.43982004,
        "sres_c_half": 92.12981868,
        "sres_c_e": -0.24272656,
        "sres_hyst": 0.04094068,
        "sres_c_tan": 0.0,
        "sres_fast": 0.89815371,
        "sres_slow": 0.85285196,
        "sres_c_tan2": 0.0,
        "npp_t_l": 0.20958114,
        "npp_t_e": -0.09892429,
        "npp_c_l": -0.19742283,
        "npp_c_half": 281.66916948,
        "npp_c_e": 0.05305683,
        "npp_hyst": 0.07758387,
        "npp_c_tan": 0.0,
        "npp_fast": 0.94904171,
        "npp_slow": 0.71402205,
        "npp_c_tan2": 0.0,
        "docn": 41.80816876,
        "docnfac": 9.50761817,
        "ocntemp": 0.0, 
        "docntemp": 0.12145903,
    }

    def __init__(self, pamset_emiconc, pamset_carbon=None):
        """
        Initialize the Carbon Cycle Model.

        Parameters
        ----------
        pamset_emiconc : dict
            Parameter set from the concentrations emission handler, containing:
            - idtm: Number of subyearly timesteps (e.g., 24 for monthly steps).
            - nystart: Start year of the simulation.
            - nyend: End year of the simulation.
        pamset_carbon : dict
            Optional carbon specific parameter set
            TODO: explain what the model is expecting
        """
        if pamset_carbon is None:
            pamset_carbon = self.carbon_cycle_model_required_pamset

        # Pass the necessary variable to the parent class
        super().__init__(pamset_emiconc, pamset_carbon=pamset_carbon)

        # Create and store an object implementing CARMEN
        carmen_instance = CarbonCycle(
            # if you want to use ESM scenario data
            # {"model": "NorESM2-LM", "scenario": "ssp126"},
            # if you want to run without any pre-loaded scenario
            {"model": self.pamset.get("esm_to_emulate"), "initial_year": self.pamset["nystart"], "final_year": self.pamset["nyend"]},
            1/self.pamset["idtm"],
            1/self.pamset["idtm"],
            npp_flag=True,
            **pamset_carbon,
        )
        self.carmen = carmen_instance


    def co2em2conc(self, _yr, em_co2_common, dtemp=0):  # pylint: disable=unused-argument
        """
        Calculate CO2 concentrations from emissions.

        Parameters
        ----------
        yr : int
            Year for which to calculate.
        em_co2_common : float
            Total CO2 emissions (GtC/yr).
        dtemp : float
            Global mean temperature difference from pre-industrial (K).

        Returns
        -------
        float
            Updated atmospheric CO2 concentration (ppm).
        """

        # TODO: the asumption here is taht this function is going to be run everytime we want a new timestep
        # TODO: but with carmen you can't really go back, which from the description of the function I'm getting the
        # TODO: sense it is not the case here.

        new_input = {
            # emis is technically the rate of yearly emissions
            "emis": em_co2_common,
            "dtocn": dtemp, # TODO: need ocean temperature passed
            "dtglb": dtemp,
        }

        # Run model one step into the future
        # CICERO-SCM actually expect this function to run a whole year, so run the appropriate number of CARMEN steps
        for _ in range(int(1/self.carmen.dtpred)):
            self.carmen.run_one_step(new_input)

        return self.carmen.catm[self.carmen.current_step]


    def get_carbon_cycle_output(
        self, years, conc_run=False, conc_series=None, dtemp_series=None
    ):
        # Report the carbon cycle quantities in year format
        self.carmen.interpolate_results(years)
        df_carbon = pd.DataFrame(
            data={
                "Biosphere carbon flux": self.carmen.land.npp - self.carmen.land.sres,
                "Biosphere carbon pool": self.carmen.land.cveg + self.carmen.land.sres,
                "Vegetation carbon pool": self.carmen.land.cveg,
                "Soil carbon pool": self.carmen.land.csoil,
                "Ocean carbon pool": self.carmen.ocean.carbon_increase,
                "Net primary production": self.carmen.land.npp,
                "Litterfall": self.carmen.land.lit,
                "Soil respiration": self.carmen.land.sres,
                "Ocean carbon flux": self.carmen.ocean.oflux,
            },
            index=years,
        )
        return df_carbon


    def reset_co2_hold(self, pamset_carbon=None):
        if pamset_carbon is not None:
            self.pamset = update_pam_if_numeric(
                self.pamset,
                pamset_new=pamset_carbon,
                can_change=list(self.carbon_cycle_model_required_pamset.keys())+["nystart", "nyend"],
            )


        carmen_instance = CarbonCycle(
            # if you want to use ESM scenario data
            # {"model": "NorESM2-LM", "scenario": "ssp126"},
            # if you want to run without any pre-loaded scenario
            {"model": self.pamset.get("esm_to_emulate"), "initial_year": self.pamset["nystart"], "final_year": self.pamset["nyend"]},
            1/self.pamset.get("idtm"),
            1/self.pamset.get("idtm"),
            npp_flag=True,
            **self.pamset,
        )
        self.carmen = carmen_instance


    def _get_co2_hold_values(self):
        """
        Get a dictionary of state variables for the carbon cycle that allows you
        to reset your model state, typically pool sizes, chemical composition
        current year etc...

        Make sure to overwrite this is if your model actually has state dependence
        """
        current_state = {
            # land box
            "cveg": self.carmen.land.cveg,
            "csoil": self.carmen.land.csoil,
            "npp": self.carmen.land.npp,
            "lit": self.carmen.land.lit,
            "sres": self.carmen.land.sres,
            "fcva": self.carmen.land.fcva,
            "fcsa": self.carmen.land.fcsa,
            "fcvs": self.carmen.land.fcvs,
            "peak_t": self.carmen.land.peak_t,
            "old_t": self.carmen.land.old_t,
            "step_index_land": self.carmen.land.timestep_ind,

            # Ocean box
            "cocean": self.carmen.ocean.carbon_increase,
            "oflux": self.carmen.ocean.oflux,
            "step_index_ocean": self.carmen.ocean.timestep_ind,

            # General carbon cycle box
            "catm": self.carmen.catm,
            "emis": self.carmen.emis,
            "step_index_general": self.carmen.current_step,
        }

        return current_state

    def _set_co2_hold_values(self, hold_dict):
        """
        Set state variables for the carbon cycle from hold_dict dictionary
        thereby resetting your model state, typically pool sizes, chemical composition
        current year etc...

        Make sure to overwrite this is if your model actually has state dependence
        """
        self.carmen.land.cveg = hold_dict["cveg"]
        self.carmen.land.csoil = hold_dict["csoil"]
        self.carmen.land.npp = hold_dict["npp"]
        self.carmen.land.lit = hold_dict["lit"]
        self.carmen.land.sres = hold_dict["sres"]
        self.carmen.land.fcva = hold_dict["fcva"]
        self.carmen.land.fcsa = hold_dict["fcsa"]
        self.carmen.land.fcvs = hold_dict["fcvs"]
        self.carmen.land.peak_t = hold_dict["peak_t"]
        self.carmen.land.old_t = hold_dict["old_t"]
        self.carmen.land.timestep_ind = hold_dict["step_index_land"]

        # Ocean box
        self.carmen.ocean.carbon_increase = hold_dict["cocean"]
        self.carmen.ocean.oflux = hold_dict["oflux"]
        self.carmen.ocean.timestep_ind = hold_dict["step_index_ocean"]

        # General carbon cycle box
        self.carmen.catm = hold_dict["catm"]
        self.carmen.emis = hold_dict["emis"]
        self.carmen.current_step = hold_dict["step_index_general"]
