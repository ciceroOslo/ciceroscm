"""
Energy budget upwelling diffusion model
"""

import logging

import numpy as np
from scipy.linalg import solve_banded

from .._utils import cut_and_check_pamset
from ..constants import (
    C_OHC,
    CONV_FAC_UDM,
    CP_UDM,
    DAYS_PER_YEAR,
    DZ,
    INDEX_700M,
    OCEAN_AREA_HEMISPHERE,
    RHO_OHC,
    RHO_SEAWATER,
    SECONDS_PER_DAY,
)
from .abstract_thermal_model import AbstractThermalModel

LOGGER = logging.getLogger(__name__)


def _band(a_array, b_array, c_array, d_array):
    """
    Calculate band

    Parameters
    ----------
    a_array : np.ndarray
              a_array through ocean layers
    b_array : np.ndarray
              b_array through ocean layers
    c_array : np.ndarray
              c_array through ocean layers
    d_array : np.ndarray
               d_array through ocean layers

    Returns
    -------
    np.ndarray
             band value through ocean layers
    """
    return solve_banded((1, 1), np.array([c_array, b_array, a_array]), d_array)


def check_pamset(pamset):
    """

    Check that parameterset has necessary values for run
    Otherwise set to default values which are defined here

    Parameters
    ----------
    pamset : dict
          Dictionary of parameters to define the physics
          of the run

    Returns
    -------
    dict
        Updated pamset with default values used where necessary
    """
    if pamset is None:
        pamset = {}
    required = {
        "rlamdo": 15.0,
        "akapa": 0.66,
        "cpi": 0.21,
        "W": 2.2,
        "beto": 6.9,
        "threstemp": 7.0,
        "lambda": 0.61,
        "mixed": 107.0,
        "foan": 0.61,
        "foas": 0.81,
        "ebbeta": 0.0,
        "fnso": 0.7531,
        "lm": 40,
        "ldtime": 12,
        "ocean_efficacy": 1.0,
    }
    pamset = cut_and_check_pamset(required, pamset, cut_warnings=True)
    pamset["rakapa"] = 1.0e-4 * pamset["akapa"]
    pamset["rlamda"] = 1.0 / pamset["lambda"]
    pamset["dt"] = 1 / pamset["ldtime"] * SECONDS_PER_DAY * DAYS_PER_YEAR
    pamset["c1"] = RHO_OHC * CP_UDM * CONV_FAC_UDM * DZ * SECONDS_PER_DAY
    pamset["fnx"] = (
        pamset["rlamda"] + pamset["foan"] * pamset["rlamdo"] + pamset["ebbeta"]
    )
    pamset["fsx"] = (
        pamset["rlamda"] + pamset["foas"] * pamset["rlamdo"] + pamset["ebbeta"]
    )
    return pamset


class UpwellingDiffusionModel(
    AbstractThermalModel
):  # pylint: disable=too-few-public-methods,too-many-instance-attributes
    """
    Upwelling Diffusion Model for ocean thermal dynamics.

    A 40-layer ocean thermal model that simulates heat diffusion
    into the ocean from the surface using upwelling and diffusion processes.
    This model provides vertical ocean temperature structure and
    accounts for hemisphere-specific ocean dynamics.

    The model inherits from AbstractThermalModel, enabling integration
    with the thermal model factory system and interchangeability
    with other thermal models.

    Parameters
    ----------
    pamset : dict, optional
        Parameter set containing model configuration. If None, default
        parameters are used. Key parameters include:

        - lambda : float
            Climate feedback parameter (W/m²/K)
        - mixed : float
            Mixed layer depth (m)
        - akapa : float
            Thermal diffusion coefficient
        - foan, foas : float
            Northern/Southern hemisphere ocean fractions
        - lm : int
            Number of ocean layers (default 40)
        - ldtime : int
            Number of time steps per year

    Attributes
    ----------
    thermal_model_required_pamset : dict
        Dictionary of required parameter names with default values
    pamset : dict
        Current parameter set used by the model, including derived parameters
        such as rakapa, rlamda, dt, c1, fnx, fsx
    dz : np.ndarray
        Layer thicknesses (m) for each ocean layer. First element is mixed
        layer depth, remaining are standard 100m layers
    tn : np.ndarray
        Northern hemisphere ocean temperatures by layer (K)
    ts : np.ndarray
        Southern hemisphere ocean temperatures by layer (K)
    varrying : dict
        Dictionary containing time-varying calculation parameters including
        matrix coefficients (acoeffn, bcoeffn, ccoeffn, acoeffs, bcoeffs, ccoeffs)
        and derived thermal terms (dtrm1n, dtrm2n, etc.)
    prev_values : dict
        Dictionary of values to preserve between yearly calculations.
        Contains 'fn', 'fs' (hemisphere forcings), and 'dtemp' (global temperature)
    dtempprev : float
        Previous global temperature change (K), used for continuity
    gamn : float
        Northern hemisphere thermal response coefficient
    gams : float
        Southern hemisphere thermal response coefficient

    See Also
    --------
    AbstractThermalModel : Base class for thermal models
    TwoLayerOceanModel : Alternative simplified thermal model
    """

    @property
    def thermal_model_required_pamset(self):
        """
        Required parameter set for upwelling diffusion model

        Returns
        -------
        dict
            Dictionary of required parameters with default values
        """
        return {
            "rlamdo": 15.0,
            "akapa": 0.66,
            "cpi": 0.21,
            "W": 2.2,
            "beto": 6.9,
            "threstemp": 7.0,
            "lambda": 0.61,
            "mixed": 107.0,
            "foan": 0.61,
            "foas": 0.81,
            "ebbeta": 0.0,
            "fnso": 0.7531,
            "lm": 40,
            "ldtime": 12,
            "ocean_efficacy": 1.0,
        }

    @classmethod
    def get_thermal_model_required_pamset(cls):
        """
        Get the required parameter set for the Upwelling Diffusion Model.

        Returns the dictionary of parameter names with their default values
        that are required for proper operation of the UDM thermal model.
        This method enables the AbstractThermalModel interface compatibility
        and supports the factory system's parameter validation.

        Returns
        -------
        dict
            Dictionary of required parameter names with default values.
            These parameters must be provided (or have defaults) for the
            model to function correctly. The dictionary includes physical
            parameters such as:

            - Climate feedback parameters (lambda, rlamdo)
            - Ocean structure (mixed, lm, ldtime)
            - Thermal properties (akapa, W, cpi)
            - Ocean fractions (foan, foas)
            - Other physics constants (beto, ebbeta, fnso, etc.)

        Notes
        -----
        This class method is required by the AbstractThermalModel interface
        and enables parameter validation in the factory system. The returned
        parameter dictionary corresponds to the `thermal_model_required_pamset`
        property.

        Examples
        --------
        >>> required_params = UpwellingDiffusionModel.get_thermal_model_required_pamset()
        >>> print(f"UDM requires {len(required_params)} parameters")
        >>> print(sorted(required_params.keys()))

        See Also
        --------
        thermal_model_required_pamset : Property containing the same information
        AbstractThermalModel.get_thermal_model_required_pamset : Parent method
        """
        # Create a temporary instance to access the property
        temp_instance = cls.__new__(cls)
        return temp_instance.thermal_model_required_pamset

    def __init__(self, pamset=None):
        """
        Initialize the Upwelling Diffusion Model.

        Sets up the ocean layer structure, validates parameters, and calculates
        derived parameters needed for the thermal dynamics simulation. The model
        can be initialized with custom parameters or will use defaults if none
        are provided.

        This constructor maintains backward compatibility with the original UDM
        interface while adding AbstractThermalModel inheritance and validation.

        Parameters
        ----------
        pamset : dict, optional
            Dictionary of physical parameters to configure the model instance.
            If None, default parameter values are used. Required parameters
            are validated through the parent AbstractThermalModel class.

            Key parameters include:
            - lambda : Climate feedback parameter (W/m²/K)
            - mixed : Mixed layer depth (m)
            - akapa : Thermal diffusion coefficient
            - foan, foas : Ocean fractions for N/S hemispheres
            - lm : Number of ocean layers (default 40)
            - ldtime : Time steps per year

        Notes
        -----
        The constructor automatically calculates several derived parameters:
        - rakapa : Scaled thermal diffusion (1e-4 * akapa)
        - rlamda : Inverse climate feedback (1/lambda)
        - dt : Time step in seconds
        - c1 : Heat capacity conversion factor
        - fnx, fsx : Hemisphere-specific thermal response factors

        The ocean layer structure is set up with variable thickness layers,
        starting with the mixed layer depth and continuing with 100m layers
        down to the deep ocean.

        Examples
        --------
        >>> # Initialize with default parameters
        >>> udm = UpwellingDiffusionModel()

        >>> # Initialize with custom parameters
        >>> params = {'lambda': 0.6, 'mixed': 120.0, 'akapa': 0.8}
        >>> udm = UpwellingDiffusionModel(params)

        >>> # Use through factory system
        >>> from ciceroscm.component_factory_functions import create_thermal_model
        >>> thermal_class = create_thermal_model('default')
        >>> udm = thermal_class(params)
        """
        # Call parent constructor which handles parameter validation
        super().__init__(pamset)

        # Add derived parameters that check_pamset used to calculate
        self.pamset["rakapa"] = 1.0e-4 * self.pamset["akapa"]
        self.pamset["rlamda"] = 1.0 / self.pamset["lambda"]
        self.pamset["dt"] = 1 / self.pamset["ldtime"] * SECONDS_PER_DAY * DAYS_PER_YEAR
        self.pamset["c1"] = RHO_OHC * CP_UDM * CONV_FAC_UDM * DZ * SECONDS_PER_DAY
        self.pamset["fnx"] = (
            self.pamset["rlamda"]
            + self.pamset["foan"] * self.pamset["rlamdo"]
            + self.pamset["ebbeta"]
        )
        self.pamset["fsx"] = (
            self.pamset["rlamda"]
            + self.pamset["foas"] * self.pamset["rlamdo"]
            + self.pamset["ebbeta"]
        )

        # Setting up dz height difference between ocean layers
        self.dz = np.ones(self.pamset["lm"]) * DZ
        self.dz[0] = self.pamset["mixed"]
        self.varrying = {}
        self.setup_ebud()

        # Intialising temperature values
        self.tn = np.zeros(self.pamset["lm"])
        self.ts = np.zeros(self.pamset["lm"])
        # Dict of values to keep from one year to the next
        self.prev_values = {
            "fn": 0.0,
            "fs": 0.0,
            "dtemp": 0.0,
        }

        self.dtempprev = 0.0

    def get_gam_and_fro_factor_ns(self, northern_hemisphere):
        """
        Get correct gam and fro variables

        Get correct gam and fro variables depending on
        whether Northern or Southern hemisphere is
        considered

        Parameters
        ----------
        northern_hemisphere : bool
                           Whether northern or southern hemisphere
                           is being considered

        Returns
        -------
        float
             The correct gam_fro_factor
        """
        blm = self.pamset["ebbeta"] / self.pamset["rlamdo"]
        if northern_hemisphere:
            gam1 = self.gamn
            gam2 = self.gams
            fro1 = self.pamset["foan"]
        else:
            gam1 = self.gams
            gam2 = self.gamn
            fro1 = self.pamset["foas"]
        factor = (
            self.pamset["rlamdo"]
            * (1.0 - fro1 * gam2 / (gam2 * gam1 - blm * blm))
            * self.pamset["dt"]
            / (self.pamset["c1"] * self.dz[0])
        )
        return factor

    def coeff(self, wcfac, gam_fro_fac):
        """
        Calculate a, b c coefficient arrays for hemisphere

        This method has been changed since the fortran version
        to have the coefficients a be the coefficients for the
        layer underneath, b be the coefficient for current layer
        and c be the coefficient for the layer above in such a
        way that they can represent a banded matrix and solved for
        as such

        Parameters
        ----------
        wcfac : float
             wc factor is W per fraction of the year
        gam_fro_fac : float
                      Is the factor combo of gam and fro
        Returns
        -------
        list
            Containing the a, b and c coefficients over the ocean layers
        """
        lm = self.pamset["lm"]
        a = np.zeros(lm)
        b = np.zeros(lm)
        c = np.zeros(lm)
        rakapafac = 2 * self.pamset["rakapa"] * self.pamset["dt"]

        c[1] = -rakapafac / (
            self.dz[0] * (0.0 * self.dz[0] + self.dz[1])
        )  # Can the 0.*dz(0) term be dropped here?
        b[0] = 1.0 - c[1] + gam_fro_fac - wcfac / self.dz[0]
        a[0] = -rakapafac / (self.dz[1] ** 2) + wcfac / self.dz[1]
        a[1 : lm - 1] = -rakapafac / (self.dz[2:] * (self.dz[1 : lm - 1] + self.dz[2:]))
        c[2:] = (
            -rakapafac / (self.dz[1 : lm - 1] * (self.dz[1 : lm - 1] + self.dz[2:]))
            - wcfac / self.dz[1 : lm - 1]
        )
        b[1 : lm - 1] = 1.0 - a[: lm - 2] - c[2:]
        b[lm - 1] = (
            1.0 - a[lm - 2] + wcfac / self.dz[lm - 1]
        )  # Her var det brukt i selvom vi var utenfor loekka, litt uklart hva som er ment...
        return a, b, c

    def setup_ebud2(self, temp_1n, temp_1s):
        """
        Set up coefficients and more for the two hemispheres

        Set up coefficients and various parameters
        for the two hemispheres. To be redone every timestep

        Parameters
        ----------
        temp_1n : float
               Northern hemisphere surface temperature
        temp_1s : float
               Southern hemisphere surface temperature
        """
        ocean_efficacy = self.pamset.get("ocean_efficacy", 1.0)

        # Northern hemisphere:
        if self.pamset["threstemp"] == 0:  # pylint: disable=compare-to-zero
            wcfac = (
                self.pamset["W"] / (SECONDS_PER_DAY * DAYS_PER_YEAR) * self.pamset["dt"]
            )
        else:
            wcfac = (
                self.pamset["W"]
                / (SECONDS_PER_DAY * DAYS_PER_YEAR)
                * (1 - 0.3 * temp_1n / self.pamset["threstemp"])
                * self.pamset["dt"]
            )

        ##NEW efficacy logic
        # 1. Get the standard, physical coefficients.
        a_phys, b_phys, c_phys = self.coeff(wcfac, self.get_gam_and_fro_factor_ns(True))

        # 2. Create the EFFECTIVE coefficients for the solver as a copy.
        a_eff, b_eff, c_eff = a_phys.copy(), b_phys.copy(), c_phys.copy()

        # 3. Apply the ASYMMETRIC scaling to the mixed-layer/deep-ocean interface.
        if ocean_efficacy != 1.0:
            # Scale the influence of the deep ocean ON the mixed layer's budget (Row 0).
            c_eff[1] = c_phys[1] * ocean_efficacy

            # Correct the mixed layer's main diagonal for the scaled flux.
            # The influence of the mixed layer ON the deep ocean (a_eff[0]) remains UNCHANGED.
            b_eff[0] = b_phys[0] + c_phys[1] - c_eff[1]

        # 4. Store the effective coefficients for the NH solver.
        self.varrying["acoeffn"], self.varrying["bcoeffn"], self.varrying["ccoeffn"] = (
            a_eff,
            b_eff,
            c_eff,
        )

        # Pre-calculate RHS factors as before.
        self.varrying["dtrm1n"] = (
            1.0
            - self.pamset["cpi"] * wcfac / self.dz[0]
            - self.pamset["beto"] * self.pamset["dt"] / (self.pamset["c1"] * self.dz[0])
        )
        self.varrying["dtmnl2"] = (
            wcfac * self.pamset["cpi"] / self.dz[self.pamset["lm"] - 1]
        )

        # Southern hemisphere:
        if self.pamset["threstemp"] == 0:  # pylint: disable=compare-to-zero
            wcfac = (
                self.pamset["W"] / (SECONDS_PER_DAY * DAYS_PER_YEAR) * self.pamset["dt"]
            )
        else:
            wcfac = (
                self.pamset["W"]
                / (SECONDS_PER_DAY * DAYS_PER_YEAR)
                * (1 - 0.3 * temp_1s / self.pamset["threstemp"])
                * self.pamset["dt"]
            )

        a_phys, b_phys, c_phys = self.coeff(
            wcfac, self.get_gam_and_fro_factor_ns(False)
        )
        a_eff, b_eff, c_eff = a_phys.copy(), b_phys.copy(), c_phys.copy()

        if ocean_efficacy != 1.0:
            c_eff[1] = c_phys[1] * ocean_efficacy
            b_eff[0] = b_phys[0] + c_phys[1] - c_eff[1]

        self.varrying["acoeffs"], self.varrying["bcoeffs"], self.varrying["ccoeffs"] = (
            a_eff,
            b_eff,
            c_eff,
        )

        self.varrying["dtrm1s"] = (
            1.0
            - self.pamset["cpi"] * wcfac / self.dz[0]
            - self.pamset["fnso"]
            * self.pamset["beto"]
            * self.pamset["dt"]
            / (self.pamset["c1"] * self.dz[0])
        )
        self.varrying["dtmsl2"] = (
            wcfac * self.pamset["cpi"] / self.dz[self.pamset["lm"] - 1]
        )

    def setup_ebud(self):
        """
        Set up energy budget before run

        Find various paramaters at the start of the run,
        to get the energy budget ready to run
        """
        fnsa = 1.0  # Can it be something else?
        c1fac = self.pamset["dt"] / (self.pamset["c1"] * self.dz[0])

        blm = self.pamset["ebbeta"] / self.pamset["rlamdo"]
        self.gamn = (
            self.pamset["foan"] + self.pamset["rlamda"] / self.pamset["rlamdo"] + blm
        )
        self.gams = (
            self.pamset["foas"]
            + self.pamset["rlamda"] / self.pamset["rlamdo"]
            + fnsa * blm
        )

        # Northern hemisphere
        self.varrying["dtrm2n"] = (
            self.pamset["beto"]
            + self.pamset["foas"]
            * self.pamset["ebbeta"]
            / (self.gams * self.gamn - fnsa * blm**2)
        ) * c1fac
        self.varrying["dtrm3n"] = (
            self.gams / (self.gams * self.gamn - fnsa * blm**2) * c1fac
        )
        self.varrying["dtrm4n"] = blm / (self.gams * self.gamn - fnsa * blm**2) * c1fac

        # Southern hemisphere
        self.varrying["dtrm2s"] = (
            self.pamset["fnso"] * self.pamset["beto"]
            + self.pamset["foan"]
            * fnsa
            * self.pamset["ebbeta"]
            / (self.gams * self.gamn - fnsa * blm**2)
        ) * c1fac
        self.varrying["dtrm3s"] = (
            self.gamn / (self.gams * self.gamn - fnsa * blm**2) * c1fac
        )
        self.varrying["dtrm4s"] = (
            fnsa * blm / (self.gams * self.gamn - fnsa * blm**2) * c1fac
        )

        self.varrying["dtmnl3"] = (
            self.pamset["dt"]
            * self.pamset["beto"]
            / (self.pamset["c1"] * self.dz[self.pamset["lm"] - 1])
        )
        self.varrying["dtmnl1"] = 1.0 - self.varrying["dtmnl3"]
        self.varrying["dtmsl3"] = self.pamset["fnso"] * self.varrying["dtmnl3"]
        self.varrying["dtmsl1"] = 1.0 - self.varrying["dtmsl3"]
        self.setup_ebud2(0, 0)

    def energy_budget(
        self, forc_nh, forc_sh, fn_volc, fs_volc
    ):  # pylint: disable=too-many-locals, too-many-statements
        """
        Perform energy budget calculation for a single year.

        This is the main computational method of the Upwelling Diffusion Model.
        It calculates the thermal response of the ocean-atmosphere system to
        radiative forcing, including temperature changes, radiative imbalances,
        and ocean heat content changes for both hemispheres.

        The method implements the core physics of upwelling and diffusion in
        the ocean, accounting for hemisphere-specific forcing and volcanic
        perturbations. It returns comprehensive diagnostics including surface
        and air temperatures, ocean heat content, and radiative imbalances.

        Parameters
        ----------
        forc_nh : float
            Northern hemisphere radiative forcing (W/m²)
        forc_sh : float
            Southern hemisphere radiative forcing (W/m²)
        fn_volc : array
            Northern hemisphere volcanic forcing perturbation (W/m²)
        fs_volc : array
            Southern hemisphere volcanic forcing perturbation (W/m²)

        Returns
        -------
        dict
            Dictionary containing comprehensive thermal model outputs:

            Temperature Changes:
            - 'dtemp' : Global mean temperature change (K)
            - 'dtempnh' : Northern hemisphere temperature change (K)
            - 'dtempsh' : Southern hemisphere temperature change (K)
            - 'dtemp_air' : Global air temperature change (K)
            - 'dtempnh_air' : Northern hemisphere air temperature (K)
            - 'dtempsh_air' : Southern hemisphere air temperature (K)
            - 'dtemp_sea' : Global sea surface temperature change (K)
            - 'dtempnh_sea' : Northern hemisphere sea temperature (K)
            - 'dtempsh_sea' : Southern hemisphere sea temperature (K)

            Radiative Imbalances:
            - 'RIB' : Global radiative imbalance (W/m²)
            - 'RIBN' : Northern hemisphere radiative imbalance (W/m²)
            - 'RIBS' : Southern hemisphere radiative imbalance (W/m²)

            Ocean Heat Content (in 10²² J units):
            - 'OHCTOT' : Total ocean heat content change
            - 'OHC700' : Ocean heat content change down to 700m depth

        Notes
        -----
        The energy budget calculation involves:
        1. Processing hemisphere-specific and volcanic forcing
        2. Solving the ocean thermal diffusion equation using finite differences
        3. Calculating air and sea surface temperature responses
        4. Computing radiative imbalances and ocean heat content changes
        5. Combining hemisphere results into global averages

        The method maintains internal state (ocean temperatures by layer) that
        evolves between calls, making it suitable for time-series simulations.

        Examples
        --------
        >>> udm = UpwellingDiffusionModel({'lambda': 0.6})
        >>> # Apply 1 W/m² forcing uniformly
        >>> result = udm.energy_budget(1.0, 1.0, [0.0], [0.0])
        >>> print(f"Temperature change: {result['dtemp']:.3f} K")
        >>> print(f"Ocean heat content: {result['OHCTOT']:.3f} x10²² J")

        See Also
        --------
        _get_ocean_heat_content : Calculate ocean heat content diagnostics
        AbstractThermalModel.energy_budget : Parent interface method
        """
        # --- At the start of the year, store the initial temperature profiles ---
        tn_start = self.tn.copy()
        ts_start = self.ts.copy()
        temp1n = 0.0
        temp1s = 0.0

        tempn = 0.0
        temps = 0.0
        tempn_air = 0.0
        temps_air = 0.0
        tempn_sea = 0.0
        temps_sea = 0.0
        lm = self.pamset["lm"]
        templ = np.zeros(lm)

        dtyear = 1.0 / self.pamset["ldtime"]
        dn = np.zeros(lm)
        ds = np.zeros(lm)
        for im in range(self.pamset["ldtime"]):
            volc_idx = im % len(fn_volc)
            if self.pamset["threstemp"] != 0:  # pylint: disable=compare-to-zero
                self.setup_ebud2(temp1n, temp1s)

            dqn = (
                (im + 1) * forc_nh * dtyear
                + (1 - (im + 1) * dtyear) * self.prev_values["fn"]
                + fn_volc[volc_idx]
            )
            dqs = (
                (im + 1) * forc_sh * dtyear
                + (1 - (im + 1) * dtyear) * self.prev_values["fs"]
                + fs_volc[volc_idx]
            )
            dn[0] = (
                self.varrying["dtrm1n"] * self.tn[0]
                + self.varrying["dtrm2n"] * self.ts[0]
                + self.varrying["dtrm3n"] * dqn
                + self.varrying["dtrm4n"] * dqs
            )
            ds[0] = (
                self.varrying["dtrm1s"] * self.ts[0]
                + self.varrying["dtrm2s"] * self.tn[0]
                + self.varrying["dtrm3s"] * dqs
                + self.varrying["dtrm4s"] * dqn
            )
            dn[1 : lm - 1] = self.tn[1 : lm - 1] + self.pamset["beto"] * self.pamset[
                "dt"
            ] / (self.pamset["c1"] * self.dz[1 : lm - 1]) * (
                self.ts[1 : lm - 1] - self.tn[1 : lm - 1]
            )
            ds[1 : lm - 1] = self.ts[1 : lm - 1] + self.pamset["fnso"] * self.pamset[
                "beto"
            ] * self.pamset["dt"] / (self.pamset["c1"] * self.dz[1 : lm - 1]) * (
                self.tn[1 : lm - 1] - self.ts[1 : lm - 1]
            )

            dn[lm - 1] = (
                self.varrying["dtmnl1"] * self.tn[lm - 1]
                + self.varrying["dtmnl2"] * self.tn[0]
                + self.varrying["dtmnl3"] * self.ts[lm - 1]
            )
            ds[lm - 1] = (
                self.varrying["dtmsl1"] * self.ts[lm - 1]
                + self.varrying["dtmsl2"] * self.ts[0]
                + self.varrying["dtmsl3"] * self.tn[lm - 1]
            )

            # Where are these being initialised? Ok, I think
            self.tn = _band(
                self.varrying["acoeffn"],
                self.varrying["bcoeffn"],
                self.varrying["ccoeffn"],
                dn,
            )
            # Replace _band with _solve_tridiagonal
            self.ts = _band(
                self.varrying["acoeffs"],
                self.varrying["bcoeffs"],
                self.varrying["ccoeffs"],
                ds,
            )
            # print("self.aceoffn: %s self.varrying["bcoeffn"]: %s self.ccoeffn %s"%(self.varrying["acoeffn"], self.varrying["bcoeffn"], self.varrying["ccoeffn"]))
            # print("self.aceoffs: %s self.bcoeffs: %s self.ccoeffs %s"%(self.varrying["acoeffs"], self.bcoeffs, self.ccoeffs))
            temp1n = self.tn[0]
            temp1s = self.ts[0]
            # print("temp1n: %.5e temp1s %.5e"%(temp1n, temp1s))
            templ = (
                templ + 0.5 * (self.tn + self.ts) / self.pamset["ldtime"]  # 12.0
            )  # used to be 12, now ldtime

            tempan = (
                dqn
                + self.pamset["foan"] * self.pamset["rlamdo"] * temp1n
                + self.pamset["ebbeta"]
                * (dqs + self.pamset["foas"] * self.pamset["rlamdo"] * temp1s)
                / self.pamset["fsx"]
            )
            tempan = tempan / (
                self.pamset["fnx"] - self.pamset["ebbeta"] ** 2 / self.pamset["fsx"]
            )
            tempas = (
                dqs
                + self.pamset["foas"] * self.pamset["rlamdo"] * temp1s
                + self.pamset["ebbeta"]
                * (dqn + self.pamset["foan"] * self.pamset["rlamdo"] * temp1n)
                / self.pamset["fnx"]
            )
            tempas = tempas / (
                self.pamset["fsx"] - self.pamset["ebbeta"] ** 2 / self.pamset["fnx"]
            )
            tmpn = self.pamset["foan"] * temp1n + (1.0 - self.pamset["foan"]) * tempan
            tmps = self.pamset["foas"] * temp1s + (1.0 - self.pamset["foas"]) * tempas

            # x1=1638.+float(years_since_start)+float(im-1)/12.

            tempn = tempn + tmpn / self.pamset["ldtime"]  # Previously 12.0
            temps = temps + tmps / self.pamset["ldtime"]  # Previously 12.0

            tempn_air = tempn_air + tempan / self.pamset["ldtime"]  # Previously 12.0
            temps_air = temps_air + tempas / self.pamset["ldtime"]  # Previously12.0

            tempn_sea = tempn_sea + temp1n / self.pamset["ldtime"]  # Previously 12.0
            temps_sea = temps_sea + temp1s / self.pamset["ldtime"]  # Previously 12

        dtemp = (tempn + temps) / 2.0  # Global temp chg)

        # Updating previous values for next year
        self.prev_values["fn"] = forc_nh
        self.prev_values["fs"] = forc_sh
        self.prev_values["dtemp"] = dtemp

        # Getting Ocean temperature:
        ocean_res = self.ocean_temperature()

        # Getting anomalous radiation
        anomalous_radiation = self._calculate_anomalous_radiation(tn_start, ts_start)

        ribn = (
            forc_nh
            + np.mean(fn_volc)
            - self.pamset["rlamda"] * tempn
            - anomalous_radiation
        )
        ribs = (
            forc_sh
            + np.mean(fs_volc)
            - self.pamset["rlamda"] * temps
            - anomalous_radiation
        )

        # Returning results_dict
        return {
            "dtemp": dtemp,
            "dtempnh": tempn,
            "dtempsh": temps,
            "tn": self.tn,
            "ts": self.ts,
            "dtemp_air": (tempn_air + temps_air) / 2.0,
            "dtempnh_air": tempn_air,
            "dtempsh_air": temps_air,
            "dtemp_sea": (tempn_sea + temps_sea) / 2.0,
            "dtempnh_sea": tempn_sea,
            "dtempsh_sea": temps_sea,
            "RIBN": ribn,
            "RIBS": ribs,
            "RIB": (ribn + ribs) / 2.0,
            "OHC700": ocean_res["OHC700"],
            "OHCTOT": ocean_res["OHCTOT"],
            "anomalous_radiation": anomalous_radiation,
        }

    def _calculate_anomalous_radiation(self, tn_start, ts_start):
        """
        Calculate the anomalous radiation across year

        This should be called at the end of the year, using
        tn_start and ts_start from the beginning of the year

        Parameters
        ----------
        tn_start : np.ndarray
            Temperature in the Northern Hemisphere ocean layers
            start of the year
        ts_start : np.ndarray
            Temperature in the Southern Hemisphere ocean layers
            start of the year

        Returns
        -------
            float
        Anomalous radiation
        """
        # 1. Calculate the temperature change in each layer over the year.
        delta_tn = self.tn - tn_start
        delta_ts = self.ts - ts_start

        # 2. Calculate the change in heat content per unit of ocean area for each layer (in J/m^2).
        #    We use a standard value for the volumetric heat capacity of seawater.
        heat_capacity_volumetric = 4.184e6  # Joules per m^3 per Kelvin

        # Area-weight the temperature change, then multiply by heat capacity and layer thickness.
        delta_heat_content_profile = (
            heat_capacity_volumetric
            * self.dz
            * (delta_tn * self.pamset["foan"] + delta_ts * self.pamset["foas"])
        )

        # 3. Calculate the annual average heat uptake N (in W/m^2) by summing the
        #    heat content change of the DEEP OCEAN layers and dividing by the seconds in a year.
        year_in_seconds = DAYS_PER_YEAR * SECONDS_PER_DAY
        ann_heat_uptake = np.sum(delta_heat_content_profile[1:]) / year_in_seconds

        # 4. Calculate the TRUE Top-of-Atmosphere Radiative Imbalance (RIB).
        ocean_efficacy = self.pamset.get("ocean_efficacy", 1.0)

        # The anomalous radiation to space is (E-1)*ann_heat_uptake.
        anomalous_radiation = (ocean_efficacy - 1.0) * ann_heat_uptake

        return anomalous_radiation

    def ocean_temperature(self):
        """
        Compute the ocean temperature total and at 700 m depth

        Compute the ocean temperature total and at 700 m depth
        to get ocean heat content for total and down to 700 m

        Returns
        -------
        dict
            Containing the ocean heat content to 700 m and total
            with keys OHC700 and OHCTOT
        """
        havtemp = (
            RHO_SEAWATER
            * C_OHC
            * OCEAN_AREA_HEMISPHERE
            * self.dz
            * (self.tn * self.pamset["foan"] + self.ts * self.pamset["foas"])
        )

        # Finding the max layer down to 700m
        max_layer = int(INDEX_700M - self.dz[0] // DZ)
        frac = (1 + self.dz[0] // DZ) - self.dz[0] / DZ

        return {
            "OHC700": np.sum(havtemp[:max_layer]) + frac * havtemp[max_layer],
            "OHCTOT": np.sum(havtemp[:]),
        }
