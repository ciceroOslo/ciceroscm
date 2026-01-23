Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_, and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

The changes listed in this file are categorised as follows:

    - Added: new features
    - Changed: changes in existing functionality
    - Deprecated: soon-to-be removed features
    - Removed: now removed features
    - Fixed: any bug fixes
    - Security: in case of vulnerabilities.

[Unreleased]
---------------------------
### Added
- Deployment workflow to publish package to PyPI on manual trigger

[Version 2.0.1]
---------------------------
### Added

### Changed

### Fixed

- Ensure code can run when ozone is not changing, avoiding divide by zero.
- Ensure calculation of airborne fraction gives nans when cumulative emissions are zero.
- Fix bug affecting backward compatibility with not sending explicit carbon_cycle argument to concentrations_emissions_handler

### Removed

[Version 2.0.0]
---------------------------

### Changed
- Simplified Tropospheric Ozone forcing calculation, using the same calculation throughout both before and after emstart 
    (previous implementation relied on fossil fuel CO2 concentration likely for historical reasons of missing data).

### Added

- Option to chunk configuration distributions to allow for millions of member runs without killing memory
- Option to output ozone and stratospheric water vapour effective radiative forcing in parallel run
- _config_distro can now write distribution files with full setup metadata, and Distroconfig can read in files with or with out shared metadata.
-  New structural switch functionality to allow for different structural configurations to be run for both carbon cycle and thermal model.
- A simplified box carbon model and a simplified two-layer ocean model as structural options.
- Option to run with regionally split aerosol species, also including functionality to make gaspamfile that supports this from forcing txt file.

### Fixed

- Skip precaclulation of empty concentrations matrix if emstart >= nyend
- Carbon cycle output functioninng for concentrations runs
- Parallel concentrations run pass back-calculated emissions series when prompted for Emissions|CO2
- Fix so standalone concentrations_emissions_handler will run without explicitly sending carbon cycle model argument, ensuring backward compatibility.

### Removed
- Removed old calibrator code, as calibration is now done in a separate repository.
- fnso (fraction of Northern hemisphere ocean to Southern hemisphere ocean) parameter removed from upwelling diffusion model. This parameter should be derived from the ratio between foan and foas, and not be something that can be set separately.

[Version 1.5.0]
---------------------------

### Changed

- Change so that non-changing gases are precalculated before a run in concentrations_emissions_handler.py to speed up and avoid looping over every component and looking up in pandas every timestep
- Removed use of scmdata, and changed parallel output to be pandas to offload extra dependencies
- Moved carbon cycle functionality out of concentrations_emissions_handler.py to separate module carbon_cycle_mod.py
- Changed methodology for parametrised decay functions for biotic decay (rb_function) and mixed layer pulse response rs_function, imposing a few more consistency criteria on them and makeing them changable by a dictionary for each with keys "coeffs" and "timescales" with values that should be lists or np.arrays of coefficient and timescale values 
- Upgrade infrasturucture for development
- Output data / output files no longer include empty columns for species that don't make sense (landuse albedo emissions, concentrations for aerosols like BC etc)
- Aerosol forcing is no longer scaled to a reference year, but is per change in annual emissions, this also changes aerosol forcing parameter values and units
- Added NMVOC, NH3 and NOx as forcing aerosols
- Changed default gaspam unit for SO2 to Tg_SO2 rather than Tg_S
- Carbon cycle parameters now passed in a separate dictionary
- Simplifying _config_distro


### Added

- Functionality to feedback current temperature change to the carbon cycle and have simple temperature feedback 
- Functionality to back calculate CO2 emissions from concentrations
- Functionality to calculate carbon fluxes to ocean and land
- Functionality to calculate airborne fraction
- Functionality to output new carbon cycle information
- Functionality to output new carbon cycle information in parallel runs
- Changeable parameter for the depth of the mixed_layer ocean seen by the carbon cycle
- Newer gaspam files
- Light check for emissions and gaspam file compatibility
- Changeable ocean efficacy parameter to the upwelling diffusion model
- Carbon cycle temperature feedbacks for mixed_layer ocean depth (sigmoid decline), 
- Land primary production (optimum increase and threshold dampening decline) and
- Ocean solubility (exponential scaling to limit)


### Fixed

- Removed `jupyter nbextension` from `Makefile` causing install errors
- Locked numpy in requirements to below version 2 to work with scmdata for now
- Fix error in flat N2O assumption in natural precalculate_natural_emissions.py script
- Change logging messages for misisng parameters from warning to info
- Fix bug in file cutting for read_inputfile in input_handler
- Temporarily mute docs generation with sphinx on github as it is not working, plan to work in different method
- Temporarily mute distribution generation, plan to work this back in with more updated package manager

[Versions 1.1.1 and v1.1.0]
---------------------------

### Added

- Support for flat natural emissions for NO2 and CH4 from gaspamfile as defaults
- Update script to make natural emissions for NO2 and CH4 using ODE solution
- Updates to default vulcanoe and solar forcing timeseries
- Rationalising _band implementation using library for speed and readibility
- Clearing out unused density realated methods
- Made the fertilisation factor beta_f part a changeable parameter
- Parallelisation wrapper for parallel runs
- Config distributions support to make parameter distributions
- Distrobution runner to run over distributions in parallel
- Calibrator to make calibrated set of configurations based on calibration data
- Made automatic plots more polished looking
- Taking out unused parts from gaspam-file
