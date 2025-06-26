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

## [Unreleased]

### Changed
- Removed use of scmdata, and changed parallel output to be pandas to offload extra dependencies
- Moved carbon cycle functionality out of concentrations_emissions_handler.py to separate module carbon_cycle_mod.py
- Changed methodology for parametrised decay functions for biotic decay (rb_function) and mixed layer pulse response rs_function, imposing a few more consistency criteria on them and makeing them changable by a dictionary for each with keys "coeffs" and "timescales" with values that should be lists or np.arrays of coefficient and timescale values 
- Upgrade infrasturucture for development

### Added
- Functionality to feedback current temperature change to the carbon cycle and have simple temperature feedback 
- Functionality to back calculate CO2 emissions from concentrations
- Functionality to calculate carbon fluxes to ocean and land
- Functionality to calculate airborne fraction
- Functionality to output new carbon cycle information
- Functionality to output new carbon cycle information in parallel runs
- Changeable parameter for the depth of the mixed_layer ocean seen by the carbon cycle


### Fixed

- Removed `jupyter nbextension` from `Makefile` causing install errors
- Locked numpy in requirements to below version 2 to work with scmdata for now
- Fix error in flat N2O assumption in natural precalculate_natural_emissions.py script
- Change logging messages for misisng parameters from warning to info
- Fix bug in file cutting for read_inputfile in input_handler
- Temporarily mute docs generation with sphinx on github as it is not working, plan to work in different method
- Temporarily mute distribution generation, plan to work this back in with more updated package manager

## [Versions 1.1.1 and v1.1.0]
  
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