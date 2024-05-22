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
- Moved carbon cycle functionality out of concentrations_emissions_handler.py to separate module carbon_cycle_mod.py

### Added 
- Functionality to back calculate CO2 emissions from concentrations
- Functionality to calculate carbon pool fluxes to ocean and land
- Functionality to calculate airborne fraction
- Functionality to output new carbon cycle information 

### Fixed

- Removed `jupyter nbextension` from `Makefile` causing install errors
- Locked numpy in requirements to below version 2 to work with scmdata for now

## [Versions 1.1.1 and v1.1.0]


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