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

### Added

- Support for flat natural emissions for NO2 and CH4 from gaspamfile as defaults
- Updates to default vulcanoe and solar forcing timeseries
- Rationalising _band implementation using library for speed and readibility
- Clearing out unused density realated methods
- Made the fertilisation factor beta_f part a changeable parameter
- Parallelisation wrapper for parallel runs
- Config distributions support to make parameter distributions
- Distrobution runner to run over distributions in parallel
- Calibrator to make calibrated set of configurations based on calibration data
