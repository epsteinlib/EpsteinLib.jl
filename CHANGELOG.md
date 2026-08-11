# Changelog

## [0.6.1] - unreleased

### Added
- Julia wrappers: `epsteinzetaaniso`, `epsteinzetaanisoreg`, supporting optional arguments.
- Unit tests for the anisotropic variants.
- Enhanced guards against non-matching dimensions in input values.

## [0.5.1] - 2026-04-27

### Changed
- Updated to Epsteinlib_jll `v0.5.1`, see [EpsteinLib changelog](https://github.com/epsteinlib/epsteinlib/blob/main/CHANGELOG.md) for details.

### Added
- Standalone minimal example in `README.md` uses the registered wrapper, e.g. `Pkg.add("EpsteinLib")`.

## [0.5.0] - 2025-12-04

_First public release_

### Added
- Julia wrappers: `epsteinzeta`, `epsteinzetareg`, supporting optional arguments.
- Unit test suite for core functionality.
