# Changelog
All notable changes to NEXD 2D will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.4] - 2020-04-22
### Added
- Full waveform inversion
- Interface for meshs created by gmsh

### Changed
- Definition of fractures / linear slip interfaces now requires the relaxation frequency instead of the thickness and elastic moduli.
- Definitions of reference frequencies for viscoelastic attenuation are changed to frequencies instead of angular frequencies.

### Removed
- Flag 'debug' in parfile


## [0.3] - 2019-12-13
### Added
- Porous material saturated by one or two fluids
- Subsampling for seismograms (sampling rate of seismograms can now be an integer multiple of the simulation timestep)
- Rotation angles can now be set for each receiver individually.
- Duration of simulation can now be set in terms of total simulated time in addition to the number of timesteps.
- Sources can have a source time other than 0.


## [0.2.1] - 2019-03-20
### Fixed
- In the case of several sources in the same partition of the mesh, only one source was active.


## [0.2] - 2018-09-20
### Added
- Fractures in terms of linear slip interfaces can be included in the model.
- Manual

### Removed
- Weak form of PDE


## [0.1] - 2017-03-23
### Added
- Wave propagation in elastic and viscoelastic materials
