# NEXD 2D

NEXD is a software package for high order simulation of seismic waves using the nodal discontinuous Galerkin method.
NEXD 2D is the 2D solver of this software package.

To deal with highly complex heterogeneous models, here the Nodal Discontinuous Galerkin Method (NDG) is used to calculate synthetic seismograms. The advantage of this method is that complex mesh geometries can be computed by using triangular or tetrahedral elements for domain discretization together with a high order spatial approximation of the wave field. The simulation tool NEXD is developed which has the capability of simulating elastic and anelastic wave fields for seismic experiments for one-, two- and three- dimensional settings. Porous materials and fractures (in terms of linear slop interfaces) are working for the one- and two-dimensional versions. External models provided by e.g. Trelis/Cubit can be used for parallelized computations on triangular or tetrahedral meshes. For absorbing boundary conditions either a fluxes based approach or a Nearly Perfectly Matched Layer (NPML) can be used. In addition full waveform inversion is supported in 2D.

## Authors and License

NEXD 2D and some of its components, as well as documentation and some examples
are available under terms of the [GNU General Public License](LICENSE) (version 3 or higher)
on [github](https://github.com/seismology-RUB/NEXD-2D).
Please find contact addresses [there](https://github.com/seismology-RUB), or visit 
http://www.rub.de/nexd in case you want to get in touch with the authors. If you 
encounter any problems installing or using the software, it will be helpful to 
open (or add to) an "issues" topic at the [github repository](https://github.com/seismology-RUB/NEXD-2D).

The main authors are Lasse Lambrecht, Andre Lamert, Wolfgang Friederich, Thomas Möller (all Ruhr-Universität Bochum, Germany) and Marc S. Boxberg (now at RWTH Aachen University, Germany).


## Documentation

A user manual can be found in the folder 'doc'.


## Examples

NEXD 2D comes with five examples: Lamb's problem (simulations/lambs), a two layer problem (simulations/2layer_pml_att), an example for including fractures (simulations/example_linearSlipInterface), an example for poroelastic wave propagation (simulations/example_poro), and an example for full waveform inversion (simulations/example_inversion). Please use these examples to familiarize yourself with the software.


## Requirements

* GNU Make
* Fortran compiler (sufficient standard)
* LAPACK libraries
* MPI libraries (e.g., OpenMPI) and METIS (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview; tested with version 4.0.3) for parallel applications


## Installation

0. If not yet done, you should download the source code of the NEXD 2D main package by cloning the master branch of the NEXD 2D repository on gitHub.com:
     ```
     git clone --depth 1 --branch master https://github.com/seismology-RUB/NEXD-2D
     ```

1. Install all software dependencies.

2. Adjust the software to your system and personal requirements by changing the [Makefile](Makefile) appropriately (e.g., change the path to your METIS installation and set your compiler).

3. Run command
     ```
     make all
     ```
   from installation path `NEXD-2D/` to compile NEXD 2D.
