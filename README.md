# mboct-fem-pkg<sup>&copy;</sup>
**mboct-fem-pkg** belongs to a suite of packages which can be used for pre- and post-processing of flexible bodies in MBDyn (https://www.mbdyn.org) with GNU-Octave (http://www.gnu.org/software/octave/) and Gmsh (http://www.gmsh.info/).
It contains a general purpose structural Finite Element toolkit for linear statics and dynamics, which can be used to generate flexible body data for MBDyn's modal element. 
In addition to that, there are functions for pre- and post-processing of MBDyn's elastohydrodynamic journal- and slider- plain bearings.

# List of features
  - Generate meshes (builtin structured hexahedral mesh generator for simple meshes, general unstructured meshes with Gmsh)
  - Assign material properties (arbitrary isotropic or anisotropic linear elastic constitutive laws, Rayleigh damping, density)
  - Assign loads (lumped nodal forces and torques, distributed loads via RBE3 elements, arbitrary pressure loads on surfaces with variable pressure, arbitrary prescribed displacements and rotations, acceleration loads)
  - Assign boundary conditions (lock selected nodal displacements and rotations in the global reference frame)
  - Impose constraints (arbitrary linear relationships between a set of nodes, connect two flexible surfaces with compatible or incompatible meshes with or without sliding between surfaces)
  - Assemble finite element matrices and right hand side vectors (mass matrix, stiffness matrix, damping matrix, lumped or consistent load vector, invariants of inertia needed for multibody dynamics analysis)
  - Solve linear static problems
  - Solve eigenvalue problems
  - Compute deformations and stresses (stress tensor, von Mises stress)
  - Import and export meshes and results (Gmsh, EOSSP, APDL)
  - Create reduced order models (multibody dynamics analysis in MBDyn assuming small elastic deformations but arbitrary rigid body motions)
  - Expand modal results from MBDyn from several flexible bodies, scale elastic deformations, compute stresses and perform post-processing of combined meshes from several flexible bodies in Gmsh
  - Generate compliance matrices for elastohydrodynamic journal and slider bearings

Copyright<sup>&copy;</sup> 2019-2020

[Reinhard](mailto:octave-user@a1.net)

# Installation

## GNU Octave installation
  - Follow the instructions on (http://www.gnu.org/software/octave/) to install GNU Octave.  
  - Make sure, that `mkoctfile` is installed.  
    `mkoctfile --version` 

## MBDyn installation:
  - Clone the source tree of MBDyn.  
    `git clone https://public.gitlab.polimi.it/DAER/mbdyn.git -b develop`
  - Compile and install MBDyn.  
    `cd mbdyn`  
    `./bootstrap.sh`  
    `./configure CXXFLAGS=-O3 --enable-octave --enable-autodiff --with-static-modules --with-umfpack`  
    `make`  
    `make install`

## GNU Octave package installation:
  - Make sure that the GNU Octave nurbs package is installed.  
    `octave --eval 'pkg install -forge nurbs'`
  - Install the following packages from github.  
    `for pkg in octave numerical mbdyn fem; do`    
        `git clone https://github.com/octave-user/mboct-${pkg}-pkg.git && make -C mboct-${pkg}-pkg install_local`	  
    `done`

## Gmsh installation:
  - Follow the instructions on (http://www.gmsh.info/) to install Gmsh.  
  - Make sure Gmsh is on your path (e.g. `export PATH=/opt/gmsh/bin:${PATH}`)

## Usage
  - Run Octave.  
    `octave`
  - At the Octave prompt load the package.   
    `pkg load mboct-fem-pkg`
  - At the Octave prompt execute a demo.  
    `demo fem_cms_export`
	