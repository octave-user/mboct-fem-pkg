# mboct-fem-pkg<sup>&copy;</sup>
**mboct-fem-pkg** belongs to a suite of packages which can be used for pre- and post-processing of flexible bodies in MBDyn (https://www.mbdyn.org) with GNU-Octave (http://www.gnu.org/software/octave/) and Gmsh (http://www.gmsh.info/).
It contains a general purpose Finite Element toolkit for linear structural and linear thermal problems, which can be used to generate flexible body data for MBDyn's modal element. 
In addition to that, there are functions for pre- and post-processing of MBDyn's elastohydrodynamic journal- and slider- plain bearings.

# List of features
  - Generate meshes (builtin structured hexahedral mesh generator for simple meshes, general unstructured meshes with Gmsh)
  - Assign material properties (arbitrary isotropic or anisotropic linear elastic constitutive laws, Rayleigh damping, density, specific heat capacity, thermal conductivity, thermal expansion)
  - Assign loads (lumped nodal forces and torques, distributed loads via RBE3 elements, arbitrary pressure loads on surfaces with variable pressure, arbitrary prescribed displacements and rotations, acceleration loads, thermal strain, thermal convection, heat source, prescribed temperatures)
  - Assign boundary conditions (lock selected nodal displacements and rotations in the global reference frame)
  - Impose constraints (arbitrary linear relationships between a set of nodes, connect two surfaces with compatible or incompatible meshes for structural or thermal problems)
  - Assemble finite element matrices and right hand side vectors (mass matrix, stiffness matrix, damping matrix, lumped or consistent structural load vector, heat capacity matrix, thermal conductivity matrix, thermal load vector and invariants of inertia needed for multibody dynamics analysis)
  - Solve linear static structural problems
  - Solve structural eigenvalue problems
  - Compute deformations and stresses (stress tensor, von Mises stress)
  - Solve linear stationary thermal problems
  - Solve linear transient thermal problems
  - Import and export meshes and results (Gmsh, EOSSP, APDL)
  - Create reduced order models (multibody dynamics analysis in MBDyn assuming small elastic deformations but arbitrary rigid body motions)
  - Expand modal results from MBDyn from several flexible bodies, scale elastic deformations, compute stresses and perform post-processing of combined meshes from several flexible bodies in Gmsh
  - Generate compliance matrices for elastohydrodynamic journal and slider bearings

Copyright<sup>&copy;</sup> 2019-2021

[Reinhard](mailto:octave-user@a1.net)

# Installation
 The following code is an example how mboct-fem-pkg can be installed on an Ubuntu system:
 
   `sudo apt-get install octave liboctave-dev libsuitesparse-dev libarpack2-dev libmumps-seq-dev libmetis-dev octave-nurbs gmsh libnlopt-dev libmkl-full-dev coreutils`

  `git clone -b develop https://public.gitlab.polimi.it/DAER/mbdyn.git`

  `pushd mbdyn`

  `./bootstrap.sh`

  `./configure --with-static-modules --enable-octave --enable-sparse_autodiff --enable-autodiff --disable-Werror CPPFLAGS=-I/usr/include/suitesparse --with-arpack --with-umfpack --without-metis`

  `make`

  `sudo make install`

  `popd`

  `git clone -b master https://github.com/octave-user/mboct-octave-pkg.git`

  `make -C mboct-octave-pkg install_local`

  `git clone -b master https://github.com/octave-user/mboct-numerical-pkg.git`

  `make -C mboct-numerical-pkg install_local`

  `git clone -b master https://github.com/octave-user/mboct-mbdyn-pkg.git`

  `make -C mboct-mbdyn-pkg install_local`

  `git clone -b master https://github.com/octave-user/mboct-fem-pkg.git`

  `make -C mboct-fem-pkg install_local`
    
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

## NLOpt installation:
  - Follow the instructions on https://nlopt.readthedocs.io to install nlopt.

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
	
