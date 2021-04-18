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

# Supported element types and application
  - 8 node linear hexahedral 3D element (structural, thermal)
  - 6 node linear pentahedral 3D element (structural, thermal)
  - 4 node linear tetrahedral 3D element (structural, thermal)
  - 20 node quadratic hexahedral 3D element (structural, thermal)
  - 15 node quadratic pentahedral 3D element (structural, thermal)
  - 10 node quadratic tetrahedral 3D p-element (structural, thermal)
  - 10 node quadratic tetrahedral 3D h-element (structural, thermal)
  - 4 node linear quadrilateral 3D surface element (surface to node constraint, pressure load, thermal convection, heat source)
  - 3 node linear triangular 3D surface element (surface to node constraint, pressure load, thermal convection, heat source)
  - 8 node quadratic quadrilateral 3D surface element (surface to node constraint, pressure load, thermal convection, heat source)
  - 6 node quadratic triangular 3D surface p-element (surface to node constraint, pressure load, thermal convection, heat source)
  - 6 node quadratic triangular 3D surface h-element (surface to node constraint, pressure load, thermal convection, heat source)
  - 2 node Timoshenko beam element (structural)
  - Load distributing Nastran like RBE3 element (structural)
  - Generic joint element (structural, thermal)
  - Surface to node constraint elements for coupling between unrelated meshes (structural, thermal)
  - Lumped force/torque element (structural)

Copyright<sup>&copy;</sup> 2019-2021

[Reinhard](mailto:octave-user@a1.net)

# Installation
  ## Ubuntu 20.04
  The following code is an example how mboct-fem-pkg can be installed on an Ubuntu system:

  `sudo apt-get install octave liboctave-dev libsuitesparse-dev libarpack2-dev libmumps-seq-dev libmetis-dev octave-nurbs gmsh libnlopt-dev libmkl-full-dev coreutils`

  `git clone -b develop https://public.gitlab.polimi.it/DAER/mbdyn.git`

  `pushd mbdyn`

  `./bootstrap.sh`

  `./configure --with-static-modules --enable-octave --enable-sparse_autodiff --enable-autodiff --disable-Werror CXXFLAGS="-O3 -march=native" CPPFLAGS=-I/usr/include/suitesparse --with-arpack --with-umfpack --without-metis`

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

  ## openSUSE 15.2
  The following code is an example how mboct-fem-pkg can be installed on an openSUSE system.

  `sudo zypper install octave octave-devel octave-forge-nurbs nlopt-devel suitesparse-devel autoconf automake libtool git arpack-ng-devel ginac-devel`

  `wget http://gmsh.info/bin/Linux/gmsh-4.8.3-Linux64.tgz`

  `tar -xvf gmsh-4.8.3-Linux64.tgz`

  `sudo install gmsh-4.8.3-Linux64/bin/gmsh /usr/local/bin`

  `git clone -b develop https://public.gitlab.polimi.it/DAER/mbdyn.git`

  `pushd mbdyn`

  `./bootstrap.sh`

  `./configure --enable-octave --with-umfpack --with-arpack --without-metis --with-static-modules --enable-sparse_autodiff --with-lapack --with-suitesparseqr CPPFLAGS=-I/usr/include/suitesparse CXXFLAGS=-Wall -ggdb3 -march=native -O3 --disable-Werror`

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
  
# Applications
  - Deformation of an elasto-hydrodynamic lubricated diaphragm plain bearing by Norman Owen Freund 1995 (https://www.youtube.com/watch?v=YE0gnTt35WA, https://www.youtube.com/watch?v=tipxGDXe1mI, https://www.youtube.com/watch?v=akkelq04mrU)
  
  - Rotor dynamics test case using MBDyn (https://www.youtube.com/watch?v=VohVTeggqI4)	
