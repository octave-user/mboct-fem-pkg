# mboct-fem-pkg<sup>&copy;</sup>
**mboct-fem-pkg** belongs to a suite of packages which can be used for pre- and post-processing of flexible bodies in MBDyn (https://www.mbdyn.org) with GNU-Octave (http://www.gnu.org/software/octave/) and Gmsh (http://www.gmsh.info/).
It contains a general purpose Finite Element toolkit for linear structural, linear thermal, linear acoustic and fluid structure interaction problems. It is supporting also nonlinear multibody dynamics via MBDyn and can be used to generate flexible body data for MBDyn's modal element.
In addition to that, there are functions for pre- and post-processing of MBDyn's elastohydrodynamic journal- and slider- plain bearings.

# List of features
  - Import and export meshes from and to Gmsh
  - Assign material properties:
           Structural: isotropic and arbitrary an-isotropic linear elastic constitutive laws
                       Rayleigh damping
                       density
           Thermal: specific heat capacity
                    thermal conductivity
           Acoustic: speed of propagation of sound
                     fluid density
                     fluid viscosity
                     
  - Assign loads and boundary conditions:
           Structural: nodal forces and moments
                       distributed forces and moments via RBE3
                       pressure loads
                       thermal expansion and arbitrary pre-strains
                       prescribed nodal displacements
                       arbitrary linear constraints for nodal displacements
                       bonded and sliding connections between meshes
                       gravity loads
                       prescribed angular velocity of a moving reference frame
                       prescribed angular acceleration of a moving reference frame
           Thermal: imposed nodal temperatures
                    imposed heat source
                    imposed convection
           Acoustic: imposed acoustic pressure
                     imposed particle velocity
                     imposed acoustic impedance
                     radiation boundary conditions
                     fluid-structure interfaces

  - Assemble finite element matrices and right hand side vectors
  - Solve linear static- and dynamic structural problems also for rotating and preloaded structures
  - Solve undamped and damped structural eigenmode problems
  - Solve linear buckling problems
  - Compute stresses and strains
  - Create reduced order models (multibody dynamics analysis in MBDyn assuming small elastic deformations but arbitrary rigid body motions)
  - Expand modal results from MBDyn from several flexible bodies, scale elastic deformations, compute stresses and perform post-processing of combined meshes from several flexible bodies via Gmsh
  - Generate compliance matrices for elastohydrodynamic journal and slider bearings
  - Solve linear stationary and transient thermal problems
  - Solve linear harmonic acoustic and fluid structure interaction problems
  - Solve acoustic eigenmode problems
  - Compute particle velocities, acoustic intensity and acoustic radiation


# Supported element types and application
  - 8 node linear hexahedral 3D element
  - 6 node linear pentahedral 3D element
  - 4 node linear tetrahedral 3D element
  - 20 node quadratic hexahedral 3D element
  - 15 node quadratic pentahedral 3D element
  - 10 node quadratic tetrahedral 3D p-element
  - 10 node quadratic tetrahedral 3D h-element
  - 20 node cubic tetrahedral 3D h-element
  - 4 node linear quadrilateral 3D surface element
  - 3 node linear triangular 3D surface element
  - 8 node quadratic quadrilateral 3D surface element
  - 6 node quadratic triangular 3D surface p-element
  - 6 node quadratic triangular 3D surface h-element
  - 2 node Timoshenko beam element (structural)
  - Load distributing Nastran like RBE3 element
  - Generic joint element
  - Surface to node constraint elements for coupling of unrelated meshes
  - Lumped force/torque element

Copyright<sup>&copy;</sup> 2019-2022

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