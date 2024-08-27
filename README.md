# mboct-fem-pkg<sup>&copy;</sup>
**mboct-fem-pkg** belongs to a suite of packages which can be used for pre- and post-processing of flexible bodies in MBDyn (https://www.mbdyn.org) with GNU-Octave (http://www.gnu.org/software/octave/) and Gmsh (http://www.gmsh.info/).
It contains a general purpose Finite Element toolkit for linear structural, linear thermal, linear acoustic and fluid structure interaction problems. It is supporting also nonlinear multibody dynamics via MBDyn and can be used to generate flexible body data for MBDyn's builtin modal element. Modal elements are based on the assumption of linear elastic material, small elastic deformations but arbitrary superimposed rigid body motions. In case of nonlinear material and/or large deformations, mboct-fem-pkg can be used to create input data and to load output data for MBDyn's builtin 3D solid elements.
In addition to that, mboct-fem-pkg has also functions for pre- and post-processing of MBDyn's builtin elastohydrodynamic journal- and slider- plain bearings.

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
  - Perform pre- and postprocessing of 3D flexible bodies subject to nonlinear material and large deformations
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

Copyright<sup>&copy;</sup> 2019-2024

[Reinhard](mailto:octave-user@a1.net)

# Installation
  See the workflow file [ubuntu-latest.yml](https://github.com/octave-user/mboct-fem-pkg/blob/master/.github/workflows/ubuntu-latest.yml) as an example on how to install the package.

# Docker images
  Docker images may be pulled by one of the following commands:
  - docker pull ghcr.io/octave-user/mboct-fem-pkg:master
  - docker pull octaveuser/mboct-fem-pkg:latest

# Run mboct-fem-pkg within docker
  - On any operating system with support for docker, the container can be executed by using the following command:
  docker run -it ghcr.io/octave-user/mboct-fem-pkg:master
  - If the user wishes to start the container in graphics mode and the host system is based on Ubuntu, then the following commands should be used:

    xhost + local:docker

    docker run -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY -h $HOSTNAME -v $HOME/.Xauthority:/home/ubuntu/.Xauthority  -it ghcr.io/octave-user/mboct-fem-pkg:master --gui

# Applications
  - Deformation of an elasto-hydrodynamic lubricated diaphragm plain bearing by Norman Owen Freund 1995 ([#1](https://www.youtube.com/watch?v=YE0gnTt35WA), [#2](https://www.youtube.com/watch?v=tipxGDXe1mI), [#3](https://www.youtube.com/watch?v=akkelq04mrU))

  - Test case for double elasto-hydrodynamic journal plain bearing ([#1](https://youtu.be/kCVneVwXYbc), [#2](https://youtu.be/eienVfAFyfk))

  - Rotor dynamics test case using MBDyn ([#1](https://www.youtube.com/watch?v=VohVTeggqI4))

  - Large deflection analysis of a cantilever beam using MBDyn and Gmsh ([#1](https://youtu.be/j8D821HVXDc))

  - Nonlinear elasticity of a twisted bar using MBDyn and Gmsh ([#1](https://youtu.be/D2OZHT9luQs))

  - Cook's membrane benchmark with large deformations and hyperelastic material ([#1](https://youtu.be/EAgejp4jQ00))

  - Noise radiation of a disc surrounded by air ([#1](https://youtu.be/I8R0HouG2Ck), [#2](https://youtu.be/b4oc1F3Wc3I))

  - Propagation of acoustic waves in a flexible pipe ([#1](https://youtu.be/N7NdN70kHRQ))