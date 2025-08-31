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

Copyright<sup>&copy;</sup> 2019-2025

[Reinhard](mailto:octave-user@a1.net)

# Function reference
  - The function reference is automatically generated from the source code by means of Octave's [generate_html](https://octave.sourceforge.io/generate_html/index.html) package. See [overview.html](https://octave-user.github.io/mboct-fem-pkg/mboct-fem-pkg/overview.html).

# Installation
  See the workflow file [simple.yml](https://github.com/octave-user/mboct-fem-pkg/blob/master/.github/workflows/simple.yml) as an example on how to install the package.

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

# Youtube channel
  You can find several videos at the [youtube channel](https://www.youtube.com/@nonlinearmultibodydynamics6802) of mboct-fem-pkg.

# Applications
  - modal analysis of a spherical shell [deformation](https://youtu.be/SLhWWGIIvXU)

  - modal analysis of a flexible ring [deformation](https://youtu.be/wcXl09-jTU8)

  - modal analysis of a plate [deformation](https://youtu.be/e2J_W5Y31VU)

  - modal analysis of an elbow structure [deformation](https://youtu.be/Un6r08ueieM)

# Elasto hydrodynamic lubricated bearings (MBDyn)
  - Deformation of an elasto-hydrodynamic lubricated diaphragm plain bearing by Norman Owen Freund 1995 [deformation](https://www.youtube.com/watch?v=YE0gnTt35WA), [pressure](https://youtu.be/akkelq04mrU), [fractional film content](https://www.youtube.com/watch?v=tipxGDXe1mI)

  - Double elastohydrodynamic contact [deformation](https://youtu.be/kCVneVwXYbc), [pressure](https://youtu.be/eienVfAFyfk)

# Modal elements (MBDyn)
  - Rotor dynamics test case using MBDyn [rotor dynamics](https://www.youtube.com/watch?v=VohVTeggqI4)

  - Flexible pendulum [deformation](https://youtu.be/ggihEM9Dtn4)

  - eOSSP BMW i4 crankshaft (deformation and Von Mises stress) - [deformation rotating reference frame](https://youtu.be/DGmlRYLxtuA), [fluid pressure](https://youtu.be/JzW3mUuxoHE), [deformation inertial reference frame](https://youtu.be/ZV7OLzjZxpE)

  - Flexible body dynamics of a double pendulum [deformation](https://youtu.be/I4FE0zunTDo)

# Solid elements (MBDyn)
  - Large deflection analysis of a cantilever beam using MBDyn and Gmsh [deformation](https://youtu.be/j8D821HVXDc)

  - Nonlinear elasticity of a twisted bar using MBDyn and Gmsh [deformation](https://youtu.be/D2OZHT9luQs)

  - Cook's membrane benchmark with large deformations and hyperelastic material [deformation](https://youtu.be/EAgejp4jQ00)

  - Rolling flexible ring made of hyperelastic material [deformation](https://youtu.be/rxQP8V4U0dE)

  - Nonlinear hyperviscoelasticity using MFront's Signorini model [deformation](https://youtu.be/I8HENx5mszA)

  - Flexible four-bar linkage [deformation](https://youtu.be/d4i5AYPxsG4)

  - Torsion of a beam meshed with 20-node tetrahedral elements [deformation](https://youtu.be/OGmW8F54Vc0)

  - Model of a tire with internal pressure and contact between tire and ground [deformation](https://youtu.be/KTvsuRnZuGo)

  - 1/4 pipe under radial load [deformation](https://youtu.be/0pGlOycfWRw)

  - Bending of a cantilever beam based on K.J. Bathe [deformation](https://youtu.be/RXjAefPeG6Y)

# Heat equation
  - 1D heat equation based on Matthew J. Hancock [temperature](https://youtu.be/TJgq5sdHM8Q) [temperature profile](https://youtu.be/qgsaxq3t6EQ)

# Helmholtz equation
  - Noise radiation of a disc surrounded by air [deformation](https://youtu.be/I8R0HouG2Ck), [sound pressure](https://youtu.be/b4oc1F3Wc3I)

  - Propagation of acoustic waves in a flexible pipe [deformation](https://youtu.be/N7NdN70kHRQ)

  - 2D wave equation based on Jont Allen [sound pressure](https://youtu.be/vdebCp1RNoE)

  - Reflection and transmission of sound at an interface [pressure](https://youtu.be/1cdsKLKrI3Y)

  - Spherical wave equation with perfectly matched layers based on Jont Allen and Radu Cimpeanu [sound pressure](https://youtu.be/qhwL9JmkvEA)

  - Spherical wave equation with perfectly matched layers based on Jonathan Deakin [pressure](https://youtu.be/ypimZzAQUF8)

  - 1D dissipative wave equation with perfectly matched layers based on Jonathan Deakin [pressure](https://youtu.be/NEDQKOQtcG4), [pressure waveform](https://youtu.be/de8Kx8bF44I)

  - 1D wave equation with medium interfaces [pressure waveform](https://youtu.be/n6p0i72I_dQ)

  - 1D fluid structure interaction [sound pressure](https://www.youtube.com/watch?v=txpj6ah2WAw)