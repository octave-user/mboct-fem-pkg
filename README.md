# mboct-fem-pkg<sup>&copy;</sup>
**mboct-fem-pkg** belongs to a suite of packages which can be used for pre- and postprocessing of flexible bodies in MBDyn (https://www.mbdyn.org) with GNU-Octave (http://www.gnu.org/software/octave/) and Gmsh (http://www.gmsh.info/).
It contains a general purpose structural Finite Element toolkit for linear statics and dynamics, which can be used to generate flexible body data for MBDyn's modal element. 
In addition to that, there are functions for pre- and postprocessing of MBDyn's elastohydrodynamic journal- and slider- plain bearings.

Copyright<sup>&copy;</sup> 2019-2020

[Reinhard](mailto:octave-user@a1.net)

### MBDyn installation:
  - clone the source tree of MBDyn
`git clone https://public.gitlab.polimi.it/DAER/mbdyn.git -b develop`
  - compile and install MBDyn
`cd mbdyn && ./bootstrap.sh && ./configure --enable-octave --enable-autodiff --with-static-modules && make install`

## GNU Octave installation
  - Follow the instructions on (http://www.gnu.org/software/octave/) to install GNU Octave.

### GNU Octave package installation:
  - make sure that the GNU Octave nurbs package is installed
      `octave --eval 'pkg install -forge nurbs'`
  - install the following packages from github
      `for pkg in octave numerical mbdyn fem; do`
	  `git clone https://github.com/octave-user/mboct-${pkg}-pkg.git && make -C ${pkg} install_local`
      `done`

### Gmsh installation:
  - Follow the instructions on (http://www.gmsh.info/) to install Gmsh
  - Make sure Gmsh is on your path (e.g. `export PATH=/opt/gmsh/bin:${PATH}`)

### Usage
  - run Octave
    `octave`
  - at the Octave prompt load the package
    `pkg load mboct-fem-pkg`
  - at the Octave prompt execute a demo
    `demo fem_cms_export`
