# syntax=docker/dockerfile:1

## Copyright (C) 2024(-2024) Reinhard <octave-user@a1.net>

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

#####################################################################
## ABOUT THIRD PARTY SOFTWARE
#####################################################################

## Any docker image which is built from this Dockerfile may contain
## additional third party software in binary form. All third party
## software is distributed under their respective license and terms
## of use.

## Some of those copyright statements are included below.

## Despite of extensive research by the author,
## there is no guarantee that this section is really comprehensive.
## So, just in case that the reader may find any third-party software
## inside such an image, but the respective license statement is not yet
## included here, please inform the author about the license text
## to be included!

#####################################################################
## GNU-Octave
#####################################################################

## Octave is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <https://www.gnu.org/licenses/>.

#####################################################################
## NLopt
#####################################################################

## NLopt combines several free/open-source nonlinear optimization
## libraries by various authors.  See the COPYING, COPYRIGHT, and README
## files in the subdirectories for the original copyright and licensing
## information of these packages.

## The compiled NLopt library, i.e. the combined work of all of the
## included optimization routines, is licensed under the conjunction of
## all of these licensing terms.  Currently, the most restrictive terms
## are for the code in the "luksan" directory, which is licensed under
## the GNU Lesser General Public License (GNU LGPL), version 2.1 or
## later (see luksan/COPYRIGHT).

## That means that the compiled NLopt library is governed by the terms of
## the LGPL.

## ---------------------------------------------------------------------------

## Other portions of NLopt, including any modifications to the abovementioned
## packages, are licensed under the standard "MIT License:"

## Copyright (c) 2007-2011 Massachusetts Institute of Technology

## Permission is hereby granted, free of charge, to any person obtaining
## a copy of this software and associated documentation files (the
## "Software"), to deal in the Software without restriction, including
## without limitation the rights to use, copy, modify, merge, publish,
## distribute, sublicense, and/or sell copies of the Software, and to
## permit persons to whom the Software is furnished to do so, subject to
## the following conditions:

## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
## WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#####################################################################
## Trilinos
#####################################################################

## ************************************************************************
##
##            Trilinos: An Object-Oriented Solver Framework
##                 Copyright (2001) Sandia Corporation
##
##
## Copyright (2001) Sandia Corporation. Under the terms of Contract
## DE-AC04-94AL85000, there is a non-exclusive license for use of this
## work by or on behalf of the U.S. Government.  Export of this program
## may require a license from the United States Government.
##
## 1. Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.
##
## 2. Redistributions in binary form must reproduce the above copyright
## notice, this list of conditions and the following disclaimer in the
## documentation and/or other materials provided with the distribution.
##
## 3. Neither the name of the Corporation nor the names of the
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
## EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
## CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
## EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
## PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
## PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
## NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## NOTICE:  The United States Government is granted for itself and others
## acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
## license in this data to reproduce, prepare derivative works, and
## perform publicly and display publicly.  Beginning five (5) years from
## July 25, 2001, the United States Government is granted for itself and
## others acting on its behalf a paid-up, nonexclusive, irrevocable
## worldwide license in this data to reproduce, prepare derivative works,
## distribute copies to the public, perform publicly and display
## publicly, and to permit others to do so.
##
## NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
## OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
## ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
## RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
## INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
## THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
##

#####################################################################
## MBDyn
#####################################################################

## MBDyn (C) is a multibody analysis code.
## http://www.mbdyn.org
##
## Copyright (C) 1996-2023
##
## Pierangelo Masarati     <pierangelo.masarati@polimi.it>
## Paolo Mantegazza        <paolo.mantegazza@polimi.it>
##
## Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
## via La Masa, 34 - 20156 Milano, Italy
## http://www.aero.polimi.it
##
## Changing this copyright notice is forbidden.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation (version 2 of the License).
##
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##
## ----------------------------------------------------------------

#####################################################################
## NetCDF
#####################################################################

## Copyright 2018 Unidata

## Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

## 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

## 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

## 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## The NetCDF Copyright.

## \page copyright Copyright

## Copyright 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
## 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 University
## Corporation for Atmospheric Research/Unidata.

## Portions of this software were developed by the Unidata Program at the
## University Corporation for Atmospheric Research.

## Access and use of this software shall impose the following obligations
## and understandings on the user. The user is granted the right, without
## any fee or cost, to use, copy, modify, alter, enhance and distribute
## this software, and any derivative works thereof, and its supporting
## documentation for any purpose whatsoever, provided that this entire
## notice appears in all copies of the software, derivative works and
## supporting documentation.  Further, UCAR requests that the user credit
## UCAR/Unidata in any publications that result from the use of this
## software or in any product that includes this software, although this
## is not an obligation. The names UCAR and/or Unidata, however, may not
## be used in any advertising or publicity to endorse or promote any
## products or commercial entity unless specific written permission is
## obtained from UCAR/Unidata. The user also understands that
## UCAR/Unidata is not obligated to provide the user with any support,
## consulting, training or assistance of any kind with regard to the use,
## operation and performance of this software nor to provide the user
## with any updates, revisions, new versions or "bug fixes."

## THIS SOFTWARE IS PROVIDED BY UCAR/UNIDATA "AS IS" AND ANY EXPRESS OR
## IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL UCAR/UNIDATA BE LIABLE FOR ANY SPECIAL,
## INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
## FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
## NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
## WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.


#####################################################################
## Intel MKL
#####################################################################

## Intel Simplified Software License (Version October 2022)
## Use and Redistribution. You may use and redistribute the software, which is
## provided in binary form only, (the "Software"), without modification, provided the
## following conditions are met:
## * Redistributions must reproduce the above copyright notice and these terms of use
## in the Software and in the documentation and/or other materials provided with
## the distribution.
## * Neither the name of Intel nor the names of its suppliers may be used to endorse
## or promote products derived from this Software without specific prior written
## permission.
## * No reverse engineering, decompilation, or disassembly of the Software is
## permitted, nor any modification or alteration of the Software or its operation
## at any time, including during execution.
## No other licenses. Except as provided in the preceding section, Intel grants no
## licenses or other rights by implication, estoppel or otherwise to, patent,
## copyright, trademark, trade name, service mark or other intellectual property
## licenses or rights of Intel.
## Third party software. "Third Party Software" means the files (if any) listed in
## the "third-party-software.txt" or other similarly-named text file that may be
## included with the Software. Third Party Software, even if included with the
## distribution of the Software, may be governed by separate license terms, including
## without limitation, third party license terms, open source software notices and
## terms, and/or other Intel software license terms. These separate license terms
## solely govern Your use of the Third Party Software.
## DISCLAIMER. THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
## WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT ARE
## DISCLAIMED. THIS SOFTWARE IS NOT INTENDED FOR USE IN SYSTEMS OR APPLICATIONS WHERE
## FAILURE OF THE SOFTWARE MAY CAUSE PERSONAL INJURY OR DEATH AND YOU AGREE THAT YOU
## ARE FULLY RESPONSIBLE FOR ANY CLAIMS, COSTS, DAMAGES, EXPENSES, AND ATTORNEYS'
## FEES ARISING OUT OF ANY SUCH USE, EVEN IF ANY CLAIM ALLEGES THAT INTEL WAS
## NEGLIGENT REGARDING THE DESIGN OR MANUFACTURE OF THE SOFTWARE.
## LIMITATION OF LIABILITY. IN NO EVENT WILL INTEL BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
## NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
## PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
## WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
## No support. Intel may make changes to the Software, at any time without notice,
## and is not obligated to support, update or provide training for the Software.
## Termination. Your right to use the Software is terminated in the event of your
## breach of this license.
## Feedback. Should you provide Intel with comments, modifications, corrections,
## enhancements or other input ("Feedback") related to the Software, Intel will be
## free to use, disclose, reproduce, license or otherwise distribute or exploit the
## Feedback in its sole discretion without any obligations or restrictions of any
## kind, including without limitation, intellectual property rights or licensing
## obligations.
## Compliance with laws. You agree to comply with all relevant laws and regulations
## governing your use, transfer, import or export (or prohibition thereof) of the
## Software.
## Governing law. All disputes will be governed by the laws of the United States of
## America and the State of Delaware without reference to conflict of law principles
## and subject to the exclusive jurisdiction of the state or federal courts sitting
## in the State of Delaware, and each party agrees that it submits to the personal
## jurisdiction and venue of those courts and waives any objections. THE UNITED
## NATIONS CONVENTION ON CONTRACTS FOR THE INTERNATIONAL SALE OF GOODS (1980) IS
## SPECIFICALLY EXCLUDED AND WILL NOT APPLY TO THE SOFTWARE.

FROM ubuntu:latest

LABEL org.opencontainers.image.title="mboct-fem-pkg"
LABEL org.opencontainers.image.vendor="Reinhard"
LABEL com.docker.desktop.extension.icon="https://yt3.googleusercontent.com/HNK0rqhf4dCLnlP99g8c1odImPBZLM_QxLRb7yGBGDlDCbiiReQmoRZK4fuSU0XLC4yGApVS5Q=s160-c-k-c0x00ffffff-no-rj"
LABEL org.opencontainers.image.description="This package belongs to a suite of packages which can be used for pre- and postprocessing of flexible bodies in MBDyn (www.mbdyn.org) with GNU-Octave. It contains a general purpose structural Finite Element toolkit for linear statics and dynamics, which can be used to generate flexible body data for MBDyn's modal element and also hydrodynamic lubricated journal and slider plain bearings."
LABEL com.docker.extension.screenshots=[{"alt"="Fourbar","url"="https://i.ytimg.com/an_webp/d4i5AYPxsG4/mqdefault_6s.webp?du=3000&sqp=CI6PubUG&rs=AOn4CLDf2JH21dC1U_B8UytYHSNBP_LV9g"},{"alt"="Twisted bar","url"="https://i.ytimg.com/an_webp/I8HENx5mszA/mqdefault_6s.webp?du=3000&sqp=CICNubUG&rs=AOn4CLCr0TUUL4UZoynMQxoO7JMz27qNRg"},{"alt"="Cook's membrane","url"="https://i.ytimg.com/vi/rxQP8V4U0dE/hqdefault.jpg?sqp=-oaymwE2CPYBEIoBSFXyq4qpAygIARUAAIhCGAFwAcABBvABAfgBsASAAuADigIMCAAQARgTIB4ofzAP&rs=AOn4CLC7YAoIVIleUQWKwJ00c3tNb-1CBg"}]
LABEL org.opencontainers.image.source=https://github.com/octave-user
LABEL com.docker.extension.publisher-url="https://github.com/octave-user"
LABEL com.docker.extension.additional-urls=[{"title":"MBDyn","url":"https://www.mbdyn.org"},{"title":"MBDyn-GitLab","url":"https://public.gitlab.polimi.it/DAER/mbdyn"},{"title":"Fourbar","url":"https://www.youtube.com/watch?v=d4i5AYPxsG4"},{"title":"Twisted bar","url":"https://www.youtube.com/watch?v=I8HENx5mszA"},{"title":"Rolling ring","url":"https://www.youtube.com/watch?v=rxQP8V4U0dE"},{"title":"Cook's membrane","url":"https://www.youtube.com/watch?v=EAgejp4jQ00"},{"title":"videos","url":"https://www.youtube.com/@nonlinearmultibodydynamics6802"}]
LABEL org.opencontainers.image.licenses=GPL3

ENV SRC_DIR=/usr/local/src/
ENV LICENSE_DIR=/usr/local/share/license/
ENV BUILD_DIR=/tmp/build/
ENV TESTS_DIR=/tmp/tests/
ENV MBD_NUM_TASKS=4
ENV RUN_TESTS='\<mbdyn\>|\<octave\>|\<oct-pkg\>'
ENV RUN_CONFIGURE=no
ENV CXX=g++
ENV CC=gcc
ENV FC=gfortran
ENV F77=gfortran

WORKDIR ${SRC_DIR}
WORKDIR ${BUILD_DIR}

COPY . .

RUN sed 's/Types: deb/Types: deb deb-src/g' -i /etc/apt/sources.list.d/ubuntu.sources

RUN rm -f /etc/apt/apt.conf.d/docker-clean; echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' > /etc/apt/apt.conf.d/keep-cache

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get -yq update && apt-get -yq build-dep octave && \
    apt-get -yq install mercurial git trilinos-all-dev libopenmpi-dev \
    libnlopt-dev libhdf5-dev libginac-dev libatomic-ops-dev wget \
    libnetcdf-c++4-dev parallel cmake clang++-18

ENV GMSH_URL="http://www.gmsh.info/bin/Linux/"
ENV GMSH_VERSION="stable"
ENV GMSH_TAR="gmsh-${GMSH_VERSION}-Linux64.tgz"

WORKDIR ${SRC_DIR}/gmsh
WORKDIR ${BUILD_DIR}/gmsh

RUN --mount=type=cache,target=${BUILD_DIR}/gmsh,sharing=locked <<EOT bash
    if ! test -f "${GMSH_TAR}"; then
      if ! wget "${GMSH_URL}${GMSH_TAR}"; then
        exit 1
      fi
    fi

    if ! test -d gmsh-*.*.*-Linux64/bin; then
      if ! tar -zxvf "${GMSH_TAR}"; then
        exit 1
      fi
    fi

    if ! install gmsh-*.*.*-Linux64/bin/gmsh /usr/local/bin; then
      exit 1
    fi

    if ! gmsh --version; then
      exit 1
    fi
    cp "${GMSH_TAR}" "${SRC_DIR}/gmsh"
EOT

WORKDIR ${SRC_DIR}/tfel
WORKDIR ${BUILD_DIR}/tfel

ENV TFEL_REPO=https://github.com/thelfer/tfel.git
ENV LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

RUN --mount=type=cache,target=${BUILD_DIR}/tfel,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/tfel/.git; then
      if ! git clone ${TFEL_REPO} ${BUILD_DIR}/tfel; then
        exit 1
      fi
    fi

    if ! test -d build_dir; then
      mkdir build_dir
    fi

    cd build_dir

    if ! test -f Makefile; then
      if ! cmake -S .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-18 -DCMAKE_C_COMPILER=clang-18; then
        exit 1
      fi
    fi

    if ! make -j${MBD_NUM_TASKS}; then
      exit 1
    fi

    if echo tfel | awk "BEGIN{m=0;} /${RUN_TESTS}/ {m=1;} END{if(m==0) exit(1);}"; then
      if ! make -j${MBD_NUM_TASKS} check; then
        exit 1
      fi
    fi

    if ! make install; then
      exit 1
    fi

    if test -f install_manifest.txt; then
      cat install_manifest.txt | awk '/\/lib.*\.so$/' | xargs chmod +x
    fi
EOT

ENV MGIS_REPO=https://github.com/thelfer/MFrontGenericInterfaceSupport.git

WORKDIR ${SRC_DIR}/mgis
WORKDIR ${BUILD_DIR}/mgis

RUN --mount=type=cache,target=${BUILD_DIR}/mgis,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/mgis/.git; then
      if ! git clone ${MGIS_REPO} ${BUILD_DIR}/mgis; then
        exit 1
      fi
    fi

    if ! test -d build_dir; then
      mkdir build_dir
    fi

    cd build_dir

    if ! test -f Makefile; then
      if ! cmake -DCMAKE_BUILD_TYPE=Release -S ..; then
        exit 1
      fi
    fi

    if ! make -j${MBD_NUM_TASKS} all; then
      exit 1
    fi

    if echo mgis | awk "BEGIN{m=0;} /${RUN_TESTS}/ {m=1;} END{if(m==0) exit(1);}"; then
      if ! make -j${MBD_NUM_TASKS} check; then
        exit 1
      fi
    fi

    if ! make install; then
      exit 1
    fi

    if test -f install_manifest.txt; then
      cat install_manifest.txt | awk '/\/lib.*\.so$/' | xargs chmod +x
    fi
EOT

WORKDIR ${SRC_DIR}/gallery
WORKDIR ${BUILD_DIR}/gallery

ENV GALLERY_REPO=https://github.com/thelfer/MFrontGallery.git

RUN --mount=type=cache,target=${BUILD_DIR}/gallery,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/gallery/.git; then
      if ! git clone ${GALLERY_REPO} ${BUILD_DIR}/gallery; then
        exit 1
      fi
    fi

    if ! test -d build_dir; then
      mkdir build_dir
    fi

    cd build_dir

    if ! test -f Makefile; then
      if ! cmake -S .. -DCMAKE_BUILD_TYPE=Release -Denable-generic=ON -Denable-generic-behaviours=ON; then
        echo "Warning: cmake failed"
      fi
    fi

    if ! make -j${MBD_NUM_TASKS} all; then
      exit 1
    fi

    if echo gallery | awk "BEGIN{m=0;} /${RUN_TESTS}/ {m=1;} END{if(m==0) exit(1);}"; then
      if ! make -j${MBD_NUM_TASKS} check; then
        exit 1
      fi
    fi

    if ! make install; then
      exit 1
    fi

    if test -f install_manifest.txt; then
      cat install_manifest.txt | awk '/\/lib.*\.so$/' | xargs chmod +x
    fi
EOT

WORKDIR ${SRC_DIR}/octave
WORKDIR ${BUILD_DIR}/octave

RUN --mount=type=cache,target=${BUILD_DIR}/octave,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/octave/.hg; then
      if ! hg clone https://www.octave.org/hg/octave ${BUILD_DIR}/octave; then
        exit 1
      fi
    fi

    hg pull
    hg update
    hg checkout stable

    if ! test -x ./configure; then
      ./bootstrap
    fi

    if ! test "${RUN_CONFIGURE}" = no; then
      rm -f Makefile
    fi

    if ! test -f Makefile; then
      ./configure CXXFLAGS="-O3 -Wall -march=native"
    fi

    if ! make -j${MBD_NUM_TASKS} all; then
      exit 1
    fi

    if echo octave | awk "BEGIN{m=0;} /${RUN_TESTS}/ {m=1;} END{if(m==0) exit(1);}"; then
      if ! make check; then
        exit 1
      fi
    fi

    if ! make install; then
      exit 1
    fi

    if make dist-bzip2; then
      cp octave-*.tar.bz2 ${SRC_DIR}/octave
    fi
EOT

WORKDIR ${LICENSE_DIR}/mkl
RUN wget --output-document=mkl-license.pdf https://cdrdv2.intel.com/v1/dl/getContent/749362

WORKDIR ${BUILD_DIR}

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get -yq install libmkl-full-dev

# WORKDIR ${SRC_DIR}/siconos
# WORKDIR ${BUILD_DIR}/siconos

# ENV SICONOS_REPOSITORY="https://github.com/siconos/siconos.git"
# ENV SICONOS_BRANCH="master"

# RUN --mount=type=cache,target=${BUILD_DIR}/siconos,sharing=locked <<EOT bash
#     if ! test -d ${BUILD_DIR}/siconos/.git; then
#       git clone -b "${SICONOS_BRANCH}" "${SICONOS_REPOSITORY}" "${BUILD_DIR}/siconos"
#     fi

#     git pull --force

#     if ! test -d build_dir; then
#       mkdir build_dir
#     fi

#     cd build_dir

#     if ! test -f Makefile; then
#       export CPPFLAGS="`python3-config --cflags` ${CPPFLAGS}"
#       export LDFLAGS="`python3-config --ldflags` ${LDFLAGS}"
#       export LIBS="`python3-config --libs` ${LIBS}"
#       cmake -S .. -DSICONOS_SYSTEM_WIDE_INSTALL=ON -DWITH_TESTING=OFF  -DWITH_PYTHON_WRAPPER=OFF -DWITH_CMAKE_BUILD_TYPE=Release ## -DUSER_OPTIONS_FILE=../config_samples/siconos_lite.cmake
#     fi

#     make -j${MBD_NUM_TASKS} all
#     make install
# EOT

ENV MBD_FLAGS="-Ofast -Wall -march=native -mtune=native"
ENV MBD_CPPFLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/ompi/mpi/cxx -I/usr/include/x86_64-linux-gnu/openmpi -I/usr/include/trilinos -I/usr/include/suitesparse -I/usr/include/mkl -I/usr/local/include/MGIS -I/usr/local/include/MFront"
ENV MBD_CXXFLAGS="-std=c++20"
ENV MBD_ARGS_WITH="--with-mfront --with-static-modules --with-arpack --with-umfpack --with-klu --with-arpack --with-lapack --without-metis --with-mpi --with-trilinos --with-pardiso --with-suitesparseqr --with-qrupdate"
ENV MBD_ARGS_ENABLE="--enable-octave --enable-multithread --disable-Werror"

WORKDIR ${SRC_DIR}/mbdyn
WORKDIR ${BUILD_DIR}/mbdyn

RUN --mount=type=cache,target=${BUILD_DIR}/mbdyn,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/mbdyn/.git; then
      if ! git clone -b develop https://public.gitlab.polimi.it/DAER/mbdyn.git ${BUILD_DIR}/mbdyn; then
        exit 1
      fi
    fi

    git pull --force

    if ! test -x ./configure; then
      ./bootstrap.sh
    fi

    if ! test "${RUN_CONFIGURE}" = no; then
      rm -f Makefile
    fi

    if ! test -f Makefile; then
      ./configure CPPFLAGS="${MBD_CPPFLAGS}" LDFLAGS="${MBD_FLAGS}" CXXFLAGS="${MBD_FLAGS} ${MBD_CXXFLAGS}" CFLAGS="${MBD_FLAGS}" FCFLAGS="${MBD_FLAGS}" F77FLAGS="${MBD_FLAGS}" ${MBD_ARGS_WITH} ${MBD_ARGS_ENABLE}
    fi

    if ! make -j${MBD_NUM_TASKS} all; then
      exit 1
    fi

    if echo mbdyn | awk "BEGIN{m=0;} /${RUN_TESTS}/ {m=1;} END{if(m==0) exit(1);}"; then
      if ! make test; then
        exit 1
      fi
    fi

    if ! make install; then
      exit 1
    fi

    if ! /usr/local/mbdyn/bin/mbdyn --version; then
      exit 1
    fi

    if make dist; then
      cp mbdyn-*.tar.gz ${SRC_DIR}/mbdyn
    fi
EOT

WORKDIR ${SRC_DIR}/octave-pkg
WORKDIR ${BUILD_DIR}/octave-pkg

RUN --mount=type=cache,target=${BUILD_DIR}/octave-pkg,sharing=locked <<EOT bash
    octave -q --eval 'pkg install -verbose -global -forge nurbs;pkg install -verbose -global -forge netcdf'

    if ! test -d mboct-octave-pkg; then
      git clone -b master 'https://github.com/octave-user/mboct-octave-pkg.git'
    fi

    pushd mboct-octave-pkg && git pull --force && popd

    make CXXFLAGS="${MBD_FLAGS}" -j${MBD_NUM_TASKS} -C 'mboct-octave-pkg' dist install_global

    if ! test -d mboct-numerical-pkg; then
      git clone -b master 'https://github.com/octave-user/mboct-numerical-pkg.git'
    fi

    pushd mboct-numerical-pkg && git pull --force && popd

    make CXXFLAGS="${MBD_FLAGS}" -j${MBD_NUM_TASKS} -C 'mboct-numerical-pkg' dist install_global

    if ! test -d mboct-mbdyn-pkg; then
      git clone -b master 'https://github.com/octave-user/mboct-mbdyn-pkg.git'
    fi

    pushd mboct-mbdyn-pkg && git pull --force && popd

    make CXXFLAGS="${MBD_FLAGS}" -j${MBD_NUM_TASKS} -C 'mboct-mbdyn-pkg' dist install_global

    if ! test -d mboct-fem-pkg; then
      git clone -b master 'https://github.com/octave-user/mboct-fem-pkg.git'
    fi

    pushd mboct-fem-pkg && git pull --force && popd

    make CXXFLAGS="${MBD_FLAGS}" -j${MBD_NUM_TASKS} -C 'mboct-fem-pkg' dist install_global

    pushd mboct-fem-pkg/src

    if ! test -x configure; then
      ./bootstrap
    fi

    if ! test "${RUN_CONFIGURE}" = no; then
      rm -f Makefile
    fi

    if ! test -f Makefile; then
      ./configure CXXFLAGS="${MBD_FLAGS}"
    fi

    make -j${MBD_NUM_TASKS} all
    make install
    make distclean

    popd

    find . -name 'mboct-*-pkg-*.tar.gz' -exec cp '{}' ${SRC_DIR}/octave-pkg ';'
EOT

ENV OCT_PKG_LIST="netcdf:yes:master:yes:unlimited nurbs:yes:master:yes:unlimited mboct-octave-pkg:yes:master:yes:unlimited mboct-numerical-pkg:yes:master:yes:unlimited mboct-fem-pkg:yes:master:yes:unlimited mboct-mbdyn-pkg:yes:master:yes:unlimited"

WORKDIR ${TESTS_DIR}/octave-pkg-tests
WORKDIR ${BUILD_DIR}/mbdyn

RUN --mount=type=cache,target=${BUILD_DIR}/mbdyn,sharing=locked <<EOT bash
    if echo oct-pkg | awk "BEGIN{m=0;} /${RUN_TESTS}/ {m=1;} END{if(m==0) exit(1);}"; then
      if ! ${BUILD_DIR}/mbdyn/testsuite/octave_pkg_testsuite.sh --octave-pkg-test-dir ${TESTS_DIR}/octave-pkg-tests --octave-pkg-test-mode single; then
        exit 1
      fi
    fi
EOT

WORKDIR /home/ubuntu
RUN rm -rf ${BUILD_DIR} ${TESTS_DIR} ## Clean up temporary files
USER ubuntu
ENV LANG=C

## Run on Ubuntu with graphics enabled
## docker run -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY -h $HOSTNAME -v $HOME/.Xauthority:/home/ubuntu/.Xauthority --name mboct-fem-pkg1 mboct-fem-pkg:v2

# docker run \
#   --rm \
#   --network=host \
#   --env="DISPLAY" \
#   --env="HOME=$HOME" \
#   --env="XDG_RUNTIME_DIR=$XDG_RUNTIME_DIR" \
#   --user $(id -u):$(id -g) \
#   --volume="$HOME:$HOME:rw" \
#   --volume="/dev:/dev:rw" \
#   --volume="/run/user:/run/user:rw" \
#   --workdir="$HOME" \
#   ghcr.io/octave-user/mboct-fem-pkg

ENTRYPOINT [ "/usr/local/bin/octave", "--persist", "--eval", "pkg('load','mboct-fem-pkg')" ]
