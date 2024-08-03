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

## Nevertheless, despite of extensive research by the author,
## there is no guarantee that this section is really comprehensive.
## So, if the reader can find any third-party software inside such
## an image, but the respective license text is not included in this
## section, please inform the author of this Dockerfile about it!
## Then the missing license text will be added to this section by
## the author.

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

FROM ubuntu:24.04
ENV BUILD_DIR=/opt/source

WORKDIR ${BUILD_DIR}

COPY . .

RUN sed 's/Types: deb/Types: deb deb-src/g' -i /etc/apt/sources.list.d/ubuntu.sources

RUN rm -f /etc/apt/apt.conf.d/docker-clean; echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' > /etc/apt/apt.conf.d/keep-cache
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
  --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get -yq update && apt-get -yq build-dep octave && \
    apt-get -yq install mercurial git trilinos-all-dev libopenmpi-dev libnlopt-dev libhdf5-dev libginac-dev libatomic-ops-dev wget libnetcdf-c++4-dev

WORKDIR ${BUILD_DIR}/gmsh
ENV GMSH_URL="http://www.gmsh.info/bin/Linux/"
ENV GMSH_VERSION="stable"
ENV GMSH_TAR="gmsh-${GMSH_VERSION}-Linux64.tgz"

RUN --mount=type=cache,target=${BUILD_DIR}/gmsh,sharing=locked <<EOT bash
    if ! test -f "${GMSH_TAR}"; then
      wget "${GMSH_URL}${GMSH_TAR}"
    fi
    if ! test -d gmsh-*.*.*-Linux64/bin; then
      tar -zxvf "${GMSH_TAR}"
    fi
    install gmsh-*.*.*-Linux64/bin/gmsh /usr/local/bin
EOT

WORKDIR ${BUILD_DIR}/octave

RUN --mount=type=cache,target=${BUILD_DIR}/octave,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/octave/.hg; then
      hg clone https://www.octave.org/hg/octave ${BUILD_DIR}/octave
    fi

    hg pull && hg update && hg checkout stable

    if ! test -x ./configure; then
      ./bootstrap
    fi

    if ! test -f Makefile; then
      ./configure CXXFLAGS="-O3 -Wall -march=native"
    fi

    make -j8 all
    make check
    make install
    make distclean
EOT

WORKDIR /home/ubuntu/licenses
RUN wget -o mkl-license.pdf https://cdrdv2.intel.com/v1/dl/getContent/749362

WORKDIR ${BUILD_DIR}

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get -yq install libmkl-full-dev

WORKDIR ${BUILD_DIR}/mbdyn

ENV XFLAGS="-Ofast -Wall -march=native -mtune=native"
ENV CPPFLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/ompi/mpi/cxx -I/usr/include/x86_64-linux-gnu/openmpi -I/usr/include/trilinos -I/usr/include/suitesparse -I/usr/include/mkl"

RUN --mount=type=cache,target=${BUILD_DIR}/mbdyn,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/mbdyn/.git; then
      git clone -b develop https://public.gitlab.polimi.it/DAER/mbdyn.git ${BUILD_DIR}/mbdyn
    fi
    git pull --force
    if ! test -x ./configure; then
      ./bootstrap.sh
    fi
    if ! test -f Makefile; then
      ./configure CPPFLAGS="${CPPFLAGS}" --with-static-modules --enable-octave --disable-Werror LDFLAGS="${XFLAGS}" CXXFLAGS="${XFLAGS}" CFLAGS="${XFLAGS}" FCFLAGS="${XFLAGS}" F77FLAGS="${XFLAGS}" --enable-multithread --with-arpack --with-umfpack --with-klu --with-arpack --with-lapack --without-metis --with-mpi --with-trilinos --with-pardiso --with-suitesparseqr --with-qrupdate
    fi
    make -j8 all
    make test
    make install
    make distclean
EOT

WORKDIR ${BUILD_DIR}/octave-pkg

RUN --mount=type=cache,target=${BUILD_DIR}/octave-pkg,sharing=locked <<EOT bash
    octave --eval 'pkg install -global -forge nurbs'
    octave --eval 'pkg install -global -forge netcdf'
    if ! test -d mboct-octave-pkg; then
      git clone -b master 'https://github.com/octave-user/mboct-octave-pkg.git'
    fi
    pushd mboct-octave-pkg && git pull --force && popd
    make CXXFLAGS="${XFLAGS}" -C 'mboct-octave-pkg' install_global
    if ! test -d mboct-numerical-pkg; then
      git clone -b master 'https://github.com/octave-user/mboct-numerical-pkg.git'
    fi
    pushd mboct-numerical-pkg && git pull --force && popd
    make CXXFLAGS="${XFLAGS}" -C 'mboct-numerical-pkg' install_global
    if ! test -d mboct-mbdyn-pkg; then
      git clone -b master 'https://github.com/octave-user/mboct-mbdyn-pkg.git'
    fi
    pushd mboct-mbdyn-pkg && git pull --force && popd
    make CXXFLAGS="${XFLAGS}" -C 'mboct-mbdyn-pkg' install_global
    if ! test -d mboct-fem-pkg; then
      git clone -b master 'https://github.com/octave-user/mboct-fem-pkg.git'
    fi
    pushd mboct-fem-pkg && git pull --force && popd
    make CXXFLAGS="${XFLAGS}" -C 'mboct-fem-pkg' install_global
    pushd mboct-fem-pkg/src
    if ! test -x configure; then
      ./bootstrap
    fi

    if ! test -f Makefile; then
      ./configure CXXFLAGS="${XFLAGS}"
    fi

    make -j8 all
    make install
    make distclean

    popd
EOT


WORKDIR /home/ubuntu
USER ubuntu
ENV LANG=en_US.UTF-8

## docker run -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY -h $HOSTNAME -v $HOME/.Xauthority:/home/ubuntu/.Xauthority --name mboct-fem-pkg1 mboct-fem-pkg

ENTRYPOINT [ "/usr/local/bin/octave", "--persist", "--eval", "pkg('load','mboct-fem-pkg')" ]