# syntax=docker/dockerfile:1
FROM ubuntu:24.10
ENV BUILD_DIR=/tmp/build

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
    #make check
    make install
EOT

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
      ./configure --disable-dependency-tracking CPPFLAGS="${CPPFLAGS}" --with-static-modules --enable-octave --disable-Werror LDFLAGS="${XFLAGS}" CXXFLAGS="${XFLAGS}" CFLAGS="${XFLAGS}" FCFLAGS="${XFLAGS}" F77FLAGS="${XFLAGS}" --enable-multithread --with-arpack --with-umfpack --with-klu --with-arpack --with-lapack --without-metis --with-mpi --with-trilinos --with-pardiso --with-suitesparseqr --with-qrupdate
    fi
    make -j8 all test install
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

    make install

    popd
EOT


WORKDIR /home/ubuntu
RUN rm -rf ${BUILD_DIR}
USER ubuntu
ENV LANG=en_US.UTF-8

## docker run -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY -h $HOSTNAME -v $HOME/.Xauthority:/home/ubuntu/.Xauthority --name mboct-fem-pkg1 mboct-fem-pkg

ENTRYPOINT [ "/usr/local/bin/octave", "--persist", "--eval", "pkg('load','mboct-fem-pkg')" ]
