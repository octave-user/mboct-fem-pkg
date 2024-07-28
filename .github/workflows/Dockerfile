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
    apt-get -yq install mercurial git trilinos-all-dev libopenmpi-dev libnlopt-dev libhdf5-dev libginac-dev libatomic-ops-dev wget


# RUN apt-get -yq build-dep octave
# RUN apt-get -yq install mercurial git trilinos-all-dev libopenmpi-dev libnlopt-dev libhdf5-dev libginac-dev libatomic-ops-dev wget

WORKDIR ${BUILD_DIR}/gmsh
ENV GMSH_URL="http://www.gmsh.info/bin/Linux/"
ENV GMSH_VERSION="stable"
ENV GMSH_TAR="gmsh-${GMSH_VERSION}-Linux64.tgz"

RUN --mount=type=tmpfs,target=${BUILD_DIR}/gmsh <<EOT bash
    wget "${GMSH_URL}${GMSH_TAR}"
    tar -zxvf "${GMSH_TAR}"
    install gmsh-*.*.*-Linux64/bin/gmsh /usr/local/bin
EOT

WORKDIR ${BUILD_DIR}/octave

RUN --mount=type=tmpfs,target=${BUILD_DIR}/octave --mount=type=cache,target=/usr/local,sharing=locked <<EOT bash
    hg clone https://www.octave.org/hg/octave ${BUILD_DIR}/octave
    hg checkout stable
    ./bootstrap && ./configure && make -j8 all && make check && make install
EOT

WORKDIR ${BUILD_DIR}

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get -yq install libmkl-full-dev

WORKDIR ${BUILD_DIR}/mbdyn

RUN --mount=type=tmpfs,target=${BUILD_DIR}/mbdyn --mount=type=cache,target=/usr/local,sharing=locked <<EOT bash
    git clone -b develop https://public.gitlab.polimi.it/DAER/mbdyn.git ${BUILD_DIR}/mbdyn
    ./bootstrap.sh
    XFLAGS="-flto=8 -fuse-linker-plugin -Ofast -Wall -march=native -mtune=native"
    CPPFLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/ompi/mpi/cxx -I/usr/include/x86_64-linux-gnu/openmpi -I/usr/include/trilinos -I/usr/include/suitesparse -I/usr/include/mkl"
    ./configure CPPFLAGS="${CPPFLAGS}" --with-static-modules --enable-octave --disable-Werror CXXFLAGS="${XFLAGS}" CFLAGS="${XFLAGS}" FCFLAGS="${XFLAGS}" F77FLAGS="${XFLAGS}" --enable-multithread --with-arpack --with-umfpack --with-klu --with-arpack --with-lapack --without-metis --with-mpi --with-trilinos --with-pardiso --with-suitesparseqr --with-qrupdate
    make -j8 all test install
EOT

WORKDIR ${BUILD_DIR}/octave-pkg

RUN --mount=type=tmpfs,target=${BUILD_DIR}/octave-pkg --mount=type=cache,target=/usr/local,sharing=locked <<EOT bash
    for pkg in octave numerical mbdyn fem; do
      git clone -b master https://github.com/octave-user/mboct-${pkg}-pkg.git
      make -C mboct-${pkg}-pkg install_global
    done
EOT

USER octave-user:octave-user

EXPOSE 6000

ENTRYPOINT [ "/usr/local/bin/octave", "--persist", "--eval", "pkg('load','mboct-fem-pkg')" ]
