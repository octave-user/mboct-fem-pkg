#!/bin/bash -f

shopt -s execfail

rm -f "index.html"

wget https://ftp.gnu.org/gnu/octave/

OCTAVE_TAR=`gawk -f "octave-source.awk" "index.html"`

if test -z "${OCTAVE_TAR}"; then
    echo Octave release not found
    exit 1
fi

echo found "${OCTAVE_TAR}"

wget "https://ftp.gnu.org/gnu/octave/${OCTAVE_TAR}"

cp "${OCTAVE_TAR}" "${SRC_DIR}/octave"

tar -zxvf "${OCTAVE_TAR}"

pushd $(basename -s .tar.gz "${OCTAVE_TAR}")

case "${RUN_CONFIGURE}" in
    *octave*|all)
        rm -f Makefile
        ;;
    none)
        ;;
esac

if ! test -f Makefile; then
    ./configure CXXFLAGS="-O3 -Wall -march=native" --with-hdf5-includedir=`pkg-config --cflags-only-I hdf5-serial | sed 's/^-I//'` --with-hdf5-libdir=`pkg-config --libs-only-L hdf5-serial | sed 's/^-L//'`
fi

make -j${MBD_NUM_TASKS} all

case "${RUN_TESTS}" in
    *octave*|all)
        if ! make check; then
            exit 1
        fi
        ;;
    none)
        ;;
esac

if ! make install; then
    exit 1
fi
