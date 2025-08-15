#!/bin/bash -f

wget https://ftp.gnu.org/gnu/octave/

OCTAVE_TAR=`gawk -f "octave-source.awk" "index.html"`

if test -z "${OCTAVE_TAR}"; then
    echo Octave release not found
    exit 1
fi

echo found "${OCTAVE_TAR}"

wget "https://ftp.gnu.org/gnu/octave/${OCTAVE_TAR}"

tar -zxvf "${OCTAVE_TAR}"
