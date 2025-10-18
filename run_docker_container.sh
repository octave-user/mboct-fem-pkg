#!/bin/sh -f

## Copyright (C) 2018(-2025) Reinhard <octave-user@a1.net>

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

IMAGE="octaveuser/mboct-fem-pkg:latest"
EXTRA_ARGS=""

while ! test -z "$1"; do
    case "$1" in
        --image)
            IMAGE="$2"
            shift
            ;;
        *)
            EXTRA_ARGS="${EXTRA_ARGS} $1"
            ;;
    esac
    shift
done

xhost + local:docker

docker run --entrypoint octave \
       --user root \
       -it \
       --rm \
       --network=host \
       --workdir="$HOME" \
       --volume="/run/user:/run/user:rw" \
       --volume="$HOME:$HOME:rw" \
       --volume=/tmp/.X11-unix:/tmp/.X11-unix:rw \
       --env="HOME=$HOME" \
       --env=DISPLAY=${DISPLAY} \
       --env="XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR}" \
       --env=LIBGL_ALWAYS_INDIRECT=0 \
       --env=LIBGL_ALWAYS_SOFTWARE=1 \
       --hostname=$HOSTNAME \
       --volume $HOME/.Xauthority:/home/ubuntu/.Xauthority \
       "${IMAGE}" \
       --persist \
       --norc \
       --eval 'pkg("local_list","/home/ubuntu/octave_packages");' \
       --eval 'pkg("load","mboct-fem-pkg");' \
       ${EXTRA_ARGS}

rc=$?

xhost -

exit ${rc}
