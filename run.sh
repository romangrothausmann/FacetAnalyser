#!/bin/bash

XSOCK=/tmp/.X11-unix
XAUTH=/tmp/.docker.xauth
touch $XAUTH
xauth nlist $DISPLAY | sed -e 's/^..../ffff/' | xauth -f $XAUTH nmerge -

docker run -it \
       --volume=$XSOCK:$XSOCK:rw \
       --volume=$XAUTH:$XAUTH:rw \
       --device=/dev/dri:/dev/dri \
       --env="XAUTHORITY=${XAUTH}" \
       --env="DISPLAY" \
       --env="QT_X11_NO_MITSHM=1" \
       --user="faUser" \
       -v ${PWD}:/tmp/images \
       --workdir=/tmp/images \
       registry.gitlab.com/romangrothausmann/facetanalyser/master \
       /opt/paraview/bin/paraview $@

rm $XAUTH # remove to avoid accumulation of xauth settings
