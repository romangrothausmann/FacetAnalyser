#!/bin/bash

XSOCK=/tmp/.X11-unix
XAUTH=/tmp/.docker.xauth
touch $XAUTH
xauth nlist $DISPLAY | sed -e 's/^..../ffff/' | xauth -f $XAUTH nmerge - # display number of currently visible VNC display (not working without VNC via vglconnect)
## VGL_DISPLAY: https://cdn.rawgit.com/VirtualGL/virtualgl/2.5.2/doc/index.html#hd0019001
VGL_DISPLAY=:0
xauth nlist $VGL_DISPLAY | sed -e 's/^..../ffff/' | xauth -f $XAUTH nmerge -  # display number running VGL X-Server on host

docker run -it \
       --runtime=nvidia \
       --volume=$XSOCK:$XSOCK:rw \
       --volume=$XAUTH:$XAUTH:rw \
       --volume="/usr/lib/x86_64-linux-gnu/libXv.so.1:/usr/lib/x86_64-linux-gnu/libXv.so.1" \
       --env="XAUTHORITY=${XAUTH}" \
       --env="DISPLAY" \
       --env="QT_X11_NO_MITSHM=1" \
       --user="faUser" \
       -v ${PWD}:/tmp/images \
       --workdir=/tmp/images \
       registry.gitlab.com/romangrothausmann/facetanalyser/vgl \
       vglrun /opt/paraview/bin/paraview $@

rm $XAUTH # remove to avoid accumulation of xauth settings
