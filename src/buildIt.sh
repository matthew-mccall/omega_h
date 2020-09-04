#!/bin/bash -x
runNum=$1
cd /g/g19/smith516/develop/build-omega-rocm36/src 
#src=/g/g19/smith516/develop/omega_h/src/Omega_h_swap3d_topology.cpp
src=/g/g19/smith516/develop/omega_h/src/swap.cpp
rm foo.o
cp $src ${src}.${runNum}
/opt/rocm-3.6.0/bin/hipcc \
  -I/g/g19/smith516/develop/omega_h/src \
  -I/g/g19/smith516/develop/build-omega-rocm36/src \
  -I/g/g19/smith516/develop/omega_h/tpl \
  --amdgpu-target=gfx906   -std=c++11 \
  -o foo.o \
  -c $src &> ice.${runNum}
cat ice.${runNum}
cd -
