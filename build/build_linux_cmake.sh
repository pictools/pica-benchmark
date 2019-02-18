#!/bin/sh

C_COMPILER=gcc
CXX_COMPILER=g++
LINKER=ld

NUM_CORES=$(grep processor /proc/cpuinfo | wc -l)

command_exists ()
{
    type $1 > /dev/null 2>&1;
}

if command_exists icc && command_exists icpc ; then
    C_COMPILER=icc
    CXX_COMPILER=icpc
    LINKER=icpc
fi

BUILD_DIR="unix_makefiles"
if [ ! -d $BUILD_DIR ]; then
    mkdir -p $BUILD_DIR
fi
cd $BUILD_DIR

CXX=$CXX_COMPILER CC=$C_COMPILER LD=$LINKER cmake -G "Unix Makefiles" ../..
make -j $NUM_CORES -k 2> /dev/null
if [ $? -ne 0 ]; then
  make
cd ..
