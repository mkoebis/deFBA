#!/bin/sh

# configure path settings for libsbml and c++ compilers 

#configure path setting for libSBML
export currentDir=$(pwd)
export LIBSBML=$currentDir/external/libSBML/libSBML-5.15.0-matlab

export LD_LIBRARY_PATH=$LIBSBML:$LD_LIBRARY_PATH

#make sure matlab uses the system std libs
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6:/lib/x86_64-linux-gnu/libgcc_s.so.1:/usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0


#matlab $@
/import/matlabR2016a/bin/matlab $@
