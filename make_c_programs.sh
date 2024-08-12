#!/bin/bash

###Set up DeepMainmast fast version

set -euxo pipefail

if [ -n ${CONDA_PREFIX} ]; then
  export CFLAGS="-I${CONDA_PREFIX}/include ${CFLAGS:-}"
fi

if [ ! -e "./src/DAQscore_Unet_src" ]||[ ! -e "./server_bin" ];then
    echo "Can not find src directory."
    exit
fi

pushd ./src

for pg in VESPER_Power DAQscore_Unet MainmastC_UnetAssembleMtx  MainmastC_Unet_node MainmastC_UnetAF2  MainmastC_UnetChainAssign  VESPER_Power_colab MainmastC_UnetAF2Refine
do
    dir=${pg}_src
    pushd $dir
    make clean
    make
        if [ ! -e $pg ];then
            echo "Failed to compile $pg"
            exit
        fi
        cp $pg ../../server_bin/.
    popd
done

pushd pulchra_306
cc -O3 -o ../../server_bin/pulchra pulchra.c pulchra_data.c -lm
popd

popd

#Check Files:
for pg in pulchra VESPER_Power DAQscore_Unet MainmastC_UnetAssembleMtx  MainmastC_Unet_node MainmastC_UnetAF2  MainmastC_UnetChainAssign  VESPER_Power_colab MainmastC_UnetAF2Refine
do
    if [ ! -e ./server_bin/${pg} ];then
        echo "Missing $pg. Please check ./src/${pg}_src/"
        exit
    fi
done

echo "DONE: All Programs were compiled."
