###Set up DeepMainmast fast version




if [ ! -e "./src/DAQscore_Unet_src" ]||[ ! -e "./server_bin" ];then
    echo "Can not find src directory."
    exit
fi

pushd ./src

for pg in VESPER_Power DAQscore_Unet MainmastC_UnetAssembleMtx  MainmastC_Unet_node MainmastC_UnetAF2  MainmastC_UnetChainAssign  VESPER_Power_colab
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


popd