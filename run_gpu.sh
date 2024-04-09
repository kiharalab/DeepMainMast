#!/bin/bash
echo "DeepMainMast started"
cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DeepMainMast"
#source /etc/profile.d/modules.sh
module load miniconda38
eval "$(conda shell.bash hook)"
source /apps/miniconda38/etc/profile.d/conda.sh
conda activate /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/deepmainmast

#Inputs
map_path=$1
fasta_path=$2
af_path=$3
contour=$4
output_dir=$5
TIME=$6
cpu=8
rosetta_path="/apps/rosetta/rosetta_bin_linux_2021.16.61629_bundle"
HOMO=$7
rosetta=$8
# Acitivate the conda enviroment
command_line="./dmm_gpu.sh -p ./ -m $map_path -f $fasta_path -c $contour -o $output_dir -t $TIME -T $TIME -C 8 -M 8 "
if [[ "$HOMO" -eq 1 ]]; then
        command_line+=" -H"
fi
if [[ "$af_path" != "" ]]; then
        command_line+=" -A $af_path "
fi
if [[ "$rosetta" -eq 1 ]]; then
            command_line+=" -x $rosetta_path"
fi
echo $command_line
$command_line

echo "INFO : DeepMainMast GPU Prediction Done"
