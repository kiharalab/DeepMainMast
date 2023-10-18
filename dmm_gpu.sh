#!/bin/bash
# only for server GPU running part
check_exists(){
	if [ ! -e "$1" ];then
		echo "Can not find $1"
		exit 1
	fi
}

check_not_exists(){
	if [ -e "$1" ];then
		echo "Already exists: $1"
		exit 1
	fi
}

check_empty(){
	if [ -z "$1" ];then
		echo "****"
		echo "Missing Value for $2"
		echo "****"
		show_help
		exit 1
	fi
}

function run_commands(){
	local command_list=("$@")
	local max_jobs=$MAX_Threads
	local job_count=0
	echo "max_jobs: $max_jobs"
	echo "commands: ${command_list[@]}"
	for command in "${command_list[@]}"; do
		while ((job_count>=max_jobs)); do
			sleep 5
			job_count=$(jobs -pr|wc -l)
		done
		#execute command
		echo "$command"
		eval "$command" &
		((job_count++))
	done

	wait
}

function run_subcommands(){
	local command_list=("$@")
	local max_jobs=$MAXSUB_Threads
	local job_count=0
	echo "max_jobs: $max_jobs"
	echo "commands: ${command_list[@]}"
	for command in "${command_list[@]}"; do
		while ((job_count>=max_jobs)); do
			sleep 5
			job_count=$(jobs -pr|wc -l)
		done
		#execute command
		echo "$command"
		eval "$command" &
		((job_count++))
	done

	wait
}

show_help() {
    echo "Usage: dmm_full_multithreads.sh [option] -f [FASTA file] -m [MAP file]"
    echo "DeepMaimast full-protocol"
    echo "Paratemters:"
    echo "  -h, Show help page"
    echo "  -A [AlphaFold Model]"
    echo "  -p [PATH of the programs]"
    echo "  -m [MAP file]"
    echo "  -f [Sequence File in FASTA format]"
    echo "  -c [Contour Level, density cut-off value]"
    echo "  -o [OutPut Directory Path]"
    echo "  -t [Time for main-chain tracing (sec)]"
    echo "  -T [Time for fragment assembly (sec)]"
    echo "  -C [Max number of threads]"
	echo "  -M [Max cpu cores per one thread job]"
    echo "  -r, Resume the job using an existing directory"
    echo "  -x [Rosetta program PATH], Execute RosettaCM modeling part"
    echo "  -H, Homo-oligomer targets. DeepMainmast refine chain-assignment."
	exit 1
}

#Parameters please edit-----------
## Number of Vehicles for main-chain tracing
Num_vehcles=(5 10 20 40)
## Number of Path-Sequence Alignments per path
Num_alignments=(5 10 20)
## Probability cut-off for CA atoms
Pca_cutoff=(0.3 0.4 0.5)
## Minimum number of aligned positions in the AF2 model
Aligned_positions=(30 50)
## Resolution to generate simulated density from AF2 models
sim_resolutions=(5.0 6.0 7.0)
## GPU



resume_flag=false
Buildfullatom_flag=false
MAX_Threads=8
Nsub=4
af2_mode=false
homo_mode=false

while getopts "A:p:m:f:c:o:t:T:C:rhM:x:H" option; do
	case $option in
		A)
			echo "Option -A: $OPTARG"
			af2_model_ori=$OPTARG
			check_exists $af2_model_ori
			af2_mode=true
			;;
		H)
			echo "Option -H: HomoOligomer Mode"
			homo_mode=true
			;;
		p)
			echo "Option -p: $OPTARG"
			PROGRAM_PATH=$OPTARG
			check_exists $PROGRAM_PATH
			;;
		m)
			echo "Option -m: $OPTARG"
			map=$OPTARG
			check_exists $map
			;;
		f)
			echo "Option -f: $OPTARG"
			fasta=$OPTARG
			check_exists $fasta
			;;
		c)
			echo "Option -c: $OPTARG"
			CONTOUR=$OPTARG
			;;
		o)
			echo "Option -o: $OPTARG"
			output_dir=$OPTARG
			;;
		t)
			echo "Option -t: $OPTARG"
			TIME_TRACE=$OPTARG
			;;
		T)
			echo "Option -T: $OPTARG"
			TIME_ASB=$OPTARG
			;;
		r)
			echo "Resume Existing Data"
			resume_flag=true
			;;
		C)
			echo "Option MAX Threads -C: $OPTARG"
			MAX_Threads=$OPTARG
			;;
		M)
			echo "Option Max SubThreads: $OPTARG"
			Nsub=$OPTARG
			;;
		x)
			echo "Option -x: $OPTARG"
			RosettaPath=$OPTARG
			check_exists $RosettaPath
			Buildfullatom_flag=true
			;;
		h)
			echo "**HELP**"
			show_help
			exit
			;;
		\?)
			echo "Invalid option: -$option"
			show_help
			exit
			;;
		:)
			echo "Option -$OPTARG needs an argument"
			show_help
			;;
	esac
done

echo "#PROGRAM_PATH: $PROGRAM_PATH"
echo "#map: $map"
echo "#fasta: $fasta"
echo "#CONTOUR: $CONTOUR"
echo "#output_dir: $output_dir"
echo "#TIME_TRACE: $TIME_TRACE"
echo "#TIME_Assemble: $TIME_ASB"
echo "#MAX_Threads: $MAX_Threads"
echo "#FullAtom Building: $RosettaPath"

((MAXSUB_Threads=$MAX_Threads/$Nsub)) #each sub-thread use 4 cores

check_empty "$PROGRAM_PATH" "-p"
check_exists $PROGRAM_PATH/server_bin/MainmastC_UnetAF2

if "${resume_flag}";then
	echo "Resume Process"
else
	check_not_exists $output_dir
fi

check_empty "$CONTOUR" "-c"
check_empty "$output_dir" "-o [output dir]"
check_empty "$TIME_TRACE" "-t [Time for Tracing (sec)]"
check_empty "$TIME_ASB" "-T [Time for Assembly (sec)]"
check_empty "$fasta" "-f [fasta formatfile]"

#Environment-------------
#PLEASE EDIT THIS SECTION ACCORDING TO 
#YOUR ENVIRONMENT
#module load miniconda38
#eval "$(conda shell.bash hook)"
## Acitivate the conda enviroment
#module load deepmainmast-multi
#------------------------
#Programs
BIN_DIR="$PROGRAM_PATH/server_bin/"
EMAP2SF=$BIN_DIR/emap2secplus/main_web.py
NODE=$BIN_DIR/MainmastC_Unet_node


###Check Rosetta programs
if "${Buildfullatom_flag}"; then
	echo "Checking $RosettaPath"
	ROSETTA3=$RosettaPath/main/
	RosettaScript=$ROSETTA3/source/bin/rosetta_scripts.static.linuxgccrelease
	check_exists $RosettaScript

	export ROSETTA3=$RosettaPath/main/
fi
#Set Up files
mkdir -p $output_dir/results/unet

UNET_DIR=$output_dir/results/unet
RESULTS_DIR=$output_dir/results
OUTF=$RESULTS_DIR
OUTCA=$OUTF/FINAL_CA_MODELs/
OUT_RANKED=$OUTF/RANKED_DATA/

if "${resume_flag}";then
	echo "INFO: Resuming..."
else
	#cp $map $RESULTS_DIR/input.map
	cp $map $RESULTS_DIR/input.mrc
	cp $fasta $RESULTS_DIR/seq.fasta
fi

SEQ=$RESULTS_DIR/seq.fasta
check_exists $SEQ




#Run Emap2sf----------------------------------
if [ ! -e $UNET_DIR/sigmoidAA_GLN.mrc ]||[ ! -e $UNET_DIR/atom_CB.mrc ];then
    #python3 $EMAP2SF --gpu 0 --mode=0 -O $RESULTS_DIR -F=$RESULTS_DIR/input.mrc -M=$BIN_DIR/emap2secplus/best_model/Unet_Protein_Atom.pth.tar --contour=$CONTOUR  --batch_size=16 --type=0
    python3 $EMAP2SF --mode=0 -O $RESULTS_DIR -F=$RESULTS_DIR/input.mrc -M=$BIN_DIR/emap2secplus/best_model/Unet_Protein_Atom.pth.tar --contour=$CONTOUR  --batch_size=16 --type=0
    echo "INFO : ATOM prediction Done"
    #python3 $EMAP2SF --gpu 0 --mode=0 -O $RESULTS_DIR -F=$RESULTS_DIR/input.mrc -M=$BIN_DIR/emap2secplus/best_model/Unet_Protein_AA.pth.tar --contour=$CONTOUR  --batch_size=16 --type=1
    python3 $EMAP2SF --mode=0 -O $RESULTS_DIR -F=$RESULTS_DIR/input.mrc -M=$BIN_DIR/emap2secplus/best_model/Unet_Protein_AA.pth.tar --contour=$CONTOUR  --batch_size=16 --type=1
    echo "INFO : AA prediction Done"
fi
#check output files
if [ ! -e $UNET_DIR/sigmoidAA_GLN.mrc ]||[ ! -e $UNET_DIR/atom_CB.mrc ];then
    echo "ERROR : Emap2sf failed"
    exit;
fi
