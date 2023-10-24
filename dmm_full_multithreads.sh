#!/bin/bash

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

function run_single_command(){
	local command_list=("$@")
	local max_jobs=1
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
    echo "  -C [Max number of CPU cores]"
	echo "  -M [Max cpu cores per one task]"
    #echo "  -r, Resume the job using an existing directory"
    echo "  -x [Rosetta program PATH], Execute RosettaCM modeling part"
    echo "  -H, Homo-oligomer targets. DeepMainmast refine chain-assignment."
    echo "  -F, Using reduced parameters for fast computing"
    echo "  -s, Using server mode. Only Generate 1 full-atom model/initial model."
	exit 1
}

resume_flag=false
Buildfullatom_flag=false
MAX_Threads=8
Nsub=4
af2_mode=false
homo_mode=false
fast_mode=false
server_mode=false

while getopts "A:p:m:f:c:o:t:T:C:rhM:x:HFs" option; do
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
		F)
			echo "Option -F: Using Reduced Parameters. FAST Mode"
			fast_mode=true
			;;
		s)
			echo "Option -s: Generate one full atom model/initial CA model"
			server_mode=true
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
			echo "Option MAX Number of CPU cores -C: $OPTARG"
			MAX_Threads=$OPTARG
			;;
		M)
			echo "Option Max SubThreads per task: $OPTARG"
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

#Parameters please edit-----------

#Full parameters
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

#Fast reduced parameters
if "${fast_mode}";then
 Num_vehcles=(5 10 20)
 Num_alignments=(5 10)
 Pca_cutoff=(0.3) #reduced
 Aligned_positions=(30 50)
 sim_resolutions=(5.0)
fi



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

#rand_tag=`python -c "import random, string; print(''.join(random.choices(string.ascii_letters + string.digits, k=8)))"`
rand_tag=''
RESULTS_DIR=$output_dir/results

#Attempt 3 times
if [ -e $RESULTS_DIR ];then
	rand_tag=`python -c "import random, string; print(''.join(random.choices(string.ascii_letters + string.digits, k=8)))"`
	RESULTS_DIR=$output_dir/results${rand_tag}
fi
if [ -e $RESULTS_DIR ];then
	rand_tag=`python -c "import random, string; print(''.join(random.choices(string.ascii_letters + string.digits, k=8)))"`
	RESULTS_DIR=$output_dir/results${rand_tag}
fi
if [ -e $RESULTS_DIR ];then
	echo "Can not generate new dir: $RESULTS_DIR"
	exit
fi

mkdir -p $RESULTS_DIR/unet
mkdir -p $RESULTS_DIR/unet
UNET_DIR=$RESULTS_DIR/unet
OUTF=$RESULTS_DIR
OUTCA=$OUTF/FINAL_CA_MODELs/
OUT_RANKED=$OUTF/RANKED_DATA/

cp $map $RESULTS_DIR/input.mrc

#Check non-standerd characters and coding
python3 $BIN_DIR/format_seq.py $fasta --OutPath $RESULTS_DIR/seq.fasta
#cp $fasta $RESULTS_DIR/seq.fasta

SEQ=$RESULTS_DIR/seq.fasta
check_exists $SEQ

##AF2 Model check
if "${af2_mode}";then
	#Renumber and assign chain ID
	RENUM=$BIN_DIR/Renum_chain.py
	af2_model=$RESULTS_DIR/AF2_renum.pdb
	if [ ! -e $af2_model ];then
		python3 $RENUM $SEQ $af2_model_ori --OutPath $af2_model
	fi
	check_exists $af2_model
fi


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

com_list=()

#C-alpha Node Detection-----------------------
for p in ${Pca_cutoff[@]}
do
	out=$RESULTS_DIR/NODE_p$p.pdb
	if [ ! -e $out ];then
		cmd="$NODE -i $UNET_DIR -s $SEQ -t $p -G > $out"
		com_list+=("$cmd")
	fi
done

#echo "$com_list"
run_commands "${com_list[@]}"
echo "INFO : C-alpha Node computing Done"

#Path tracing by VRP solver----------------###
com_list=()
ASB=$BIN_DIR/UnetPathEdge.py #Large penalties for edges with Pbb=0
for p in ${Pca_cutoff[@]}
do
	for Nchain in ${Num_vehcles[@]}
	do
		in=$RESULTS_DIR/NODE_p${p}.pdb
		out=$RESULTS_DIR/PATH_p${p}Nch${Nchain}.csv
		log=$RESULTS_DIR/PATH_p${p}Nch${Nchain}.log

		if [ ! -e $log ];then
		    cmd="python $ASB $in --Nchain $Nchain --OutCSV $out --SecTrace $TIME_TRACE  &> $log"
			com_list+=("$cmd")
		fi
	done
done
run_commands "${com_list[@]}"
echo "INFO : Path Tracing by VRP Done"

##Generate INP file
com_list=()
PG=$BIN_DIR/MainmastC_UnetAF2
ASB=$BIN_DIR/Assemble_Iter.py
for p in ${Pca_cutoff[@]}
do
	for Nchain in ${Num_vehcles[@]}
	do
        in=$OUTF/PATH_p${p}Nch${Nchain}.csv
		UnetDir=$UNET_DIR

		if [ ! -e $in ];then
			echo "Missing.... $in"
			continue
		fi

		for Nali in ${Num_alignments[@]}
		do
			INP=$OUTF/INP_p${p}Nch${Nchain}Nali${Nali}.txt
			if [ ! -e $INP ];then
			    cmd="$PG -i $UnetDir -s $SEQ -G -N $in -l 9 -r 1.5 -c 1 -T $Nali > $INP"
			 	com_list+=("$cmd")  
			fi
			if "${af2_mode}";then
				for Rsize in ${Aligned_positions[@]};do
                	INP=$OUTF/INP_p${p}Nch${Nchain}Nali${Nali}R${Rsize}.txt
					if [ ! -e $INP ];then
						cmd="$PG -i $UnetDir -s $SEQ -G -N $in -l 9 -r 1.5 -c 1 -T $Nali -R $Rsize -A $af2_model > $INP"
						com_list+=("$cmd")
					fi
				done
			fi
		done
    done
done
run_commands "${com_list[@]}"
echo "INFO : generating INP file Done"

##Assemble Fragments
com_list=()
ASB=$BIN_DIR/Assemble_Iter.py
OUTF=$RESULTS_DIR

for p in ${Pca_cutoff[@]}
do
	for Nchain in ${Num_vehcles[@]}
	do
        in=$OUTF/PATH_p${p}Nch${Nchain}.csv
		UnetDir=$UNET_DIR

		if [ ! -e $in ];then
			echo "Missing.... $in"
			continue
		fi

		for Nali in ${Num_alignments[@]}
		do
			INP=$OUTF/INP_p${p}Nch${Nchain}Nali${Nali}.txt
			OUT=$OUTF/OUT_p${p}Nch${Nchain}Nali${Nali}.pdb
			log=$OUTF/OUT_p${p}Nch${Nchain}Nali${Nali}.log
			if [ -e $INP ] && [ ! -e $OUT ];then
				cmd="python $ASB --Nstock 10000 --Niter 3 $INP --OutPath $OUT --Ncpu $Nsub --SecAssemble $TIME_ASB &> $log"
			 	com_list+=("$cmd")  
			fi
			if "${af2_mode}";then
				for Rsize in ${Aligned_positions[@]};do
                	INP=$OUTF/INP_p${p}Nch${Nchain}Nali${Nali}R${Rsize}.txt
					OUT=$OUTF/OUT_p${p}Nch${Nchain}Nali${Nali}R${Rsize}.pdb
					log=$OUTF/OUT_p${p}Nch${Nchain}Nali${Nali}R${Rsize}.log
					if [ -e $INP ] && [ ! -e $OUT ];then
						cmd="python $ASB --Nstock 10000 --Niter 3 $INP --OutPath $OUT --Ncpu $Nsub --SecAssemble $TIME_ASB &> $log"
						com_list+=("$cmd")
					fi
				done
			fi
		done
    done
done
run_subcommands "${com_list[@]}"
echo "INFO : Assembling Done"

#VESPER DOCKING--------
if "${af2_mode}";then
	##Set Up simulated maps
	VESPER_DIR=$OUTF/VESPER_MODELs
	PG1=$BIN_DIR/SplitPDBchain.py
	PG2=$BIN_DIR/btpdb2mrc.py 
	PG3=$BIN_DIR/VESPER_Power_colab
	PG4=$BIN_DIR/RotateVesper.py

	mkdir -p $VESPER_DIR

	#Split PDB to chains
	python $PG1 $af2_model --OutPath=$VESPER_DIR/af2

	com_list=()
	for reso in ${sim_resolutions[@]}
	do
		for input_model in $VESPER_DIR/af2_?.pdb $VESPER_DIR/af2_??.pdb
		do
			model_id=`basename $input_model .pdb`
			map_model=$VESPER_DIR/${model_id}_R${reso}.mrc
			if [ ! -e $map_model ]&& [ -e $input_model ];then
				cmd="python $PG2 $reso $input_model $map_model"
				com_list+=("$cmd")
			fi
		done
	done
	run_commands "${com_list[@]}"

	##VESPER DOCKING
	com_list=()
	for reso in ${sim_resolutions[@]}
	do
		for input_model in $VESPER_DIR/af2_?.pdb $VESPER_DIR/af2_??.pdb
		do
			model_id=`basename $input_model .pdb`
			map_model=$VESPER_DIR/${model_id}_R${reso}.mrc
			outfile=$VESPER_DIR/${model_id}_R${reso}.out
			if [ ! -e $outfile ] && [ -e $map_model ] && [ -e $map ];then
				#cmd="$PG3 -a $map -b $map_model -t $CONTOUR -T 10.0 -c $Nsub -g 8.0 -s 1.0 > $outfile"
				cmd="$PG3 -a $map -b $map_model -t $CONTOUR -T 10.0 -c 2 -g 8.0 -s 1.0 > $outfile"
				com_list+=("$cmd")
			fi			
		done
	done
	run_single_commands "${com_list[@]}"

	##Generate Fitted Models
	for reso in ${sim_resolutions[@]}
	do
		for input_model in $VESPER_DIR/af2_?.pdb $VESPER_DIR/af2_??.pdb
		do
			model_id=`basename $input_model .pdb`
			outfile=$VESPER_DIR/${model_id}_R${reso}.out
			first_model=$VESPER_DIR/${model_id}_R${reso}_FIT_MODEL0.pdb
			if [ -e $outfile ] && [ ! -e $first_model ];then
				python3 $PG4 $outfile $input_model --OutPath $VESPER_DIR/${model_id}_R${reso}_
			fi
		done
	done

	echo "INFO : VESPER DOCKING Done"
fi

##COMBINATION===================================
echo "INFO : Combine All Models"
PG_ASB=$BIN_DIR/MainmastC_UnetAssembleMtx
PG=$BIN_DIR/Assemble_Iter.py #combine all models
comb_model2=$OUTF/COMBINEi_DMonly.pdb
log=$OUTF/COMBINEi_DMonly.log

#Generate MODLES file

modelfile=$OUTF/MODELs_DMonly.pdb
if [ ! -e $modelfile ];then
	for Nali in ${Num_alignments[@]}; do
	for p in $OUTF/OUT_p*ali$Nali.pdb; do
		cat $p
		echo "TER"
	done;done > $modelfile
fi

if "${af2_mode}";then
	modelfile=$OUTF/MODELs_all.pdb
	if [ ! -e $modelfile ];then
		for p in $OUTF/OUT_p*.pdb
		do
			cat $p
			echo "TER"
		done > $modelfile
	fi
	modelfile=$OUTF/MODELs_AFonly.pdb
	if [ ! -e $modelfile ];then
		if [ -e $OUTF/VESPER_MODELs ];then
			for p in $OUTF/OUT_p*R*.pdb
			do
				cat $p
				echo "TER"
			done > $modelfile
		else
			echo "Running without AF model output"
		fi	
	fi
	modelfile=$OUTF/MODELs_VESPER.pdb
	if [ ! -e $modelfile ];then
		if [ -e $OUTF/VESPER_MODELs ];then
			for p in $OUTF/VESPER_MODELs/*FIT_MODEL?.pdb
			do
				cat $p
				echo "TER"
			done > $modelfile
		else
			echo "Running without AF model output"
		fi
	fi
fi


#generate MTX files
PG_ASB=$BIN_DIR/MainmastC_UnetAssembleMtx
com_list=()
for method in DMonly all AFonly VESPER; do
	mtxfile=$OUTF/MTX_${method}.txt
	modelfile=$OUTF/MODELs_${method}.pdb
	if [ ! -e $mtxfile ] && [ -e $modelfile ];then
		cmd="$PG_ASB  -i $UnetDir -s $SEQ -l 10 -A $modelfile -G -c $Nsub > $mtxfile"
		com_list+=("$cmd")
	fi
done
run_subcommands "${com_list[@]}"

#COMBINE all models
com_list=()
for method in DMonly all AFonly VESPER; do
	mtxfile=$OUTF/MTX_${method}.txt
	comb_model=$OUTF/COMBi_${method}.pdb
	log=$OUTF/COMBi_${method}.log
	if [ -e $mtxfile ] && [ ! -e $comb_model ];then
		cmd="python $PG $mtxfile --Nstock 10000 --OutPath $comb_model --SecAssemble $TIME_ASB --Ncpu $Nsub > $log"
		com_list+=("$cmd")
	fi
done

run_subcommands "${com_list[@]}"

echo "INFO : Combine All Models Done"

if "${homo_mode}";then
	PG_ASB=$BIN_DIR/MainmastC_UnetChainAssign
	PG=$BIN_DIR/Assemble_ReChainID.py
	Filter=-10.0
	DiffCut=3.0
	com_list=()
	for method in DMonly all AFonly VESPER; do
		mtxfile=$OUTF/CH_MTX_${method}.txt
		modelfile=$OUTF/COMBi_${method}.pdb
		if [ ! -e $mtxfile ] && [ -e $modelfile ];then
			cmd="$PG_ASB -i $UnetDir -s $SEQ -A $modelfile -H -f $Filter -r $DiffCut > $mtxfile"
			com_list+=("$cmd")
		fi
	done
	run_commands "${com_list[@]}"

	#chain assignment
	com_list=()
	for method in DMonly all AFonly VESPER; do
		mtxfile=$OUTF/CH_MTX_${method}.txt
		output_file=$OUTF/COMBi_${method}_rechain.pdb
		log=$OUTF/COMBi_${method}_rechain.log
		if [ -e $mtxfile ] && [ ! -e $output_file ];then
			cmd="python $PG $mtxfile --HomWeight 1.0 --Nstock 10000 --OutPath $output_file --SecAssemble $TIME_ASB --Ncpu $Nsub > $log"
			com_list+=("$cmd")
		fi
	done
	run_subcommands "${com_list[@]}"
	echo "INFO : Chain-Assignments Done"
fi

##Make Copies
mkdir -p $OUTCA

for method in DMonly all AFonly VESPER; do
	ca_file=$OUTF/COMBi_${method}.pdb
	if "${homo_mode}";then
		ca_file=$OUTF/COMBi_${method}_rechain.pdb
	fi
	outfile=$OUTCA/COMBi_${method}.pdb
	if [ -e $ca_file ];then
		cp $ca_file $outfile
	fi

done

echo "INFO : DAQ scoring for CA models"
DAQ=$BIN_DIR/DAQscore_Unet
com_list=()
for method in DMonly all AFonly VESPER; do
	comb_model=$OUTCA/COMBi_${method}.pdb
	DAQ_OUT=$OUTCA/COMBi_${method}.daq
	if [ -e $comb_model ] && [ ! -e $DAQ_OUT ];then
		cmd="$DAQ -i $UnetDir -Q $comb_model > $DAQ_OUT"
		com_list+=("$cmd")
	fi
done
run_commands "${com_list[@]}"

#WindowSize 9
#WINDOW_PG=$BIN_DIR/DAQwindow.py 
#for method in DMonly all AFonly VESPER; do
#	DAQ_OUT=$OUTCA/COMBINEi_${method}.daq
#	if [ -e $DAQ_OUT ];then
##		python3 $WINDOW_PG $DAQ_OUT 9 --OutPath $OUTCA/COMBINEi_${method}_
#	fi
#done

echo "INFO : DAQ scoring Done"

##DOT Scoring
PG=$BIN_DIR/VESPER_Power_colab
for method in DMonly all AFonly VESPER; do
	ca_model=$OUTCA/COMBi_${method}.pdb
	ca_dot=$OUTCA/COMBi_${method}.dot
	new_map=$OUTCA/COMBi_${method}.mrc
	if [ -e $ca_model ] && [ ! -e $ca_dot ];then
		#simulate map
		python3 $BIN_DIR/btpdb2mrc.py 5.0 $ca_model $new_map 
        $PG -a $map -b $new_map -t 0.001 -T 10.0 -c 4 -g 8.0 -e -s 2.0 > $ca_dot
	fi
done
##Ranking CA models
echo "INFO : DOT scoring Done"
echo "INFO : CA modeling Done"

python $BIN_DIR/ca_ranker.py $OUTCA

##Window19 models
PG=$BIN_DIR/DAQwindow.py
for mid in {1..4};do #Max 4 CA models
	daq_file=$OUTCA/rank${mid}.daq
	if [ -e $daq_file ];then
		python $PG $daq_file 9 --OutPath $OUTCA/rank${mid}_
	fi
done

if ! "${Buildfullatom_flag}";then
	echo "INFO : DeepMainmast Computation Done: CA models"
	check_exists $OUTCA/rank1_daq_score_w9.pdb
	cp $OUTCA/rank1_daq_score_w9.pdb $output_dir/DeepMainmast${rand_tag}.pdb
	echo "DONE" >$output_dir/done.out
	exit
fi

##Rosetta CM part---------------------
##SetUp Rosetta files
echo "INFO : Start RosettaCM"
PG=$BIN_DIR/RosettaCM.py
for method in DMonly all AFonly VESPER; do
	ca_model=$OUTCA/COMBi_${method}.pdb
	output_dir=$OUTCA/CM_${method}/
	if [ -e $ca_model ] && [ ! -e $output_dir ];then
		python3 $PG $SEQ $ca_model $map --OutPath=$output_dir --XMLPath=$BIN_DIR --PulchraPath=$BIN_DIR/pulchra
	fi
done

com_list=()
for method in DMonly all AFonly VESPER; do
	ca_model=$OUTCA/COMBi_${method}.pdb
	output_dir=$OUTCA/CM_${method}/
	if [ -e $output_dir ] && [ -e $output_dir/A_setup.sh ] && [ ! -e $output_dir/1tmpA_thread.pdb ];then
		cmd="bash $output_dir/A_setup.sh $output_dir/"
		com_list+=("$cmd")
	fi
done
run_commands "${com_list[@]}"

#Run RosettaCM, generate five models for each method
com_list=()
for method in DMonly all AFonly VESPER; do
	ca_model=$OUTCA/COMBi_${method}.pdb
	output_dir=$OUTCA/CM_${method}/
	log=$OUTCA/CM_${method}.log
	if [ -e $output_dir ] && [ -e $output_dir/A_setup.sh ] && [ -e $output_dir/C_rosettaCM.sh ] && [ ! -e $log ];then
		cmd="bash $output_dir/C_rosettaCM.sh $output_dir 5 > $log"
		if "${server_mode}";then #Generate one model
			cmd="bash $output_dir/C_rosettaCM.sh $output_dir 1 > $log"
		fi
		com_list+=("$cmd")
	fi
done
run_commands "${com_list[@]}"

#Evaluation

echo "INFO : DAQ scoring for Full-atom models"
DAQ=$BIN_DIR/DAQscore_Unet
com_list=()
for method in DMonly all AFonly VESPER; do
	for mid in 1 2 3 4 5; do
		FullAtomModel=$OUTCA/CM_${method}/S_singletgt_000${mid}.pdb
		DAQ_OUT=$OUTCA/CM_${method}/S_singletgt_000${mid}.daq
		echo "$DAQ_OUT"
		if [ -e $FullAtomModel ] && [ ! -e $DAQ_OUT ];then
			cmd="$DAQ -i $UnetDir -Q $FullAtomModel > $DAQ_OUT"
			com_list+=("$cmd")
		fi
	done
done
run_commands "${com_list[@]}"


##DOT Scoring
PG=$BIN_DIR/VESPER_Power_colab
for method in DMonly all AFonly VESPER; do
	for mid in 1 2 3 4 5; do
		FullAtomModel=$OUTCA/CM_${method}/S_singletgt_000${mid}.pdb
		FullDot=$OUTCA/CM_${method}/S_singletgt_000${mid}.dot
		new_map=$OUTCA/CM_${method}/S_singletgt_000${mid}.mrc
		if [ -e $FullAtomModel ] && [ ! -e $FullDot ];then
			#simulate map
			python3 $BIN_DIR/btpdb2mrc.py 5.0 $FullAtomModel $new_map 
        	$PG -a $map -b $new_map -t 0.001 -T 10.0 -c 4 -g 8.0 -e -s 2.0 > $FullDot
		fi
	done
done

mkdir -p $OUT_RANKED
python $BIN_DIR/full_ranker.py $OUTCA $OUT_RANKED

##Window19 models
PG=$BIN_DIR/DAQwindow.py
for mid in {1..20};do #Max 20 models
	daq_file=$OUT_RANKED/rank${mid}.daq
	if [ -e $daq_file ];then
		python $PG $daq_file 9 --OutPath $OUT_RANKED/rank${mid}_
	fi
done

#Check Files....
if "${Buildfullatom_flag}";then
	echo "INFO : DeepMainmast Computation Done with FullAtom Model"
	check_exists $OUT_RANKED/rank1_daq_score_w9.pdb
	cp $OUT_RANKED/rank1_daq_score_w9.pdb $output_dir/DeepMainmast${rand_tag}.pdb
	echo "DONE" >$output_dir/done.out
	exit
fi
