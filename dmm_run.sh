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
		echo "Missing Value for $2"
		exit 1
	fi
}

show_help() {
    echo "Usage: dmm_fast.sh [option]"
    echo "Fast DeepMaimast with limited parameter combinations"
    echo "Paratemters:"
    echo "  -h, Show help"
    echo "  -p [PATH of the programs]"
    echo "  -m [MAP file]"
    echo "  -f [Sequence File in FASTA format]"
    echo "  -c [Contour Level, density cut-off value]"
    echo "  -o [OutPut Directory]"
    echo "  -t [Time for main-chain tracing (sec)]"
    echo "  -T [Time for fragment assembly (sec)]"
    echo "  -C [Max number of threads]"
}

resume_flag=false
MAX_Threads=8

while getopts "p:m:f:c:o:t:T:C:rh" option; do
	case $option in
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
			echo "Option MAX Thread -C: $OPTARG"
			MAX_Threads=$OPTARG
			;;
		h)
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

echo "PROGRAM_PATH: $PROGRAM_PATH"
echo "map: $map"
echo "fasta: $fasta"
echo "CONTOUR: $CONTOUR"
echo "output_dir: $output_dir"
echo "TIME_TRACE: $TIME_TRACE"
echo "TIME_Assemble: $TIME_ASB"
echo "MAX_Threads: $MAX_Threads"

check_empty "$PROGRAM_PATH" "-p"
check_exists $PROGRAM_PATH/server_bin/MainmastC_UnetAF2

if [ $resume_flag == false ];then
	check_not_exists $output_dir
else
	echo "Resume Process"
fi

check_empty "$CONTOUR" "-c"
check_empty "$output_dir" "-o [output dir]"
check_empty "$TIME_TRACE" "-t [Time for Tracing (sec)]"
check_empty "$TIME_ASB" "-T [Time for Assembly (sec)]"
check_empty "$fasta" "-f [fasta formatfile]"

#af2_model=$5
#exit
#Environment-------------
#PLEASE EDIT THIS SECTION ACCORDING TO 
#YOUR ENVIRONMENT
#source /etc/profile.d/modules.sh
module load miniconda38
eval "$(conda shell.bash hook)"
# Acitivate the conda enviroment
module load deepmainmast-multi
#------------------------
#Programs
BIN_DIR="$PROGRAM_PATH/server_bin/"
EMAP2SF=$BIN_DIR/emap2secplus/main_web.py
NODE=$BIN_DIR/MainmastC_Unet_node


#Set Up files
mkdir -p $output_dir/results/unet

UNET_DIR=$output_dir/results/unet
RESULTS_DIR=$output_dir/results

if [ $resume_flag==false ];then
	cp $map $RESULTS_DIR/input.map
	cp $map $RESULTS_DIR/input.mrc
	cp $fasta $RESULTS_DIR/seq.fasta
fi

SEQ=$RESULTS_DIR/seq.fasta

check_exists $SEQ

#Add Faster resampling Step by Chimera


#Run Emap2sf----------------------------------
if [ ! -e $UNET_DIR/sigmoidAA_GLN.mrc ]||[ ! -e $UNET_DIR/atom_CB.mrc ];then
    python3 $EMAP2SF --gpu 0 --mode=0 -O $RESULTS_DIR -F=$RESULTS_DIR/input.mrc -M=$BIN_DIR/emap2secplus/best_model/Unet_Protein_Atom.pth.tar --contour=$CONTOUR  --batch_size=16 --type=0
    echo "INFO : ATOM prediction Done"
    python3 $EMAP2SF --gpu 0 --mode=0 -O $RESULTS_DIR -F=$RESULTS_DIR/input.mrc -M=$BIN_DIR/emap2secplus/best_model/Unet_Protein_AA.pth.tar --contour=$CONTOUR  --batch_size=16 --type=1
    echo "INFO : AA prediction Done"
fi
#check output files
if [ ! -e $UNET_DIR/sigmoidAA_GLN.mrc ]||[ ! -e $UNET_DIR/atom_CB.mrc ];then
    echo "ERROR : Emap2sf failed"
    exit;
fi


#C-alpha Node Detection-----------------------
for p in 0.3
do
	out=$RESULTS_DIR/NODE_p$p.pdb
	if [ ! -e $out ];then
		hostname > $out
		$NODE -i $UNET_DIR -s $SEQ -t $p -G > $out
	fi
done
echo "INFO : C-alpha Node computing Done"


#Path tracing by VRP solver----------------
ASB=$BIN_DIR/UnetPathEdge.py #Large penalties for edges with Pbb=0
#for p in 0.3 0.4 0.5
for p in 0.3 #x3 faster
do
	for Nchain in 5 10 20 40
	do
		in=$RESULTS_DIR/NODE_p${p}.pdb
		out=$RESULTS_DIR/PATH_p${p}Nch${Nchain}.csv
		log=$RESULTS_DIR/PATH_p${p}Nch${Nchain}.log

		if [ ! -e $log ];then
		    hostname > $log
		    python $ASB $in --Nchain $Nchain --OutCSV $out --SecTrace $TIME_TRACE  &> $log &
		fi
		
	done
done

wait

echo "INFO : Path Tracing by VRP Done"
#exit


PG=$BIN_DIR/MainmastC_UnetAF2
Ncpu=`echo $MAX_Threads/4|bc` #Number of CPU cores
ASB=$BIN_DIR/Assemble_Iter.py
OUTF=$RESULTS_DIR
#Total 4 jobs 5cores = 20 threads
for p in 0.3 #x3 faster
do
	for Nchain in 5 10 20 40
	do

        in=$OUTF/PATH_p${p}Nch${Nchain}.csv
		#UnetDir=results/$ID/unet
		UnetDir=$UNET_DIR

		if [ ! -e $in ]||[ ! -e $SEQ ]||[ ! -e $UnetDir ];then
			echo "Missing....$1 $in"
			continue
		fi


		for Nali in 10 #SpeedUp
		do
			INP=$OUTF/INP_p${p}Nch${Nchain}Nali${Nali}.txt
			OUT=$OUTF/OUT_p${p}Nch${Nchain}Nali${Nali}.pdb
			log=$OUTF/OUT_p${p}Nch${Nchain}Nali${Nali}.log
			if [ ! -e $OUT ];then
			    hostname > $OUT
			    $PG -i $UnetDir -s $SEQ -G -N $in -l 9 -r 1.5 -c $Ncpu -T $Nali > $INP
			    python $ASB --Nstock 10000 --Niter 3 $INP --OutPath $OUT --Ncpu $Ncpu --SecAssemble $TIME_ASB &> $log &
			fi
		done
    done
done

wait

echo "INFO : Assembling Done"


##COMBINATION===================================
echo "INFO : Combine All Models"
PG_ASB=$BIN_DIR/MainmastC_UnetAssembleMtx
PG=$BIN_DIR/Assemble_Iter.py #combine all models
comb_model2=$OUTF/COMBINEi_DMonly.pdb
log=$OUTF/COMBINEi_DMonly.log
modelfile2=$OUTF/MODELs_DMonly.pdb
Ncpu=$MAX_Threads #

if [ ! -e $modelfile2 ];then
	for p in $OUTF/OUT_p*ali10.pdb
	do
		cat $p
		echo "TER"
	done > $modelfile2
fi

mtxfile2=$OUTF/MTX_DMonly.txt

if [ ! -e $mtxfile2 ];then
	$PG_ASB -i $UNET_DIR -s $SEQ -l 10 -A $modelfile2 -G > $mtxfile2
fi

if [ ! -e $comb_model2 ];then
	python3 $PG $mtxfile2 --Nstock 10000 --OutPath $comb_model2 --SecAssemble $TIME_ASB --Ncpu $Ncpu &>$log
fi

echo "INFO : Combine All Models Done"
echo "INFO : DAQ scoring"

DAQ=$BIN_DIR/DAQscore_Unet
DAQ_OUT=$OUTF/COMBINEi_DMonly_daq.pdb

if [ -e $comb_model2 ];then
	$DAQ -i $UnetDir -Q $comb_model2 > $DAQ_OUT
else
    echo "INFO : ERROR"
fi

#WindowSize 9
WINDOW_PG=$BIN_DIR/DAQwindow.py 
python3 $WINDOW_PG $DAQ_OUT 9 --OutPath $OUTF/COMBINEi_DMonly_

echo "INFO : DAQ scoring Done"
echo "INFO : DeepMainmast Done"

#copy job.yml
#cp job.yml $output_dir/.
if [ -e $output_dir/results/daq_score_w9.pdb ]; then
     echo "DONE" >$output_dir/done.out
else
    echo "FAIL" >$output_dir/fail.out
fi
