# deepmainmast_multithreads (Full pipleline of DeepMainmast using multiple CPU cores)

#### 0. Compile codes in dmm_fast
- Please compile all codes in dmm_fast/ by
<pre>
bash ./make_c_program.sh
</pre>
- If you use a specific conda environment, please specify it at L212-219 in dmm_full_multithreads.sh.
<pre>
#Environment-------------
#PLEASE EDIT THIS SECTION ACCORDING TO 
#YOUR ENVIRONMENT
module load miniconda38
eval "$(conda shell.bash hook)"
# Acitivate the conda enviroment
module load deepmainmast-multi
#------------------------
</pre>

#### 1. What is deepmainmast_multithreads?
- deepmainmast_multithreads performs full pipeline of the DeepMainmast.
- It can use AlphaFold2 model to enhance the performance.
- If specified, It can use full-atom refinement by RosettaCM.
- It can use multiple CPU cores at the same time.

#### 2. Configure parameters for running DeepMainMast_multithreads

You can control the execution of DeepMainMast_multithreads pipeline using in dmm_full_multithreads.sh

Refer the following details to set each parameter under Configuration parameters section. 

- -p [PATH] : Specify the path of 'dmm_fast/'. <br>
If you pulled DeepMainmast in /home/XXXX by  <pre>git pull https://github.com/kiharalab/DeepMainMast/ </pre> you can specify the path as <pre>-p /home/XXXX/DeepMainmast/dmm_fast/</pre>

- -m [MAP file] : the path of the EM map file.

- -t [Tracing Time (sec)] : You can specify the number of seconds to be used per iteration for the main-chain tracing steps using -t [sec] parameter. It is recommended to start with 1800 for your new map of multi-chain structure. 

- -T [Fragment Assemble Time] : You can specify the number of seconds to be used per iteration for the fragment assembly steps using -T [sec] parameter. It is recommended to start with 600 for your new map of multi-chain structure.

- -c [Contour Level] : Specify the recommended contour level for your map.

- -f [FASTA file] : Specify the target sequence file. Please use the standard 20 Amini Acid Typoes. Please remove "X" or other characters which are not included in the 20 Amini Acid Types. The format should be FASTA format.
<pre>
>ChainA
AAAAAAAAAA
>ChainB
BBBBBBBBBB
</pre>
- -C [Number of threads] : Specify the maximum number of threads. For main-chain path tracing, the program uses this number of threads.
- -M [Maximum number of CPU cores per one thread] : For the fragment-assembly step, the program uses this number of CPU cores per one thread. For example, -C 24 -M 6 executes 4 threads with 6 CPU cores for each.

- -o [Output directory path] : Specify the output directory. The command makes the specified directory.

- -r : resume previous jobs in the Output directory.

- -x [Rosetta proram PATH] : (example: /home/username/rosetta_bin_linux_XXX_bundle/) If this option is not used, the program generate only CA-atom models.

- -A [AlphaFold2 Model]: Please renumber the residue numbers in the model based on the target sequence.

- -H : Homo-oligomer mode. If specified, the program performs the chain-ID reassignment.

#### 3. Run DeepMainMast_fast
Run the multi-chain pipeline using the following command

<pre>
./dmm_full_multithreads.sh -p /home/XXXX/DeepMainmast/dmm_fast -m emd_3073.mrc -f input.fasta -C 12 -M 4 -c 0.1 -o ./OutPutDir -A AF2model.pdb -x /home/XXXX/rosetta_bin_linux_XXXX/
</pre>

#### 4. Get Output of DeepMainMast 
Your results will be generated under the output directory (specified by -o [PATH]).  

#### 5. Output Files
The DeepMainmast protocol generates the following output files in the output directory.

#### Examples of output files (https://kiharalab.org/emsuites/deepmainmast/examples/)
- Multimer target EMD-3073 (https://kiharalab.org/emsuites/deepmainmast/examples/3073_dmmfull.zip)
##### Emap2sf output (Predicting local propaties of the protein structure): unet/
###### Atom type prediction:
<pre>
atom_BG.mrc		Background
atom_C.mrc		Backbone Carbon
atom_CA.mrc		Calpha
atom_CB.mrc		Cbeta
atom_N.mrc		Backbone Nitrogen
atom_O.mrc		Backbone Oxgen
atom_Others.mrc	Other atoms, Side-chain Atoms
</pre>

###### Amino Acid Type Prediction: unet/
<pre>
sigmoidAA_XXX.mrc 20 amino acid type
</pre>
##### Predicted Local Dense Points:
<pre>
NODE_p0.3.pdb	NODE_p0.4.pdb	NODE_p0.5.pdb
</pre>
##### Computed Paths using VRP Solver:
<pre>
PATH_p*Nch*.csv
</pre>
##### Computed Fragment Library:
<pre>
INP_p*Nch*Nali*.txt
</pre>
##### Computed Fragment Library with the AF2 model (if provided):
<pre>
INP_p*Nch*Nali*R*.txt
</pre>
##### Assembled fragments in PDB format:
<pre>
For each fragment library (INP*.txt), DeepMainmast generates one output (OUT*.pdb)
OUT_p*Nch*Nali*.pdb
</pre>

##### Assembled fragments in PDB format:
<pre>
For each fragment library (INP*.txt), DeepMainmast generates one output (OUT*.pdb)
OUT_p*Nch*Nali*.pdb
</pre>
##### Map-model fitting results using VESPER (if AF2 model is provided): VESPER_MODELs/
<pre>
af2_A_R*.out			output file of VESPER computation
R*_A_FIT_MODEL*.pdb		Fitted models
</pre>

##### Input files for Assembling C-alpha Models:
###### Input1: concatenated models
<pre>
MODELs_AFonly.pdb		Concatenated all OUT*.pdb files
MODELs_DMonly.pdb		Concatenated OUT*.pdb files without the AF2 data
MODELs_VESPER.pdb		Concatenated fitted models in VESPER_MODELs/
MODELs_all.pdb			Concatenated all OUT*.pdb files
</pre>

###### Input2: matrix files
<pre>
MTX_AFonly.txt
MTX_DMonly.txt
MTX_VESPER.txt
MTX_all.txt
</pre>
##### Output: Assembled Calpha Models
<pre>
COMBINEi_AFonly.pdb
COMBINEi_DMonly.pdb
COMBINEi_VESPER.pdb
COMBINEi_all.pdb
</pre>
###### For Homo-oligomer target, DeepMainmast refines the chain ID assignment (If specified):
<pre>
COMBINEi_AFonly_rechain.pdb
COMBINEi_DMonly_rechain.pdb
COMBINEi_VESPER_rechain.pdb
COMBINEi_all_rechain.pdb
</pre>
##### Ranked C-alpha Models using DAQ and DOT scores:
<pre>
FINAL_CA_MODEL/
rank1.pdb
rank2.pdb
rank3.pdb
rank4.pdb
</pre>
##### Full-atom modeling results:
<pre>
FINAL_CA_MODELs
CM_AFonly/
CM_DMonly/
CM_VESPER/
CM_all/
</pre>
##### Ranked Full-atom Models using DAQ and DOT score:
<pre>
RANKDED_DATA/
#Full-atom models
rank1.pdb
rank2.pdb
rank3.pdb
rank4.pdb
...
#full-atom models with DAQ score
rank1_daq_score_w9.pdb
rank2_daq_score_w9.pdb
rank3_daq_score_w9.pdb
rank4_daq_score_w9.pdb
...
#Ranking table
ranking.csv
</pre>

##### Ranked Full-atom Models with DAQ(AA) scores:
rank{N}_daq_score_w9.pdb (N=1,2,..) contains two MODELs (MODEL1 and MODEL2). In MODEL1, b-factor values represent DAQ(AA) score with 19 residues sliding window. In MODEL2, amino acids with DAQ(CA) score below -0.5 are excluded, and amino acids at locations with DAQ(AA) score below -0.5 are replaced with "UNK."





# deepmainmast_fast (simplified version)
### Running DeepMainMast on your EM map
#### 1. How is the deepmainmast_fast faster than original version?

- deepmainmast_fast uses only one cutoff value 0.3 of Ca-atom probability. (Original version uses 0.3, 0.4 and 0.5)
- deepmainmast_fast uses 4 parameters in the main-chain path tracing. (Original version uses 25 parameters)
- deepmainmast_fast computes 10 alignments for each main-chain path. (Original version computes 5 and 10 alignments separately)
- deepmainmast_fast uses limited number of parameter combinations (1/36 of the original)
- deepmainmast_fast does not perform Full-atom model refinement.


#### 2. Configure parameters for running DeepMainMast_fast

You can control the execution of DeepMainMast_fast pipeline using 7 parameters in dmm_run.sh

Refer the following details to set each parameter under Configuration parameters section. 

- -p [PATH] : Specify the path of 'dmm_fast/'. <br>
If you pulled DeepMainmast in /home/XXXX by  <pre>git pull https://github.com/kiharalab/DeepMainMast/ </pre> you can specify the path as <pre>-p /home/XXXX/DeepMainmast/dmm_fast/</pre>

- -m [MAP file] : the path of the EM map file.

- -t [Tracing Time (sec)] : You can specify the number of seconds to be used per iteration for the main-chain tracing steps using -t [sec] parameter. It is recommended to start with 1800 for your new map of multi-chain structure. 

- -T [Fragment Assemble Time] : You can specify the number of seconds to be used per iteration for the fragment assembly steps using -T [sec] parameter. It is recommended to start with 600 for your new map of multi-chain structure.

- -c [Contour Level] : Specify the recommended contour level for your map.

- -f [FASTA file] : Specify the target sequence file. The format should be FASTA format, like 
<pre>
>ChainA
AAAAAAAAAA
>ChainB
BBBBBBBBBB
</pre>

- -C [Number of threads] : Specify the maximum number of threads.

- -o [Output directory path] : Specify the output directory. The command makes the specified directory.



#### 3. Run DeepMainMast_fast
Run the multi-chain pipeline using the following command

<pre>
./dmm_run.sh -p /home/XXXX/DeepMainmast/dmm_fast -m emd_3073.mrc -f input.fasta -C 12 -c 0.1 -o ./OutPutDir
</pre>

#### 4. Get Output of DeepMainMast 
Your results will be generated under the output directory (specified by -o [PATH]).  

#### 5. Output Files
The DeepMainmast protocol generates the following output files in the output directory.

#### Examples of output files (https://kiharalab.org/emsuites/deepmainmast/examples/)
- Multimer target EMD-3073 (https://kiharalab.org/emsuites/deepmainmast/examples/3073_dmmfast.zip)
##### Emap2sf output (Predicting local propaties of the protein structure): unet/
###### Atom type prediction:
<pre>
atom_BG.mrc		Background
atom_C.mrc		Backbone Carbon
atom_CA.mrc		Calpha
atom_CB.mrc		Cbeta
atom_N.mrc		Backbone Nitrogen
atom_O.mrc		Backbone Oxgen
atom_Others.mrc	Other atoms, Side-chain Atoms
</pre>

###### Amino Acid Type Prediction: unet/
<pre>
sigmoidAA_XXX.mrc 20 amino acid type
</pre>
##### Predicted Local Dense Points:
<pre>
NODE_p0.3.pdb	NODE_p0.4.pdb	NODE_p0.5.pdb
</pre>
##### Computed Paths using VRP Solver:
<pre>
PATH_p*Nch*.csv
</pre>
##### Computed Fragment Library:
<pre>
INP_p*Nch*Nali*.txt
</pre>
##### Computed Fragment Library with the AF2 model (if provided):
<pre>
INP_p*Nch*Nali*R*.txt
</pre>
##### Assembled fragments in PDB format:
<pre>
For each fragment library (INP*.txt), DeepMainmast generates one output (OUT*.pdb)
OUT_p*Nch*Nali*.pdb
</pre>


##### Input files for Assembling C-alpha Models:
###### Input1: concatenated models
<pre>
MODELs_DMonly.pdb		Concatenated OUT*.pdb files without the AF2 data
</pre>

###### Input2: matrix files
<pre>
MTX_DMonly.txt
</pre>
##### Output: Assembled Calpha Models

<pre>
COMBINEi_DMonly.pdb : CA model
COMBINEi_DMonly_daq_score_w9.pdb : CA model with DAQ(AA) score.
</pre>
COMBINEi_DMonly_daq_score_w9.pdb contains two MODELs (MODEL1 and MODEL2). In MODEL1, b-factor values represent DAQ(AA) score with 19 residues sliding window. In MODEL2, amino acids with DAQ(CA) score below -0.5 are excluded, and amino acids at locations with DAQ(AA) score below -0.5 are replaced with "UNK."

## Cite this work 

>"Integrated Protocol of Protein Structure Modeling for cryo-EM with Deep Learning and Structure Prediction, Genki Terashi, Xiao Wang, Devashish Prasad, and Daisuke Kihara, In submission (2023)"