# DeepMainMast

<a href="https://github.com/marktext/marktext/releases/latest">
   <img src="https://img.shields.io/badge/DeepMainMast-v1.0.0-green">
   <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Mac%20-green">
   <img src="https://img.shields.io/badge/Language-python3-green">
   <img src="https://img.shields.io/badge/dependencies-tested-green">
   <img src="https://img.shields.io/badge/licence-GNU-green">
</a>  

<img src="DeepMainMast_Logo.png" height=380 width=300/>

DeepMainMast is a computational tool using deep learning to automatically build full protein complex structure from cryo-EM map.  

Copyright (C) 2022 Genki Terashi, Xiao Wang, Devashish Prasad, Tsukasa Nakamura, Daisuke Kihara, and Purdue University. 

License: GPL v3. (If you are interested in a different license, for example, for commercial use, please contact us.) 

Contact: Daisuke Kihara (dkihara@purdue.edu)

For technical problems or questions, please reach to Genki Terashi (gterashi@purdue.edu).

## Citation:

Genki Terashi, Xiao Wang, Devashish Prasad, Tsukasa Nakamura & Daisuke Kihara. Nature Methods, 2023.
[Paper]()
```
@article{Terashi2023deepmainmast,   
  title={Integrated Protocol of Protein Structure Modeling for cryo-EM with Deep Learning and Structure Prediction},   
  author={Genki Terashi, Xiao Wang, Devashish Prasad, Tsukasa Nakamura, and Daisuke Kihara},    
  journal={Nature Methods},    
  year={2023}    
}   
```

# Online Platform:

## Server(Recommended): https://github.com/kiharalab/DeepMainMast
<details>
We have four publicly available platforms, which basically offer similar functionality.
Input: cryo-EM map+sequence file. Output: modeled protein structure. The input and output are the same across all platforms.
</details>

### Google Colab: https://colab.research.google.com/github/kiharalab/DeepMainMast/blob/main/DeepMainMast_Multi_chain.ipynb
<details> 
   Step-by-step instructions are available. For free user, colab has 4-hour running time limit and may not work for large structures (>=1000 residues).
</details>

### Code Ocean: https://codeocean.com/capsule/0749800
<details> 
   Free online platform for easy usage. For academic users, CodeOcean has 10-hour running time limit per month.
</details>

### Local installation with source code at Github
<details>
Full code is available here and it is easier for user to modify to develop their own tools.
</details>

### Project website: https://kiharalab.org/emsuites
### Detailed pipeline instructions can be found https://kiharalab.org/emsuites/deepmainmast.php

## Introduction
<details>
   <summary>DeepMainMast is a computational tool using deep learning to automatically build full protein complex structure from cryo-EM map. </summary>
Structure modeling from maps is an indispensable step for studying proteins and their complexes with cryogenic electron microscopy (cryo-EM). Although the resolution of determined cryo-EM maps has generally improved, there are still many cases where tracing protein main-chains is difficult, even in maps determined at a near atomic resolution. Here, we have developed a protein structure modeling method, called DeepMainmast, which employs deep learning to capture the local map features of amino acids and atoms to assist main-chain tracing. Moreover, since Alphafold2 demonstrates high accuracy in protein structure prediction, we have integrated complementary strengths of de novo density tracing using deep learning with Alphafold2’s structure modeling to achieve even higher accuracy than each method alone. Additionally, the protocol is able to accurately assign chain identity to the structure models of homo-multimers.![image](https://github.com/kiharalab/dmm/assets/50850224/f04a5a46-7cd0-4961-a1c4-e7e044ea8870)

</details>

## Overall Protocol 
<details>
(1) Detecting amino-acid types and atom types using deep learning (Emap2sf). The image on the left shows the detected atom types (Ca atom: green, carbon: orange, and nitrogen: light blue). The image on the right shows the detected amino acid types in different colors. <br>
(2) Tracing Ca path and assigning the target sequence using the Vehicle Routing Problem Solver and the Dynamic Programming algorithm. Different parameter combinations are used. <br>
(3) Assembling Ca fragments using the Constraint Problem (CP) Solver. Colors indicate chain IDs. <br>
(4) Combining Ca models built under different parameter combinations using the CP Solver. Colors indicate the direction of chains from blue to red for the N-terminal to the C-terminal residues. <br>
(5) Full-atom building and refinement using PULCHRA and Rosetta-CM. <br>
(6) Scoring generated full-atom models based on the DAQ(AA) score and the DOT score.<br>
 
<p align="center">
  <img src="dmm_pipeline.png" alt="DeepMainMast framework" width="70%">
</p>
</details>

## Installation
<details>

### System Requirements
CPU: >=8 cores <br>
Memory (RAM): >=50Gb. For maps with more than 3,000 residues, memory space should be higher than 200GB. <br>
GPU: any GPU supports CUDA with at least 12GB memory. <br>
GPU is required for DeepMainMast and no CPU version is available for CryoREAD since it is too slow.

## Installation Instructions
### 1. [`Install git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) 
### 2. Clone the repository in your computer 
```
git clone  https://github.com/kiharalab/DeepMainMast.git && cd DeepMainMast
```

### 3. Build dependencies.   
You have two options to install dependency on your computer:
#### 3.2 Install with anaconda (Recommended)
##### 3.2.1 [`install anaconda`](https://www.anaconda.com/download). 
##### 3.2.2 Install dependency in the command line
Make sure you are in the DeepMainMast directory and then run 
```
conda env create -f environment.yml
```
Each time when you want to run this software, simply activate the environment by
```
conda activate deepmainmast
conda deactivate(when you want to exit) 
```

### 4. Compile C packages
Run the following command:
```
bash make_c_programs.sh
```


</details>

## Examples and commands
In data/, there are three examples (3j9sA, 1461 and 2513).
### Single-Chain Modeling (3j9sA)
#### Input files in data/3j9sA
+ Target sequence file (-f): 3j9sA.fasta
+ MAP file (-m): 3j9sA.mrc
+ AlphaFold2 Model (optional -A): 3j9sA_af2.pdb
#### Parameters (example)
+ Maximum number of CPU cores to use (-C): 20
+ Maximum number of CPU cores per one thread (-M): 5
+ Density contour level (-c): 0.01
+ Computational Time Limit for PATH tracing per thread (-t): 600 sec = 10 mins
+ Computational Time Limit for Fragment assembly (-T): 600 sec = 10 mins
#### Command for Calpha PATH tracing
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.01 -o [OUTPUT PATH] -t 600 -T 600 -C 20 -M 5 -m ./data/3j9sA/3j9sA.mrc -f ./data/3j9sA/3j9sA.fasta
```
#### Command for Calpha PATH tracing using AlphaFold2 Model
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.01 -o [OUTPUT PATH] -t 600 -T 600 -C 20 -M 5 -m ./data/3j9sA/3j9sA.mrc -f ./data/3j9sA/3j9sA.fasta -A ./data/3j9sA_af2.pdb
```
#### Command for Calpha PATH tracing using AlphaFold2 Model and Full-Atom Model building&Refinement by Rosetta
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.01 -o [OUTPUT PATH] -t 600 -T 600 -C 20 -M 5 -m ./data/3j9sA/3j9sA.mrc -f ./data/3j9sA/3j9sA.fasta -A ./data/3j9sA_af2.pdb -x [ROSETTA PROGRAM PATH]
```

### Multi-Chain Modeling (EMD-1461, Homooligomer)
#### Input files in data/1461
+ Target sequence file (-f): emd_1461.fasta
+ MAP file (-m): emd_1461.mrc
+ AlphaFold2 Model (optional -A): emd_1461_af2.pdb
#### Parameters (example)
+ Maximum number of CPU cores to use (-C): 18
+ Maximum number of CPU cores per one thread (-M): 6
+ Density contour level (-c): 0.3
+ Computational Time Limit for PATH tracing per thread (-t): 1200 sec = 20 mins
+ Computational Time Limit for Fragment assembly (-T): 600 sec = 10 mins
+ Including Homo-oliogomer: -H option
#### Command for Calpha PATH tracing
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.3 -o [OUTPUT PATH] -t 1200 -T 600 -C 18 -M 6 -m ./data/1461/emd_1461.mrc -f ./data/1461/emd_1461.fasta  -H -A ./1461/emd_1461_af2.pdb -x [ROSETTA PROGRAM PATH]
```
#### Command for Calpha PATH tracing using AlphaFold2 Model
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.3 -o [OUTPUT PATH] -t 1200 -T 600 -C 18 -M 6 -m ./data/1461/emd_1461.mrc -f ./data/1461/emd_1461.fasta  -H -A ./1461/emd_1461_af2.pdb 
```
#### Command for Calpha PATH tracing using AlphaFold2 Model and Full-Atom Model building&Refinement by Rosetta
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.3 -o [OUTPUT PATH] -t 1200 -T 600 -C 18 -M 6 -m ./data/1461/emd_1461.mrc -f ./data/1461/emd_1461.fasta  -H -A ./1461/emd_1461_af2.pdb -x [ROSETTA PROGRAM PATH]
```

### Multi-Chain Modeling (EMD-2513, HeteroOligomer)
#### Input files in data/2513
+ Target sequence file (-f): emd_2513.fasta
+ MAP file (-m): emd_2513.mrc
+ AlphaFold2 Model (optional -A): emd_2513_af2.pdb
#### Parameters (example)
+ Maximum number of CPU cores to use (-C): 18
+ Maximum number of CPU cores per one thread (-M): 6
+ Density contour level (-c): 0.01
+ Computational Time Limit for PATH tracing per thread (-t): 1200 sec = 20 mins
+ Computational Time Limit for Fragment assembly (-T): 600 sec = 10 mins
#### Command for Calpha PATH tracing
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.01 -o [OUTPUT PATH] -t 1200 -T 600 -C 18 -M 6 -m ./data/2513/emd_2513.mrc -f ./data/2513/emd_2513.fasta  -A ./2513/emd_2513_af2.pdb -x [ROSETTA PROGRAM PATH]
```
#### Command for Calpha PATH tracing using AlphaFold2 Model
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.01 -o [OUTPUT PATH] -t 1200 -T 600 -C 18 -M 6 -m ./data/2513/emd_2513.mrc -f ./data/2513/emd_2513.fasta  -A ./2513/emd_2513_af2.pdb 
```
#### Command for Calpha PATH tracing using AlphaFold2 Model and Full-Atom Model building&Refinement by Rosetta
```
./dmm_full_multithreads.sh -p [PROGRAM PATH (./)] -c 0.01 -o [OUTPUT PATH] -t 1200 -T 600 -C 18 -M 6 -m ./data/2513/emd_2513.mrc -f ./data/2513/emd_2513.fasta  -A ./2513/emd_2513_af2.pdb -x [ROSETTA PROGRAM PATH]
```
