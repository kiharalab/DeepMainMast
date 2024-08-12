#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "struct.h"


int chkcmdline( int argc, char **argv,CMD *cmd){
        char **p;
	int num=0;
	void errmsg();
        if(argc < 2){
		errmsg();
                return(FALSE);
        }
	//default values
        p=argv;
	cmd->map_t=0.00;
	cmd->Nthr=2;
	cmd->dreso=2.00;
	//cmd->MergeDist=0.50;
	cmd->MergeDist=1.00;//2021
	cmd->Filter=0.10;
	cmd->Mode=0;
	cmd->LocalR=10.0;
	cmd->Dkeep=0.5;
	cmd->Nround=5000;
	cmd->Nnb=30;
	cmd->Ntabu=100;
	cmd->Nsim=10;
	cmd->Allow=1.01;
	cmd->Pcut=0.50;
	cmd->Nbeam=20;
	cmd->Wp=1.0;
	cmd->Wd=1.0;
	cmd->Wca=2.0;
	cmd->FragLen = 9;
	cmd->zcut = 1.0;
	cmd->Amode=false;
	cmd->RMSD=2.0;
	cmd->Niter=20;
	cmd->MinReg=20;
	cmd->Weight=10.0;
        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i':
			strcpy(cmd->dirname,*(++p));
                	--argc; break;
		 case 'p':
			strcpy(cmd->pfilename,*(++p));
                	--argc; break;
		 case 's':
			strcpy(cmd->sfilename,*(++p));
                	--argc; break;
		 case 'A':
			cmd->Amode=true;
			strcpy(cmd->Afilename,*(++p));
                	--argc; break;
		 case 'N':
			//cmd->Amode=true;
			strcpy(cmd->Nfilename,*(++p));
			printf("#PathFile: %s\n",cmd->Nfilename);
                	--argc; break;
		 case 'D':
			cmd->DummyMode=true;
                	break;
		 case 'r':
			cmd->RMSD=atof(*(++p)); 
                	--argc; break;
		 case 'c':
			cmd->Nthr=atoi(*(++p)); 
			--argc; break;
		 case 't':
			cmd->map_t=atof(*(++p)); 
			--argc; break;
		 case 'g':
			cmd->dreso=atof(*(++p)); 
			--argc; break;
		 case 'f':
                        cmd->Filter=atof(*(++p));
                        --argc; break;
		 case 'm':
                        cmd->MergeDist=atof(*(++p));
                        --argc; break;
		 case 'R':
                        cmd->MinReg=atoi(*(++p));
                        --argc; break;
		 case 'k':
                        cmd->Dkeep=atof(*(++p));
                        --argc; break;
		 case 'w':
                        cmd->Weight=atof(*(++p));
                        --argc; break;
		 case 'W':
                        cmd->Wca=atof(*(++p));
                        --argc; break;
		 case 'P':
                        cmd->Pcut=atof(*(++p));
                        --argc; break;
		 //case 'r':
                 //       cmd->Nround=atoi(*(++p));
                 //       --argc; break;
		 case 'b':
                        cmd->Nnb=atoi(*(++p));
                        --argc; break;
		 case 'l':
                        cmd->FragLen=atoi(*(++p));
                        --argc; break;
		 case 'z':
                        cmd->zcut=atof(*(++p));
                        --argc; break;
		 case 'T':
                        cmd->Niter=atoi(*(++p));
                        --argc; break;
		 //case 's':
                 //       cmd->Nsim=atoi(*(++p));
                 //       --argc; break;
		 case 'a':
                        cmd->Allow=atoi(*(++p));
                        --argc; break;
		 
		 case 'V': //Visualize mode
                        cmd->Mode=4;
                        break;
		 //case 'T':
                 //       cmd->Mode=0;
                 //       break;
		 case 'G':
                        cmd->DummyMode=true;
                        break;
		 case 'M':
                        cmd->Mode=2;
                        break;
		 case 'L':
                        cmd->Mode=3;
                        break;
		 default: 
		  	fprintf(stderr,"No such a option: %s\n",*p+1); 
			errmsg(); 
			return(FALSE); 
		 break;
	  }
	 }
        }
	//option check
	printf("#Number of threads= %d\n",cmd->Nthr);
	printf("#Map Threshold= %f\n",cmd->map_t);
	printf("#Band Width= %f\n",cmd->dreso);
	printf("#Merge Dist= %f\n",cmd->MergeDist);
	printf("#Filtering Dens Cut= %f\n",cmd->Filter);
	printf("#Mode= %d 0:Tabu, 1:Graph, 2: MST\n",cmd->Mode);
	printf("#Backbone Pcut= %f\n",cmd->Pcut);
	//printf("#Weight of Backbone P= %f\n",cmd->Wp);
	printf("#Weight of Density = %f\n",cmd->Wd);
	printf("#Weight of Calpha Prob = %f\n",cmd->Wca);
	printf("#Weight of BB_N Prob = %f\n",1.00);
	printf("#Weight of BB_C Prob = %f\n",1.00);
	printf("#PTABLE: %s\n",cmd->pfilename);
	printf("#SEQ: %s\n",cmd->sfilename);
	printf("#Length of Fragment: %d\n",cmd->FragLen);
	if(cmd->Amode==true){
	 printf("#Amodel mode: %s\n",cmd->Afilename);
	 printf("#RMSD cutoff: %f\n",cmd->RMSD);
	 printf("#MinimumSizeOfRegions: %d\n",cmd->MinReg);
	}
	if(cmd->DummyMode==true){
	 puts("#Dummy fragment mode: On");
	}

        return(TRUE);
}

void errmsg(){
	puts("Usage: MainmastC_UnetAF2 -i [U-net output dir] -s [FASTA format file] -A [AF2 pdb file] -N [PATH.csv file] [(option)]");
	puts("Fast Mainmast C-lang&multi thread version");
	puts("v1.01	Add Segmentation part by using MST");
	puts("v1.02	Bug fix in malloc and Tabu-search list");
	puts("v1.03	origina1 Bug fix Ncid");
	puts("v1.031	Add Fast LDP mode");
	puts("v2.0	Using Back-bone Prob");
	puts("v3.0	Bug fix in Density weight");
	puts("v3.1	Bug fix in Extend mode");
	puts("v4.0	Path & AF model");
	puts("v4.1	Bug fix in memory leak");
	puts("---Options---");
	printf("-c [int  ] :Number of threads def=%d\n",2);
	printf("-t [float] :Threshold of the Probability map def=%.3f\n",0.00);
	printf("-g [float] : bandwidth of the gaussian filter\n");
        printf("             def=2.0, sigma = 0.5*[float]\n");
	printf("-f [float] :Filtering for representative points def=%.3f\n",0.10);
	//printf("-m [float] :After MeanShifting merge<[float] def=%.3f\n",1.00);
	//printf("-R [float] :Radius of Local MST def=%.3f\n",10.0);
	//printf("-k [float] :keep edges where d<[float] def=%.3f\n",0.5);
	//printf("-w [float] :Weight of Probability(N,CA,C) def=%.3f\n",1.0);
	//printf("-w [float] :Weight of Density def=%.3f\n",1.0);
	//printf("-W [float] :Weight of Calpha Prob def=%.3f\n",2.0);
	printf("-P [float] :Probability(N,CA,C)>P_CutOff def=%.3f\n",0.5);
	printf("-l [int  ] :Sequence Fragment Length def=%d\n",9);
	//printf("-z [float] :z-score cutoff of Sequence Fragment def=%f\n",1.0);
	//printf("-A [file ] :Another PDB model (Assemble model)\n");
	printf("-r [float] :RMSD cutoff def=2.0(Only for AlphaFold2 model)\n");
	printf("-R [int]   :Minimum Aligned position in AF2 mode def=20\n");
	printf("-D         :Dummy fragment mode\n");
	printf("-G         :SigmoidAA mode\n");
	printf("-T [int]   :Number of Alignments per path def=20\n");
	printf("-w [float] :Score weight for AF2 model def=10.0\n");
/*
	puts("---Options for Tabu Search---");
	printf("-r [int  ] :Number of rounds def=%d\n",5000);
	printf("-b [int  ] :Number of neighbors def=%d\n",30);
	printf("-l [int  ] :Size of Tabu-list def=%d\n",100);
	printf("-s [int  ] :Number of Models def=%d\n",10);
	printf("-a [float] :Total(Tree) < Total(MST)*[float] def=%f\n",1.01);
*/
	printf("Thi is Ver %.3f\n",VER);
}
