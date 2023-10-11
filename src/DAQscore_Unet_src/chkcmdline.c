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
		 case 'Q':
			cmd->Amode=true;
			strcpy(cmd->Afilename,*(++p));
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
                        cmd->LocalR=atof(*(++p));
                        --argc; break;
		 case 'k':
                        cmd->Dkeep=atof(*(++p));
                        --argc; break;
		 case 'w':
                        cmd->Wd=atof(*(++p));
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
		 //case 's':
                 //       cmd->Nsim=atoi(*(++p));
                 //       --argc; break;
		 case 'a':
                        cmd->Allow=atoi(*(++p));
                        --argc; break;
		 
		 case 'V': //Visualize mode
                        cmd->Mode=4;
                        break;
		 case 'T':
                        cmd->Mode=0;
                        break;
		 case 'G':
                        cmd->Mode=1;
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

        return(TRUE);
}

void errmsg(){
	puts("Usage: DAQscore_Unet -i [U-net output dir] -Q [PDB file][(option)]");
	puts("v1.0	Using Unet");
	puts("---Options---");
	printf("-c [int  ] :Number of threads def=%d\n",2);
	printf("-t [float] :Threshold of the Probability map def=%.3f\n",0.00);
	printf("-g [float] : bandwidth of the gaussian filter\n");
        printf("             def=2.0, sigma = 0.5*[float]\n");
	printf("-f [float] :Filtering for representative points def=%.3f\n",0.10);
	printf("-m [float] :After MeanShifting merge<[float] def=%.3f\n",1.00);
	printf("-R [float] :Radius of Local MST def=%.3f\n",10.0);
	printf("-k [float] :keep edges where d<[float] def=%.3f\n",0.5);
	//printf("-w [float] :Weight of Probability(N,CA,C) def=%.3f\n",1.0);
	printf("-w [float] :Weight of Density def=%.3f\n",1.0);
	printf("-W [float] :Weight of Calpha Prob def=%.3f\n",2.0);
	printf("-P [float] :Probability(N,CA,C)>P_CutOff def=%.3f\n",0.5);
	printf("-l [int  ] :Sequence Fragment Length def=%d\n",9);
	printf("-z [float] :z-score cutoff of Sequence Fragment def=%f\n",1.0);
	//printf("-A [file ] :Another PDB model (Assemble model)\n");
	printf("-r [float] :RMSD cutoff def=2.0(Only for AlphaFold2 model)\n");
	printf("-D         :Dummy fragment mode\n");
	
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
