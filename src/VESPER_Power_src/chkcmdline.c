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
	cmd->dreso=16.00;
	cmd->MergeDist=0.50;
	cmd->Filter=0.10;
	cmd->LocalR=10.0;
	cmd->Dkeep=0.5;
	cmd->Nround=5000;
	cmd->Nnb=30;
	cmd->Ntabu=100;
	cmd->Nsim=10;
	cmd->Allow=1.01;
	cmd->Nbeam=20;
	cmd->ssize=2.0;
	cmd->ang=15.0;
	cmd->TopN=10;
	cmd->ShowGrid=false;

	cmd->th1=0.00;
	cmd->th2=0.00;

	cmd->Mode=1;
	cmd->Amode=1;
	cmd->Emode=false;
	cmd->SimpleMode=false;

        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i':
			strcpy(cmd->filename,*(++p));
                	--argc; break;
		 case 'a':
			strcpy(cmd->file1,*(++p));
                	--argc; break;
		 case 'b':
			strcpy(cmd->file2,*(++p));
                	--argc; break;
		 case 'c':
			cmd->Nthr=atoi(*(++p)); 
			--argc; break;
		 case 't':
			cmd->th1=atof(*(++p)); 
			--argc; break;
		 case 'T':
			cmd->th2=atof(*(++p)); 
			--argc; break;
		 case 'g':
			cmd->dreso=atof(*(++p)); 
			--argc; break;
		 case 'M':
                        cmd->Amode=atoi(*(++p));
                        --argc; break;
		 case 'N':
                        cmd->TopN=atoi(*(++p));
                        --argc; break;
		 case 'R':
                        cmd->LocalR=atof(*(++p));
                        --argc; break;
		 case 'k':
                        cmd->Dkeep=atof(*(++p));
                        --argc; break;
		 case 'r':
                        cmd->Nround=atoi(*(++p));
                        --argc; break;
		 //case 'b':
                 //       cmd->Nnb=atoi(*(++p));
                 //       --argc; break;
		 case 'l':
                        cmd->Ntabu=atoi(*(++p));
                        --argc; break;
		 case 's':
                        cmd->ssize=atof(*(++p));
                        --argc; break;
		 //case 'a':
                 //       cmd->Allow=atoi(*(++p));
                 //       --argc; break;
		 //case 'T':
                 //       cmd->Mode=0;
                 //       break;
		 case 'S':
                        cmd->ShowGrid=true;
                        break;
		 case 'V':
                        cmd->Mode=1;
                        break;
		 case 'L':
                        cmd->Mode=2;
                        break;
		 case 'C':
                        cmd->Mode=3;
                        break;
		 case 'P':
                        cmd->Mode=4;
                        break;
		 case 'F':
                        cmd->Mode=5;
                        break;
		 case 'E':
                        cmd->Emode=true;
                        break;
		 case 'e':
                        cmd->Emode=true;
			cmd->SimpleMode=true;
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
	printf("#Map1 Threshold= %f\n",cmd->th1);
	printf("#Map2 Threshold= %f\n",cmd->th2);
	printf("#Band Width= %f\n",cmd->dreso);
	if(cmd->Mode==1) printf("#Vector Mode\n");
	if(cmd->Mode==2) printf("#Overlap Mode\n");
	if(cmd->Mode==3) printf("#CCC Mode\n");
	if(cmd->Mode==4) printf("#Pearson's CCC Mode Normalized by Ave\n");
	if(cmd->Mode==5) printf("#Laplacian Filter Score Mode\n");
	if(cmd->Amode==1) printf("#20.83deg->10.07deg->4.71deg search\n");
	if(cmd->Amode==2) printf("#10.07deg->4.71deg search\n");
	if(cmd->Amode==3) printf("#4.71deg search\n");
        return(TRUE);
}

void errmsg(){
	puts("Usage: VESPOR_Power -a [MAP1.mrc (large)] -b [MAP2.mrc (small)] [(option)]");
	puts("v0.10	Start as EMVEC_FIT_PowerFit");
	puts("v0.20	Pure Multi threading Process & refine 5 degree.");
	puts("v0.30	Add Overlap Mode & using single-precision");
	puts("v0.40	New rotation sampling method");
	puts("v0.41	Fix bugs in Cross Correlation Score");
	puts("v0.42     Add Pearson CC calculation");
	puts("v0.421    Fix problem at Negative density values");
	puts("v0.43     Fix problem at Translation Vector");
	puts("v0.44     Fix problem at ZERO value");
	puts("v0.50	Add scoring mode (-E no rotation, no translation)");
	puts("		Add Laplacian Filter Scoring (-F)");
	puts("---Options---");
	printf("-c [int  ] :Number of cores for threads def=%d\n",2);
	printf("-t [float] :Threshold of density map1 def=%.3f\n",0.00);
	printf("-T [float] :Threshold of density map2 def=%.3f\n",0.00);
	printf("-g [float] : bandwidth of the gaussian filter\n");
        printf("             def=16.0, sigma = 0.5*[float]\n");
	printf("-s [float] : sampling grid space def=2.0\n");
	printf("-M [int  ] : sampling Angle interval Mode 1-3 def=1\n");
	printf("             1: 20.83 degree,   648 samples and refine with 10.07 and 4.71 degree\n");
	printf("             2: 10.07 degree, 7,416 samples and refine with 4.71 degree\n");
	printf("             3: 4.71 degree, 70,728 samples only\n");
	printf("-N [int  ] : Refine Top [int] models def=10\n");
	printf("-S 	   : Show topN models in PDB format def=false\n");
	printf("-V 	   : Vector Products Mode def=true\n");
	printf("-L 	   : Overlap Mode def=false\n");
	printf("-C 	   : Cross Correlation Coefficient Mode def=false\n");
        printf("             Using normalized density Value by Gaussian Filter\n");
	printf("-P         : Pearson Correlation Coefficient Mode def=false\n");
        printf("             Using normalized density Value by Gaussian Filter and average density\n");
	printf("-F         : Laplacian Filtering Mode def=false\n");
	printf("-E         : Evaluation mode of the current position def=false\n");
	printf("-e         : Evaluation mode of the current position without computing Background Distributions def=false\n");
	printf("Thi is Ver %.3f\n",VER);
}
