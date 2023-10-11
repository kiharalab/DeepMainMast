/*
caldep + fragment + sphere
*/
//#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "struct.h"
#include "func.h"
#include "mrc.h"
#include "mrcfft.h"
//#include "scoring.h"

#define PDB_STRLEN 55

void malloc_error(char *a){
 fprintf(stderr,"malloc error in %s\n",a);
 exit(0);
}
double gettimeofday_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int readlist(char *fname,char **list){
 int num=0;
 FILE *fp;
 int len;
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(list[num],LIN,fp)!=NULL){
  len=strlen(list[num]);
  list[num][len-1]='\0';//ignore terminal \n
  num++;
 }
 fclose(fp);
 return TRUE;
}

int line_num(char *fname){
 int num=0;
 FILE *fp;
 char line[LIN];
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(line,LIN,fp)!=NULL){
  num++;
 }
 fclose(fp);
 return num;
}


CMD cmd;

int main(int argc, char **argv)
{
 double t1=gettimeofday_sec();
 double t4;
 POINTS pt;
 MRC mrc,mrc1,mrc2;
 GRAPH g;
 TREE mst;
 //Get Options
 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);

 //Set threads
 if(cmd.Nthr < omp_get_num_procs()){
  omp_set_num_threads(cmd.Nthr);
 }else{
  omp_set_num_threads(omp_get_num_procs());
 }
 

 if(readmrc(&mrc1,cmd.file1))
  return(0);
 if(readmrc(&mrc2,cmd.file2))
  return(0);

//Setup 2^x size 

 MRC mrcN1,mrcN2,mrcSTD1,mrcSTD2;

 SetUpVoxSize(&mrc1,&mrcN1,cmd.th1,cmd.ssize);
 SetUpVoxSize(&mrc1,&mrcSTD1,cmd.th1,1.0);
 SetUpVoxSize(&mrc2,&mrcN2,cmd.th2,cmd.ssize);
 SetUpVoxSize(&mrc2,&mrcSTD2,cmd.th2,1.0);

 if(mrcN1.xdim>mrcN2.xdim){
  mrcN2.xdim=mrcN2.ydim=mrcN2.zdim=mrcN1.xdim;

  mrcN2.orgxyz[0]=mrcN2.cent[0]-0.5*mrcN2.widthx*mrcN2.xdim;
  mrcN2.orgxyz[1]=mrcN2.cent[1]-0.5*mrcN2.widthx*mrcN2.xdim;
  mrcN2.orgxyz[2]=mrcN2.cent[2]-0.5*mrcN2.widthx*mrcN2.xdim;
 }else{
  mrcN1.xdim=mrcN1.ydim=mrcN1.zdim=mrcN2.xdim;
  
  mrcN1.orgxyz[0]=mrcN1.cent[0]-0.5*mrcN1.widthx*mrcN1.xdim;
  mrcN1.orgxyz[1]=mrcN1.cent[1]-0.5*mrcN1.widthx*mrcN1.xdim;
  mrcN1.orgxyz[2]=mrcN1.cent[2]-0.5*mrcN1.widthx*mrcN1.xdim;
 }


 //Common Vox
 MRC mrcCOM1,mrcCOM2;
 if(cmd.Emode==true){
	 SetUpComVoxSize(&mrc1,&mrc2,&mrcCOM1,&mrcCOM2,cmd.th1,cmd.th2,cmd.ssize);
	 puts("#Setting Common boxes");
	 if(fastVEC(&mrc1,&mrcCOM1))
	  return(0);
	 if(fastVEC(&mrc2,&mrcCOM2))
	  return(0);
	if(cmd.SimpleMode==true){
	 if(SimpleScoring(&mrcCOM1,&mrcCOM2))
  	  return(0);
	 return 0;
	}
 }




 //Big
 if(fastVEC(&mrc1,&mrcN1))
  return(0);
 //if(fastVEC(&mrc1,&mrcSTD1))
 // return(0);

 //ShowVec(&mrcN1);
 //Small
 if(fastVEC(&mrc2,&mrcN2))
  return(0);
 if(fastVEC(&mrc2,&mrcSTD2))
  return(0);

 //6D search
 //Vector based
 if(cmd.Mode==1)
  if(SearchMAPfftMT(&mrcN1,&mrcN2,&mrcSTD2,cmd.Amode,cmd.Emode,&mrcCOM1,&mrcCOM2))
   return(0);

 //Overlap based
 if(cmd.Mode==2||cmd.Mode==3||cmd.Mode==4||cmd.Mode==5)
  //if(SearchMAPfftMT_OVCC(&mrcN1,&mrcN2,cmd.ang,cmd.Mode-2))
  if(SearchMAPfftMT_OVCC(&mrcN1,&mrcN2,cmd.Mode-2,cmd.Amode,cmd.Emode))
   return(0);

 t4=gettimeofday_sec();
 printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
 return 0;

}

