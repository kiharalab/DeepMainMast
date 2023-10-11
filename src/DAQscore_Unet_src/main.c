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
#include "thread.h"
#include "rmsd.h"
#include <sys/stat.h>


//#include "scoring.h"

#define PDB_STRLEN 55

//Fragment

typedef struct{
 NODE Ca[200];
 int Start,End,Len;
 float RawScore,Score,Rmsd;
} RESULT_MODEL;

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


/*
int AddConnectModel(RESULT_MODEL **, int,int,SEQ *);
int CheckOverlap(int,int,int,int);
bool CheckConf(RESULT_MODEL *,RESULT_MODEL *);
float Dist2(NODE *,NODE *);

void TableShow(GRAPH *, SEQ *,MRC *,float);
float CheckConnect(NODE *,NODE *,MRC *,float);
*/

bool ShowDAQscore(GRAPH *, PDB *);

CMD cmd;

int main(int argc, char **argv)
{
 double t1=gettimeofday_sec();
 double t4;
 POINTS pt;
 MRC mrc,pmrc,ori_mrc,bb_mrc;
 GRAPH g;
 TREE mst;
 PTBL ptbl;
 //Get Options
 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);

 //Set threads
 if(cmd.Nthr < omp_get_num_procs()){
  omp_set_num_threads(cmd.Nthr);
 }else{
  omp_set_num_threads(omp_get_num_procs());
 }


 //Loading U-net outputs=============================
 //AAs
 MRC *AA_mrc,*ATOM_mrc;

 AA_mrc=(MRC *)malloc(sizeof(MRC)*20);
 ATOM_mrc=(MRC *)malloc(sizeof(MRC)*6);

 char mapfile[1000];
 struct stat buf;
 for(int aa=0;aa<20;aa++){
	//sprintf(mapfile,"%s/%s%s.mrc",cmd.dirname,"AA_",RES_NAMES[aa]);
	sprintf(mapfile,"%s/sigmoid%s%s.mrc",cmd.dirname,"AA_",RES_NAMES[aa]);
	printf("##Loading MAP: %s\n",mapfile);
	if(stat(mapfile,&buf)){ printf("##ERROR MAP: %s\n",mapfile); return 0; }
  if(readmrc(&AA_mrc[aa],mapfile))
   return(0);
 }

 //ATOMs========
 //Others
 sprintf(mapfile,"%s/atom_Others.mrc",cmd.dirname);
 if(stat(mapfile,&buf)){ printf("##ERROR MAP: %s\n",mapfile); return 0; }
 if(readmrc(&ATOM_mrc[0],mapfile))
   return(0);
 //N
 sprintf(mapfile,"%s/atom_N.mrc",cmd.dirname);
 if(stat(mapfile,&buf)){ printf("##ERROR MAP: %s\n",mapfile); return 0; }
 if(readmrc(&ATOM_mrc[1],mapfile))
   return(0);
 //CA
 sprintf(mapfile,"%s/atom_CA.mrc",cmd.dirname);
 if(stat(mapfile,&buf)){ printf("##ERROR MAP: %s\n",mapfile); return 0; }
 if(readmrc(&ATOM_mrc[2],mapfile))
   return(0);
 //C
 sprintf(mapfile,"%s/atom_C.mrc",cmd.dirname);
 if(stat(mapfile,&buf)){ printf("##ERROR MAP: %s\n",mapfile); return 0; }
 if(readmrc(&ATOM_mrc[3],mapfile))
   return(0);
 //O
 sprintf(mapfile,"%s/atom_O.mrc",cmd.dirname);
 if(stat(mapfile,&buf)){ printf("##ERROR MAP: %s\n",mapfile); return 0; }
 if(readmrc(&ATOM_mrc[4],mapfile))
   return(0);
 //CB
 sprintf(mapfile,"%s/atom_CB.mrc",cmd.dirname);
 if(stat(mapfile,&buf)){ printf("##ERROR MAP: %s\n",mapfile); return 0; }
 if(readmrc(&ATOM_mrc[5],mapfile))
   return(0);


 //DAQscoring
 PDB APDB;
 int Natm=CountAtom(cmd.Afilename);
 printf("#Natm= %d\n",Natm);
 MallocPdb(&APDB,Natm);
 readpdb(&APDB,cmd.Afilename,Natm);
 RESULT_MODEL **model,tmp_model;
 int Nmodel=0;
 
 if((g.node=(NODE*)malloc(sizeof(NODE)*Natm))==NULL)
  return 0;
 g.Nnode=0;
 //Compute NR node positions
 for(int i=0;i<Natm;i++){
	//Copy to node
	 g.node[g.Nnode].real_cd[0] = APDB.xyz[i][0];
	 g.node[g.Nnode].real_cd[1] = APDB.xyz[i][1];
	 g.node[g.Nnode].real_cd[2] = APDB.xyz[i][2];
	 g.Nnode++;
	
 }
 printf("#Nnode= %d\n",g.Nnode);

 //Volume and density data
 //AssignProbToNodeSimple(&ptbl,&g,&mrc);
 AssignUnetProbToNode(ATOM_mrc,AA_mrc,&g);

 //Scoring Models
/*
 for(int i=0;i<Natm;i++){
  float SumScore = 0.0;
  int AtomOnRes =  APDB.AtomOnRes[i];
  int aa_type   =  APDB.TypeResId[AtomOnRes];
  int atm_type   =  APDB.TypeAtomId[i];


  //Check Score
  float aa_sco  = g.node[i].LogAA[aa_type];
  float atm_sco = g.node[i].LogATOM[atm_type];
  printf("%d %s %d %d aa_sco= %.3f atom_sco= %.3f\n",i,APDB.TypeAtom[i],atm_type,aa_type,aa_sco,atm_sco);
 }
*/
 ShowDAQscore(&g,&APDB);
 t4=gettimeofday_sec();
 printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
 return 0;

}

float Dist2(NODE *Ca1,NODE *Ca2){
	float d2=(Ca1->real_cd[0]-Ca2->real_cd[0])
		*(Ca1->real_cd[0]-Ca2->real_cd[0])
		+(Ca1->real_cd[1]-Ca2->real_cd[1])
		*(Ca1->real_cd[1]-Ca2->real_cd[1])
		+(Ca1->real_cd[2]-Ca2->real_cd[2])
		*(Ca1->real_cd[2]-Ca2->real_cd[2]);

 return d2;
}


//OutPut as PDB format



bool ShowDAQscore(GRAPH *g, PDB *p){
 float sum_all_aa=0.0;
 float sum_ca_aa=0.0;
 float sum_bb_aa=0.0;
 float sum_side_aa=0.0;
 float sum_all_atm=0.0;
 float sum_ca_atm=0.0;
 float sum_bb_atm=0.0;
 float sum_side_atm=0.0;
 int Nall,Nca,Nbb,Nside;
 Nall=Nca=Nbb=Nside=0;
 int Natm=0;
 for(int i=0;i<g->Nnode;i++){
  int AtomOnRes =  p->AtomOnRes[i];
  int aa_type   =  p->TypeResId[AtomOnRes];
  int atm_type   =  p->TypeAtomId[i];

  //Check Score
  float aa_sco  = g->node[i].LogAA[aa_type];
  float atm_sco = g->node[i].LogATOM[atm_type];
  float tmp[3];
  //Other,N,CA,C,O,CB
  //Side-chain Other, CB
  if(atm_type==0 || atm_type==5){
   Nside++;
   sum_side_aa += aa_sco;
   sum_side_atm += atm_sco;
  }
  //bb N,CA,C,O
  if(atm_type>0 && atm_type<5){
   Nbb++;
   sum_bb_aa += aa_sco;
   sum_bb_atm += atm_sco;
  }
  //CA
  if(atm_type==2){
   Nca++;
   sum_ca_aa += aa_sco;
   sum_ca_atm += atm_sco;
  }

	tmp[0]=p->xyz[i][0];
         tmp[1]=p->xyz[i][1];
         tmp[2]=p->xyz[i][2];
  	printf("ATOM  %5d  %3s %3s %c%4d    ",Natm++
	,p->TypeAtom[i],RES_NAMES[aa_type]
       	,p->Chain[i],p->ResNum[AtomOnRes]);
	printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,aa_sco);
  	printf("%d atom_sco= %.3f aa_sco= %.3f\n",AtomOnRes,atm_sco,aa_sco);

	//20AAs
	printf("#P(AA): ");
	for(int aa=0;aa<20;aa++)
	 printf("%3s %4.2f,",RES_NAMES[aa],g->node[i].AAP[aa]);
	puts("");

	//20AAs
	printf("#Log(AA): ");
	for(int aa=0;aa<20;aa++)
	 printf("%3s %.2f,",RES_NAMES[aa],g->node[i].LogAA[aa]);
	puts("");
 }

 puts("---------SUMMARY-----------");
 printf("TOTAL SIDE-CHAIN: ATOM_SCO= %.3f AA_SCO= %.3f\n",sum_side_atm,sum_side_aa);
 printf("AVG   SIDE-CHAIN: ATOM_SCO= %.3f AA_SCO= %.3f N= %d\n",sum_side_aa/Nside,sum_side_atm/Nside,Nside);
 printf("TOTAL BACK-BONE : ATOM_SCO= %.3f AA_SCO= %.3f\n",sum_bb_atm,sum_bb_aa);
 printf("AVG   BACK-BONE : ATOM_SCO= %.3f AA_SCO= %.3f N= %d\n",sum_bb_atm/Nbb,sum_bb_aa/Nbb,Nbb);
 printf("TOTAL Calpha    : ATOM_SCO= %.3f AA_SCO= %.3f\n",sum_ca_atm,sum_ca_aa);
 printf("AVG   Calpha    : ATOM_SCO= %.3f AA_SCO= %.3f N= %d\n",sum_ca_atm/Nca,sum_ca_aa/Nca,Nca);
 printf("TOTAL All       : ATOM_SCO= %.3f AA_SCO= %.3f\n",sum_side_atm+sum_bb_atm,sum_side_aa+sum_bb_aa);
 printf("AVG   All       : ATOM_SCO= %.3f AA_SCO= %.3f N= %d\n",(sum_bb_atm+sum_side_atm)/(Nbb+Nside),(sum_side_aa+sum_bb_aa)/(Nbb+Nside),Nbb+Nside);
}
