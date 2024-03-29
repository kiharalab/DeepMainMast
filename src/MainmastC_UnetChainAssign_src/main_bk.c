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
#include "thread.h"
#include "rmsd.h"
#include <sys/stat.h>

//Add Homo mode;
//#include "scoring.h"

#define PDB_STRLEN 55

//Fragment

typedef struct{
 NODE Ca[200];
 int Start,End,Len;
 float RawScore,Score,Rmsd;
 int id;
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


void CompareFragAmodel(SEQ_NODE *,SEQFG*,int,NODE *,bool);
bool CompareFragAmodelMore(SEQ_NODE *,SEQFG*,int,NODE *,bool);
bool MakeDummyFrag(SEQ_NODE *,SEQ_NODE *,SEQFG*,int,bool);
int AddModel(RESULT_MODEL *,RESULT_MODEL **, int);
int AddConnectModel(RESULT_MODEL **, int,int,SEQ *);
int CheckOverlap(int,int,int,int);
bool CheckConf(RESULT_MODEL *,RESULT_MODEL *,SEQ *);
bool CheckHomoShape(RESULT_MODEL *,RESULT_MODEL *,SEQ *,float);
bool CheckPairs(RESULT_MODEL *,RESULT_MODEL *,RESULT_MODEL *,RESULT_MODEL *,SEQ *,float);

bool CheckContact(RESULT_MODEL *,RESULT_MODEL *,short int *,int);
float Dist2(NODE *,NODE *);

void TableShow(GRAPH *, SEQ *,MRC *,float);
float CheckConnect(NODE *,NODE *,MRC *,float);

int CountModels(PDB *,SEQ *,int);

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
	if(cmd.DummyMode==true)
         sprintf(mapfile,"%s/%s%s.mrc",cmd.dirname,"sigmoidAA_",RES_NAMES[aa]);
        else
         sprintf(mapfile,"%s/%s%s.mrc",cmd.dirname,"AA_",RES_NAMES[aa]);
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
 //O
 sprintf(mapfile,"%s/atom_CB.mrc",cmd.dirname);
 if(stat(mapfile,&buf)){ printf("##ERROR MAP: %s\n",mapfile); return 0; }
 if(readmrc(&ATOM_mrc[5],mapfile))
   return(0);


 SEQ seq;
 if(readseq(&seq,cmd.sfilename)){
  puts("Error in seq file");
  return 0;
 }

 //Model Mode!!
 PDB APDB,CPDB;
 int Ncatm=0;
 int Ncont=0;
 short int *CONTtbl;
 if(cmd.Cmode==true){
  printf("##Loading Contact PDB info....%s\n",cmd.Cfilename);
  Ncatm=CountAtom(cmd.Cfilename);
  printf("#Contact Natm= %d\n",Ncatm);
  MallocPdb(&CPDB,Ncatm);
  readpdb(&CPDB,cmd.Cfilename,Ncatm);

  CONTtbl=(short int *)calloc(seq.len*seq.len,sizeof(short int));

	for(int i=1;i<CPDB.NumOfAtom-1;i++){
	 int ResNum1 = CPDB.ResNum[i];//Start from 1
  	 int order = CPDB.AtomOnRes[i];
	 if(CPDB.TypeResId[order]!=seq.seq_code[ResNum1-1]){
	  printf("#Sequence Miss: %d != %d in %d-th ATOM\n",CPDB.TypeResId[order],seq.seq_code[ResNum1-1],i);
	  return 0;
	 }
		//Contact CHeck
		for(int j=i+5;j<CPDB.NumOfAtom-1;j++){
	 	 int ResNum2 = CPDB.ResNum[j];//Start from 1
		 if(CPDB.Chain[i] != CPDB.Chain[j])
		  continue;
		 float d1=L(CPDB.xyz[i][0],CPDB.xyz[i][1],CPDB.xyz[i][2],
			    CPDB.xyz[j][0],CPDB.xyz[j][1],CPDB.xyz[j][2]);
		 if(d1>=5.0)
	  	  continue;
		 printf("#Contact %d:%d %.3f\n",i,j,d1);
		 CONTtbl[(ResNum1-1)*seq.len+(ResNum2-1)]=1;//Contact
		}
	}
 }

 int Natm=CountAtom(cmd.Afilename);
 printf("#Natm= %d\n",Natm);
 MallocPdb(&APDB,Natm);
 readpdb(&APDB,cmd.Afilename,Natm);
 RESULT_MODEL **model,tmp_model;
 int Nmodel=0;
 int Ncnt_models = CountModels(&APDB,&seq,cmd.FragLen);
 //ADD MALLOC AND NR computation
 if((model=(RESULT_MODEL **)malloc(sizeof(RESULT_MODEL *)*Ncnt_models))==NULL)
  return 0;
 for(int i=0;i<Ncnt_models;i++)
  if((model[i]=(RESULT_MODEL *)malloc(sizeof(RESULT_MODEL)))==NULL)
   return 0;

 //convert PDB to model Only CA
 for(int i=0;i<APDB.NumOfAtom - cmd.FragLen ;i++){
  int ResNum = APDB.ResNum[i];//Start from 1
  int order = APDB.AtomOnRes[i];
		
  //Check Sequence consistency
  if(APDB.TypeResId[order]!=seq.seq_code[ResNum-1])
   continue;
   
	//Check fragment
	int Ncd=0;
	for(int j=0;j<cmd.FragLen && i+j < APDB.NumOfAtom ;j++){
	 //Check Residue Number
	 if(ResNum + j != APDB.ResNum[i+j])
	  continue;
	 tmp_model.Ca[j].real_cd[0]=APDB.xyz[i+j][0];
	 tmp_model.Ca[j].real_cd[1]=APDB.xyz[i+j][1];
	 tmp_model.Ca[j].real_cd[2]=APDB.xyz[i+j][2];
	 tmp_model.Ca[Ncd].node_id = -1;
	 Ncd++;
	}
	if(Ncd != cmd.FragLen)
	 continue;
	tmp_model.Start=ResNum-1;
	tmp_model.End=tmp_model.Start + cmd.FragLen;
	tmp_model.Score=0;//Not yet
	tmp_model.RawScore=0;//not yet
	tmp_model.Rmsd=0;
	tmp_model.Len=Ncd;
	Nmodel=AddModel(&tmp_model,model,Nmodel);

	//Check Homo-oligomer fragments
	for(int pos=0;pos<seq.len-cmd.FragLen;pos++){
		if(pos==tmp_model.Start)
		 continue;
		bool homo_check=true;
		for(int k=0;k<cmd.FragLen;k++){
		 if(seq.seq_code[pos+k]!=seq.seq_code[tmp_model.Start+k]){
		  homo_check=false;
		  break;
		 }
		}
		
		if(homo_check==true){
		 //printf("Found Homo %d %d\n",tmp_model.Start,pos);
		 tmp_model.Start=pos;
		 tmp_model.End=tmp_model.Start + cmd.FragLen;
		 Nmodel=AddModel(&tmp_model,model,Nmodel);
		 //Back
		 tmp_model.Start=ResNum-1;
		 tmp_model.End=tmp_model.Start + cmd.FragLen;
		}
	}

 }
 
 printf("#Nmodel= %d\n",Nmodel);
 return 0;
 if((g.node=(NODE*)malloc(sizeof(NODE)*Nmodel*cmd.FragLen))==NULL)
  return 0;
 g.Nnode=0;
 //Compute NR node positions
 for(int i=0;i<Nmodel;i++){
 	for(int j=0;j<model[i]->Len;j++){
	 //printf("Model%d pos%d x= %f\n",i,j,model[i]->Ca[j].real_cd[0]);
	 //Check redundant node (or position)
		bool same = false;
		for(int k=0;k<g.Nnode;k++){
		 float d2 = Dist2(&(g.node[k]),&(model[i]->Ca[j]));
		 if(d2 < 1.0){
		  same = true;
		  model[i]->Ca[j].node_id = k;
		  break;
		 }
		}
		//Add
		if(same == false){
		 g.node[g.Nnode].real_cd[0] = model[i]->Ca[j].real_cd[0];
		 g.node[g.Nnode].real_cd[1] = model[i]->Ca[j].real_cd[1];
		 g.node[g.Nnode].real_cd[2] = model[i]->Ca[j].real_cd[2];
		 model[i]->Ca[j].node_id = g.Nnode;
		 g.Nnode++;
		}
	}
 }
 printf("#Nnode= %d\n",g.Nnode);

 //Volume and density data
 //AssignProbToNodeSimple(&ptbl,&g,&mrc);
 AssignUnetProbToNode(ATOM_mrc,AA_mrc,&g);

 float DAQcutoff=cmd.Filter*cmd.FragLen;

 //Scoring Models
 for(int i=0;i<Nmodel;i++){
  float SumScore = 0.0;
 	for(int j=0;j<model[i]->Len;j++){
	 int node_id = model[i]->Ca[j].node_id;
	 int res_order = model[i]->Start + j;
	 int aa_type   =  seq.seq_code[res_order];

	 //Check Score
	 SumScore += g.node[node_id].LogAA[aa_type];

	}
  model[i]->RawScore=SumScore;
  //if(SumScore < 0.0) //Ignore negative scores
   //model[i]->RawScore=0.01;
   
  //printf("%d Sco= %.3f\n",i,SumScore);
 }
 //return 0;
 //For fragments-------------------------------------------
 //Show SEQ
 for(int p=0;p<seq.len;p++){
  //printf("SEQ %d %s\n",1+p,RES_NAMES[seq.seq_code[p]]);
  printf("SEQ %d %d %s\n",1+p,seq.chain_id[p],RES_NAMES[seq.seq_code[p]]);
 }
 //Node coordinates
 for(int m=0;m<g.Nnode;m++){
  printf("NODE %8.3f %8.3f %8.3f\n",g.node[m].real_cd[0],g.node[m].real_cd[1],g.node[m].real_cd[2]);
 }
 //Simplify
 int Nmodel2=0;
 for(int m=0;m<Nmodel;m++){
  if(model[m]->RawScore < DAQcutoff)//Remove low scores
   continue;
  if(model[m]->RawScore < 0.0)//set 0.01 to negative score
   model[m]->RawScore=0.01;
  printf("FRG %d RAWSCO= %f ZSCORE= %f RMSD= %.3f ALI= ",
  Nmodel2,model[m]->RawScore,model[m]->Score,model[m]->Rmsd);
  //m,model[m]->RawScore,model[m]->Score,model[m]->Rmsd);
  for(int p=0;p<model[m]->Len;p++){
   printf("%d,%d,",model[m]->Ca[p].node_id,model[m]->Start+1+p);
  }
  printf("\n");
  //Copy data
  if(m != Nmodel2){
   model[Nmodel2]=model[m];
  }
  Nmodel2++;
 }
 puts("#DONE FRG");
/*
 printf("1S %d E %d L %d\n",model[936]->Start,model[936]->End,model[936]->Len);
 printf("2S %d E %d L %d\n",model[939]->Start,model[939]->End,model[939]->Len);
 printf("Dist, %f\n",sqrt(Dist2(&(model[936]->Ca[9]),&(model[939]->Ca[0]))));
 printf("cd1, %.3f %.3f %.3f\n",model[936]->Ca[9].real_cd[0],model[936]->Ca[9].real_cd[1],model[936]->Ca[9].real_cd[2]);
 printf("cd2, %.3f %.3f %.3f\n",model[939]->Ca[0].real_cd[0],model[939]->Ca[0].real_cd[1],model[939]->Ca[0].real_cd[2]);
*/
 //return 0;
 //check matrix
 for(int m1=0;m1<Nmodel2;m1++){
  printf("MTX %d ",m1);
  for(int m2=0;m2<Nmodel2;m2++){
	if (m1==m2)
	 continue;
	if(CheckConf(model[m1],model[m2],&seq)==true){
	 printf("%d,",m2);
	 continue;
	}
	if(cmd.Cmode==true && CheckContact(model[m1],model[m2],CONTtbl,seq.len)==true){
	 printf("%d,",m2);
	 continue;
	}
	//new for homo-oligomer
	if(cmd.Hmode==true){
	 if(CheckHomoShape(model[m1],model[m2],&seq,cmd.RMSD)==true){
	  printf("%d,",m2);
	  continue;
	 }
	}
  }
   puts("");
 }
 //Pair fragments for Homo****
 if(cmd.Hmode==true){
 //if(false){
  printf("##Homockecking Nmodel= %d\n",Nmodel2);

  RESULT_MODEL **pair_list;
  int Npair=0;
  if((pair_list=(RESULT_MODEL **)malloc(sizeof(RESULT_MODEL *)*(int)(Nmodel2*0.25)*(int)(Nmodel2*0.25)))==NULL)
   return false;



  for(int m11=0;m11<Nmodel2;m11++)
   model[m11]->id=m11;//init id
  for(int m11=0;m11<Nmodel2;m11++){
   int chid11=seq.chain_id[model[m11]->Start];
   int pos11=model[m11]->Start - seq.chain_start[chid11];
   if(pos11%4!=0)
    continue;
   	for(int m12=m11+1;m12<Nmodel2;m12++){
   	 int chid12=seq.chain_id[model[m12]->Start];
   	 if(chid11 !=chid12)
   	  continue;
   	 int pos12=model[m12]->Start - seq.chain_start[chid12];
   	 if(pos12%4!=0)
   	  continue;

	 //no-overlap
	 if((pos11-pos12)*(pos11-pos12) < model[m11]->Len * model[m11]->Len)
	  continue;

   	 if(CheckConf(model[m11],model[m12],&seq)==true)
   	  continue;
	 
	 pair_list[2*Npair]=model[m11];
	 pair_list[2*Npair+1]=model[m12];
    	 Npair++;
	}
  }
  printf("##Npair= %d\n",Npair);
  RESULT_MODEL *f11,*f12,*f21,*f22;	
  for(int pr1=0;pr1<Npair;pr1++){
   f11 = pair_list[2*pr1];
   f12 = pair_list[2*pr1+1];
   int chid1=seq.chain_id[f11->Start];
  	for(int pr2=pr1+1;pr2<Npair;pr2++){
   	 f21 = pair_list[2*pr2];
   	 f22 = pair_list[2*pr2+1];
   	 int chid2=seq.chain_id[f21->Start];
	 if(chid1==chid2)
	  continue;
	  
	 if(CheckPairs(f11,f12,f21,f22,&seq,cmd.RMSD)==true){
		 if(CheckConf(f11,f21,&seq)==true)
		   continue;
		  if(CheckConf(f11,f22,&seq)==true)
		   continue;
		  if(CheckConf(f12,f21,&seq)==true)
		   continue;
		  if(CheckConf(f12,f22,&seq)==true)
		   continue;
		  printf("PAIR %d,%d %d,%d\n",f11->id,f12->id,f21->id,f22->id);
	  //printf("PAIR %d:%d %d-%d vs %d-%d ch: %d %d\n",pr1,pr2,f11->Start,f12->Start,f21->Start,f22->Start,chid1,chid2);
     }
	}
  }
 }
 t4=gettimeofday_sec();
 printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
 return 0;

}

bool CheckPairs(RESULT_MODEL *f11,RESULT_MODEL *f12,RESULT_MODEL *f21,RESULT_MODEL *f22, SEQ *seq, float rcut){
/*
f11--f12
    f21----f22
*/
 bool seq_flag11,seq_flag12;
 bool seq_flag21,seq_flag22;
 bool seq_flag;
 float r2 = rcut*rcut;
 int chid1 = seq->chain_id[f11->Start];
 int chid2 = seq->chain_id[f21->Start];
 int pos11=f11->Start - seq->chain_start[chid1];
 int pos12=f12->Start - seq->chain_start[chid1];
 int pos21=f21->Start - seq->chain_start[chid2];
 int pos22=f22->Start - seq->chain_start[chid2];
 int st1,st2;

 //Position check
 //if(((pos11==pos21) && (pos12==pos22)) || ((pos11==pos21) && (pos12==pos22)))

 //check homomer?
 //f11 vs f21 and f12 vs f22
 //SeqID check
 if((pos11==pos21) && (pos12==pos22)){
  	st1=f11->Start;
  	st2=f21->Start;
 	for(int i=0;i<f11->Len;i++){
         int p1=st1+i;
         int p2=st2+i;
         //AAcheck
         if(seq->seq_code[p1]!=seq->seq_code[p2]){
	  return false;
	 }
 	}
 	st1=f12->Start;
 	st2=f22->Start;
 	for(int i=0;i<f11->Len;i++){
         int p1=st1+i;
         int p2=st2+i;
         //AAcheck
         if(seq->seq_code[p1]!=seq->seq_code[p2]){
	  return false;
         }
 	}
 }else if((pos11==pos22) && (pos12==pos21)){

 	st1=f11->Start;
 	st2=f22->Start;
 	for(int i=0;i<f11->Len;i++){
         int p1=st1+i;
         int p2=st2+i;
         //AAcheck
         if(seq->seq_code[p1]!=seq->seq_code[p2]){
	  return false;
	 }
 	}
 	st1=f12->Start;
 	st2=f21->Start;
 	for(int i=0;i<f11->Len;i++){
         int p1=st1+i;
         int p2=st2+i;
         //AAcheck
         if(seq->seq_code[p1]!=seq->seq_code[p2]){
	  return false;
         }
 	}
 }else{
  return false;
 }

 //Check Distance Consistency
 //head-head, tail-tail
 float Dhh1,Dhh2,Dtt1,Dtt2;
 Dhh1=Dist2(&(f11->Ca[0]),&(f12->Ca[0]));
 Dhh2=Dist2(&(f21->Ca[0]),&(f22->Ca[0]));

 if((Dhh1-Dhh2)*(Dhh1-Dhh2)>r2)
  return true;
 int L=f11->Len-1;
 Dtt1=Dist2(&(f11->Ca[L]),&(f12->Ca[L]));
 Dtt2=Dist2(&(f21->Ca[L]),&(f22->Ca[L]));
 if((Dtt1-Dtt2)*(Dtt1-Dtt2)>r2)
  return true;


 return false;
}



void TableShow(GRAPH *g,SEQ *s,MRC *mrc,float Pcut){

 int AAcd;
 float AAsco,ATOMsco,SSsco;
 float Ph,Pe,Pc,d;

 NODE *TBL[10000];//10k points
 int Npoints=0;

 //Show coordinates in PDB format
 for(int i=0;i<g->Nnode;i++){
/*
  if(g->node[i].ATOMP[1] + g->node[i].ATOMP[2] + g->node[i].ATOMP[3] < Pcut){
   continue;
  }
  printf("ATOM  %5d  CA  %3s%6d    ",Npoints+1,RES_NAMES[0],Npoints+1);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",g->node[i].real_cd[0],g->node[i].real_cd[1],g->node[i].real_cd[2],
1.0,1.0);
*/
  TBL[Npoints]=&(g->node[i]);
  Npoints++;
 }
 printf("Npoints= %d Nseq= %d\n",Npoints,s->len);
 //Show SEQ
 //for(int p=0;p<s->len;p++){
 // printf("SEQ %d %s\n",1+p,RES_NAMES[s->seq_code[p]]);
 //}


 //Distance Matrix = Edge data
 int STOCK[2000];
 int Nstock=0;
 for(int i=0;i<Npoints;i++){
        for(int j=i+1;j<Npoints;j++){
         d=Dist2(TBL[i],TBL[j]);
         d=sqrt(d);
         if(d>10.0){ //Too far > 10A
          printf("DMTX %d %d %.0f %.0f\n",i+1,j+1,d*10.0,0.00);//Int Format*10
          continue;
         }
         //Check Connection
         float MinP = CheckConnect(TBL[i],TBL[j],mrc,Pcut);
         printf("DMTX %d %d %.0f %.0f\n",i+1,j+1,d*10.0,MinP*100);//Int Format*10
        }
 }
 //Sequence Node and LDP Node
 for(int p=0;p<s->len;p++){//Seq
        for(int i=0;i<Npoints;i++){//Node
         AAcd=s->seq_code[p];
         Ph=s->Pss[p][0];
         Pe=s->Pss[p][1];
         Pc=s->Pss[p][2];

         AAsco=TBL[i]->LogAA[AAcd];
         ATOMsco=TBL[i]->LogATOM[2];//Other,N,CA,C,O,CB
         SSsco=Ph*TBL[i]->LogSS[0];
             +Ph*TBL[i]->LogSS[1];
             +Pc*TBL[i]->LogSS[2];
         //printf("SCO %d %d %f\n",p+Npoints+1,i+1,AAsco+ATOMsco+SSsco);
         //SEQ, NODE
         printf("SCO %d %d %f %f %f\n",p+1,i+1,ATOMsco,AAsco,SSsco);
        }
 }
}


//Check Homo position and RMSD
bool CheckHomoShape(RESULT_MODEL *f1,RESULT_MODEL *f2, SEQ *seq, float rcut){
 //Same chain
 if(seq->chain_id[f1->Start]==seq->chain_id[f2->Start])
  return false;
 
 int chid1=seq->chain_id[f1->Start];
 int chid2=seq->chain_id[f2->Start];
 int chst1=seq->chain_start[chid1];
 int chst2=seq->chain_start[chid2];
 int st1 = f1->Start;
 int st2 = f2->Start;

 if(st1 - chst1 != st2 - chst2) // Different chain positions
  return false;

 //SeqID check
 for(int i=0;i<f1->Len;i++){
	 int p1=st1+i;
	 int p2=st2+i;
	 //AAcheck
	 if(seq->seq_code[p1]!=seq->seq_code[p2])
	 	return false;
 }
 //printf("###Start RMSD comp\n");
 double rmsd=0;
 double cd1[100][3],cd2[100][3];
 for(int i=0;i<f1->Len;i++){
	 cd1[i][0]=f1->Ca[i].real_cd[0];
	 cd1[i][1]=f1->Ca[i].real_cd[1];
	 cd1[i][2]=f1->Ca[i].real_cd[2];
	 
	 cd2[i][0]=f2->Ca[i].real_cd[0];
	 cd2[i][1]=f2->Ca[i].real_cd[1];
	 cd2[i][2]=f2->Ca[i].real_cd[2];
	 
	 //printf("cd1= %3.1f %3.1f %3.1f\n",cd1[i][0],cd1[i][1],cd1[i][2]);
	 //printf("cd2= %3.1f %3.1f %3.1f\n",cd2[i][0],cd2[i][1],cd2[i][2]);

	 

 }
 fast_rmsd(cd1,cd2,f1->Len,&rmsd);
 //printf("GotRMSD= %f %f\n",rmsd,rcut);
 if(rmsd>rcut)
  return true;


 return false;
}


bool CheckConf(RESULT_MODEL *f1,RESULT_MODEL *f2, SEQ *seq){
 //int mode = CheckOverlap(f1->Start,f1->End,f2->Start,f2->End);
 float dist2=0.0;
 float ClashDist = 2.5*2.5;
 float ShiftDist = 3.5*3.5;
 float SepDist = 5.0;
 //float SepDist = 4.0;
 int SeqDiff;
 //Connection Distance Check

 //Same chain
 if(seq->chain_id[f1->Start]==seq->chain_id[f2->Start]){

  //Start1 - Start2
  SeqDiff = f2->Start - f1->Start;
  dist2 = Dist2(&(f1->Ca[0]),&(f2->Ca[0]));
  if(dist2 > (SepDist*SeqDiff)*(SepDist*SeqDiff))
   return true;

  //Start1 - End2
  SeqDiff = f1->Start - (f2->End-1);
  dist2 = Dist2(&(f1->Ca[0]),&(f2->Ca[f2->Len-1]));
  if(dist2 > (SepDist*SeqDiff)*(SepDist*SeqDiff))
   return true;

  //End1 - Start2
  SeqDiff = (f1->End -1) - f2->Start;
  dist2 = Dist2(&(f1->Ca[f1->Len-1]),&(f2->Ca[0]));
  if(dist2 > (SepDist*SeqDiff)*(SepDist*SeqDiff))
   return true;

  //End1 - End2
  SeqDiff = (f1->End -1) - (f2->End -1);
  dist2 = Dist2(&(f1->Ca[f1->Len-1]),&(f2->Ca[f2->Len-1]));
  if(dist2 > (SepDist*SeqDiff)*(SepDist*SeqDiff))
   return true;
 }

 //puts("Overlap");
 for(int a=0;a<f1->Len;a++){
  int pos1 = f1->Start+a;
  for(int b=0;b<f2->Len;b++){
   int pos2 = f2->Start+b;
   //printf("%d:%d\n",a,b);
   dist2=Dist2(&(f1->Ca[a]),&(f2->Ca[b]));
   //Overlap Check
   if(pos1==pos2 && dist2 > ShiftDist)
    return true;
   else if(pos1!=pos2 && dist2 < ClashDist)//Clash Check
    return true;
  }
 }

 //puts("DONE");
 return false;
}


//Compare with AF2 model
void CompareFragAmodel(SEQ_NODE *result,SEQFG *sfg,int Nsfg,NODE *mod,bool mode){
 for(int s=0;s<Nsfg;s++){
  int N=result[s].Nali;
  int posi=sfg[s].pos;
  double cd1[100][3],cd2[100][3];
  double mtx[3][3],mov_com[3],mov_ref[3];
  int Ncd;
  printf("N= %d\n",N);
  if(N==0)
   continue;
  SEQ_NODE *res=&result[s];
	for(int i=0;i<N;i++){
	 //copy coordinates
	 Ncd=0;
	 for(int j=0;j<sfg[s].l;j++){
	  if(res->ali[i][j]==NULL)
	   continue;

	   cd1[Ncd][0]=res->ali[i][j]->real_cd[0];
	   cd1[Ncd][1]=res->ali[i][j]->real_cd[1];
	   cd1[Ncd][2]=res->ali[i][j]->real_cd[2];

	   cd2[Ncd][0]=mod[posi+j].real_cd[0];
	   cd2[Ncd][1]=mod[posi+j].real_cd[1];
	   cd2[Ncd][2]=mod[posi+j].real_cd[2];

	  Ncd++;

	 }
	 if(Ncd<4){//too short
	  res->rmsd[i]=100;
	  continue;
	 }
	 //check rmsd
	 double rmsd=0;
	 fast_rmsd(cd1,cd2,Ncd,&rmsd);
	 //printf("%d %d RMSD= %f\n",s,i,rmsd);
	 res->rmsd[i]=rmsd;
	 if(!isfinite(rmsd))
	  res->rmsd[i]=100;


	 if(mode==true && res->rmsd[i]<cmd.RMSD && res->score[i] >= cmd.zcut){
	  float rotated[3],tmp[3];
		Ncd=0;
		for(int j=0;j<sfg[s].l;j++){
	 	 if(res->ali[i][j]==NULL)
	 	  continue;
	 	 cd1[Ncd][0]=res->ali[i][j]->real_cd[0];
	   	 cd1[Ncd][1]=res->ali[i][j]->real_cd[1];
	  	 cd1[Ncd][2]=res->ali[i][j]->real_cd[2];
	  	 cd2[Ncd][0]=mod[posi+j].real_cd[0];
	  	 cd2[Ncd][1]=mod[posi+j].real_cd[1];
	  	 cd2[Ncd][2]=mod[posi+j].real_cd[2];
	  	 Ncd++;
	 	}
	  	calculate_rotation_rmsd(cd1,cd2,Ncd,mov_com,mov_ref,mtx,&rmsd);
		//printf("com %f %f %f\n",mov_com[0],mov_com[1],mov_com[2]);
		//printf("ref %f %f %f\n",mov_ref[0],mov_ref[1],mov_ref[2]);
	 	printf("##POS %d RAWSCO= %f ZSCORE= %f RMSD= %f Amodel\n",sfg[s].pos,res->raw_score[i],res->score[i],res->rmsd[i]);
	 	puts("MODEL");
		//move rotate and extend
		for(int j=0;j<sfg[s].l;j++){
		 if(res->ali[i][j]==NULL)
	 	  continue;
		  //mov -> o
		  tmp[0]=mod[posi+j].real_cd[0]-mov_com[0];
		  tmp[1]=mod[posi+j].real_cd[1]-mov_com[1];
		  tmp[2]=mod[posi+j].real_cd[2]-mov_com[2];
			for(int i1=0;i1<3;i1++){
				rotated[i1]=0.0;
				for(int i2=0;i2<3;i2++)
				 rotated[i1] +=mtx[i1][i2]*tmp[i2];
			}
		  //o -> mov -> ref
		  tmp[0]=rotated[0]+mov_com[0]+mov_ref[0];
		  tmp[1]=rotated[1]+mov_com[1]+mov_ref[1];
		  tmp[2]=rotated[2]+mov_com[2]+mov_ref[2];
		 printf("ATOM  %5d  CA  %3s%6d    ",j+1,RES_NAMES[sfg[s].seq[j]],sfg[s].pos+j+1);
  		 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,res->score[i]);
		}
 	 	puts("ENDMDL");
	 }
	}
 }
}

int CheckOverlap(int st1,int ed1,int st2,int ed2){
/*
     st1------ed1
st2-----------------ed2
or
     st1-----------------ed1
st2-----------------ed2
*/
 if(st2<=st1 && st1<ed2)
  return -1;
/*
st1---------------ed1
     st2------ed2
or
st1--------------ed1
     st2-----------------ed2
*/
 if(st1<=st2 && st2<ed1)
  return -1;

/*
st1----ed1
            st2------ed2
*/
 if(ed1<=st2)
  return 0;
/*
             st1----ed1
st2------ed2
*/
 if(ed2<=st1)
  return 1;
 return -1;
}

bool CheckGeo(RESULT_MODEL *f1,RESULT_MODEL *f2){

 //Check N and C-ter
 //N-ter
 int Cter=f1->Len-1;
 float d2;
 for(int p=0;p<f2->Len;p++){
	 d2=(f2->Ca[p].real_cd[0]-f1->Ca[0].real_cd[0])
	   *(f2->Ca[p].real_cd[0]-f1->Ca[0].real_cd[0])
	   +(f2->Ca[p].real_cd[1]-f1->Ca[0].real_cd[1])
	   *(f2->Ca[p].real_cd[1]-f1->Ca[0].real_cd[1])
	   +(f2->Ca[p].real_cd[2]-f1->Ca[0].real_cd[2])
	   *(f2->Ca[p].real_cd[2]-f1->Ca[0].real_cd[2]);
	if(d2 <1.0){
	 //printf("Nter-connect\n");
	 return true;
	}
	d2=(f2->Ca[p].real_cd[0]-f1->Ca[Cter].real_cd[0])
	   *(f2->Ca[p].real_cd[0]-f1->Ca[Cter].real_cd[0])
	   +(f2->Ca[p].real_cd[1]-f1->Ca[Cter].real_cd[1])
	   *(f2->Ca[p].real_cd[1]-f1->Ca[Cter].real_cd[1])
	   +(f2->Ca[p].real_cd[2]-f1->Ca[Cter].real_cd[2])
	   *(f2->Ca[p].real_cd[2]-f1->Ca[Cter].real_cd[2]);

	if(d2<1.0){
	 //printf("Cter-connect\n");
	 return true;
	}

 }
 return false;
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


bool Connect(RESULT_MODEL *f1,RESULT_MODEL *f2,RESULT_MODEL *f3,SEQ *seq,bool show){
 int Npoint,Cpoint;
 float d2,D2;
 float d_ng2 = 3.0*3.0;
 Npoint=-1;
 Cpoint=-1;
 //check overlap between f1 and f2
 for(int i=0;i<f1->Len;i++){
  for(int j=0;j<f2->Len;j++){
   d2=Dist2(&(f1->Ca[i]),&(f2->Ca[j]));
   if(d2 < d_ng2)
    return false;
  }
 }
 for(int i=0;i<f3->Len;i++){
  d2=Dist2(&(f1->Ca[f1->Len-1]),&(f3->Ca[i]));
  D2=Dist2(&(f2->Ca[0]),&(f3->Ca[i]));
  if(d2<=1.0)
   Npoint=i;
  if(D2<=1.0)
   Cpoint=i;
 }
 if(Npoint <0||Cpoint<0||Npoint==Cpoint)
  return false;

 //check Nca
 int Lmiss=f2->Start - f1->End+1;
 int Lfill = Cpoint - Npoint;
 if(Lfill<0)
  Lfill*=-1;
 if(Lfill != Lmiss)
  return false;

 //Keep insert coords
 NODE ins[30];
 int pp=0;
 int Lins=0;
 if(Npoint<Cpoint){
  for(int p=Npoint+1;p<Cpoint;p++){
   ins[pp].real_cd[0]=f3->Ca[p].real_cd[0];
   ins[pp].real_cd[1]=f3->Ca[p].real_cd[1];
   ins[pp].real_cd[2]=f3->Ca[p].real_cd[2];
   pp++;
  }
 }else{
  for(int p=Npoint-1;p>Cpoint;p--){
   ins[pp].real_cd[0]=f3->Ca[p].real_cd[0];
   ins[pp].real_cd[1]=f3->Ca[p].real_cd[1];
   ins[pp].real_cd[2]=f3->Ca[p].real_cd[2];
   pp++;
  }
 }
 Lins=pp;
 //Check Clash between ins and f1,f2
 for(int i=0;i<Lins;i++){
  for(int j=0;j<f1->Len;j++){
   d2=Dist2(&(ins[i]),&(f1->Ca[j]));
   if(d2 < d_ng2)
    return false;
  }
  for(int j=0;j<f2->Len;j++){
   d2=Dist2(&(ins[i]),&(f2->Ca[j]));
   if(d2 < d_ng2)
    return false;
  }
 }

 //printf("Find!! %d-%d\n",Npoint,Cpoint);
 //Show PDB format 
 if(show==true){
  printf("f1 %d-%d\n",f1->Start,f1->End);
  printf("f2 %d-%d\n",f2->Start,f2->End);
  printf("In %d-%d\n",Npoint,Cpoint);
  printf("##POS %d RAWSCO= %f ZSCORE= %f RMSD= %f Extended\n",
  f1->Start,f1->RawScore+f2->RawScore,f1->Score+f2->Score,0.00);
  puts("MODEL");
  int Natm,Nres;
  Natm=1;
  Nres=f1->Start;
	for(int p=0;p<f1->Len;p++){
	 printf("ATOM  %5d  CA  %3s%6d    ",Natm++,RES_NAMES[f1->Ca[p].AAtype],Nres+1+p);
  	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",f1->Ca[p].real_cd[0],f1->Ca[p].real_cd[1],f1->Ca[p].real_cd[2],1.0,f1->Score);
	}
	//Insert
	for(int p=0;p<Lins;p++){
	 printf("ATOM  %5d  CA  %3s%6d    ",Natm++,RES_NAMES[seq->seq_code[f1->End+p]],f1->End+1+p);
  	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f X\n",ins[p].real_cd[0],ins[p].real_cd[1],ins[p].real_cd[2],1.0,0.0);
	}
	//f2
  	Nres=f2->Start;
	for(int p=0;p<f2->Len;p++){
	 printf("ATOM  %5d  CA  %3s%6d    ",Natm++,RES_NAMES[f2->Ca[p].AAtype],Nres+1+p);
  	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",f2->Ca[p].real_cd[0],f2->Ca[p].real_cd[1],f2->Ca[p].real_cd[2],1.0,f2->Score);
	}
  puts("ENDMDL");

 }
 return true;
}



bool FindConnects(RESULT_MODEL *f1,RESULT_MODEL **f1_tbl,int n1,RESULT_MODEL *f2,RESULT_MODEL **f2_tbl,int n2,SEQ *seq,bool show){

 RESULT_MODEL *F1,*F2;
 RESULT_MODEL **F1_TBL,**F2_TBL;
 int N1,N2,Lmiss;
 if(f1->End <=f2->Start){
  F1=f1; F1_TBL=f1_tbl; N1=n1;
  F2=f2; F2_TBL=f2_tbl; N2=n2;
 }else{
  F2=f1;
  F2_TBL=f1_tbl;
  N2=n1;
  F1=f2;
  F1_TBL=f2_tbl;
  N1=n2;
 }
 Lmiss=F2->Start-F1->End+1; 
 if(Lmiss>20)
  return false;
 float Dmiss=Lmiss*3.70;
 float Dmiss2=Dmiss*Dmiss;
 float d2=Dist2(&(F1->Ca[F1->Len-1]),&(F2->Ca[0]));

 if(d2 > Dmiss2||d2==0.0)//Too far or overlap
  return false;

 for(int i=0;i<N1;i++){
  if(F1_TBL[i]->Len <Lmiss)
   continue;
  bool chk=false;
	for(int j=0;j<N2;j++){
	 if(F1_TBL[i]==F2_TBL[j]){
 	  chk=true;
	  break;
	 }
	}
  if(chk==false)
   continue;
  //printf("Check..%f\n",F1_TBL[i]->Score);
  if(Connect(F1,F2,F1_TBL[i],seq,show)==true){
	if(show==true){
   	 printf("#Connect....%d %d L=%d ",N1,N2,Lmiss);
   	 printf("%d-%d %d-%d %f : %f\n",F1->Start,F1->End,F2->Start, F2->End,d2,Dmiss2);
	}
   return true;
   break;
  }
 


 }
 return false;

}

int AddConnectModel(RESULT_MODEL **all, int N,int Ntot,SEQ *seq){
 puts("Start COmpare....");
 int Num=N;
 bool DisSim=false;
 bool FindSimilar=false;
 float d2=0;
 float d2cut=1.00;
 RESULT_MODEL *f1,*f2,*f3;
 unsigned int Nchk=0;

 RESULT_MODEL ***tmp_tbl;

 int *Ntbl;
 if((Ntbl=(int *)calloc(N,sizeof(int)))==NULL)
  return -1;
 if((tmp_tbl=(RESULT_MODEL ***)malloc(sizeof(RESULT_MODEL **)*N))==NULL)
  return -1;
 for(int i=0;i<N;i++)
  if((tmp_tbl[i]=(RESULT_MODEL **)malloc(sizeof(RESULT_MODEL *)*N))==NULL)
   return -1;

 int *NbestID,*CbestID;
 if((NbestID=(int *)malloc(sizeof(int)*N))==NULL)
  return -1;
 if((CbestID=(int*)malloc(sizeof(int)*N))==NULL)
  return -1;



 //select two fragments
 puts("Checking Connection....");
 #pragma omp parallel for schedule(dynamic,10)
 for(int fgid1=0;fgid1<Num;fgid1++){
  //for(int fgid2=fgid1+1;fgid2<Num;fgid2++){
  for(int fgid2=0;fgid2<Num;fgid2++){
	if(fgid1==fgid2)
	 continue;
	if(CheckGeo(all[fgid1],all[fgid2])==false)
	 continue;
	tmp_tbl[fgid1][Ntbl[fgid1]]=all[fgid2];
	Ntbl[fgid1]++;
	//tmp_tbl[fgid2][Ntbl[fgid2]]=all[fgid1];
	//Ntbl[fgid2]++;
  }
 }
 //Show..connections
 for(int fgid1=0;fgid1<Num;fgid1++){
  printf("#FG%d Ntbl= %d\n",fgid1,Ntbl[fgid1]);
  NbestID[fgid1]=CbestID[fgid1]=-1;
 }

 for(int fgid1=0;fgid1<Num;fgid1++){
  int Nfind=0;
  for(int fgid2=fgid1+1;fgid2<Num;fgid2++){
	int type=CheckOverlap(all[fgid1]->Start,all[fgid1]->End,all[fgid2]->Start,all[fgid2]->End);
	//ignore overlap region
	if(type==-1)
	 continue;
	
	if(FindConnects(all[fgid1],tmp_tbl[fgid1],Ntbl[fgid1],all[fgid2],tmp_tbl[fgid2],Ntbl[fgid2],seq,false)==true){
		if(type==0){
		 //Extend C-ter of fgid1
		 if(CbestID[fgid1]==-1){
		  CbestID[fgid1]=fgid2;
		 }else if(all[CbestID[fgid1]]->Score < all[fgid2]->Score){
		  CbestID[fgid1]=fgid2;
		  printf("Add %d-%d %d-%d (%d)\n",all[fgid1]->Start,all[fgid1]->End,all[fgid2]->Start,all[fgid2]->End,fgid2);
		 }
		 //Extend N-ter of fgid2
		 if(NbestID[fgid2]==-1){
		  NbestID[fgid2]=fgid1;
		 }else if(all[NbestID[fgid2]]->Score < all[fgid1]->Score){
		  NbestID[fgid2]=fgid1;
		 }
		}
		if(type==1){
		 //Extend N-ter of fgid1
		 if(NbestID[fgid1]==-1)
		  NbestID[fgid1]=fgid2;
		 else if(all[NbestID[fgid1]]->Score < all[fgid2]->Score){
		  NbestID[fgid1]=fgid2;
		 }
		 //Extend C-ter of fgid2
		 if(CbestID[fgid2]==-1){
		  CbestID[fgid2]=fgid1;
		 }else if(all[CbestID[fgid2]]->Score < all[fgid1]->Score){
		  CbestID[fgid2]=fgid1;
		 }
		}
	}
  }
 }

 //Remove redundant pairs
 for(int fgid1=0;fgid1<Num;fgid1++){
  int Nfgid = NbestID[fgid1];
  int Cfgid = CbestID[fgid1];
  if(Nfgid!=-1){
   if(CbestID[Nfgid]==fgid1){
    printf("Duplicated!!\n");
    continue;
   }

   printf("%d Nbest= %d %f ",fgid1,Nfgid,all[Nfgid]->Score);
   printf("%d-%d %d-%d\n",all[fgid1]->Start,all[fgid1]->End,all[Nfgid]->Start,all[Nfgid]->End);
   
   FindConnects(all[fgid1],tmp_tbl[fgid1],Ntbl[fgid1],all[Nfgid],tmp_tbl[Nfgid],Ntbl[Nfgid],seq,true);
  }
  if(Cfgid!=-1){
   printf("%d Cbest= %d %f ",fgid1,Cfgid,all[Cfgid]->Score);
   printf("%d-%d %d-%d\n",all[fgid1]->Start,all[fgid1]->End,all[Cfgid]->Start,all[Cfgid]->End);
   FindConnects(all[fgid1],tmp_tbl[fgid1],Ntbl[fgid1],all[Cfgid],tmp_tbl[Cfgid],Ntbl[Cfgid],seq,true);
  }
 }
 return Num;

}





//Model (NODEs) results

int AddModel(RESULT_MODEL *in,RESULT_MODEL **all, int N){
 //puts("Start COmpare....");
 int Num=N;
 bool DisSim=false;
 bool FindSimilar=false;
 float d2=0;
 float d2cut=1.00;
 for(int i=0;i<Num;i++){
  if(in->Len != all[i]->Len)
   continue;
  if(in->Start != all[i]->Start)
   continue;
  if(in->End != all[i]->End)
   continue;

  	//Check positions
	DisSim=false;
	for(int p=0;p<in->Len;p++){
	 d2=(in->Ca[p].real_cd[0]-all[i]->Ca[p].real_cd[0])
	   *(in->Ca[p].real_cd[0]-all[i]->Ca[p].real_cd[0])
	   +(in->Ca[p].real_cd[1]-all[i]->Ca[p].real_cd[1])
	   *(in->Ca[p].real_cd[1]-all[i]->Ca[p].real_cd[1])
	   +(in->Ca[p].real_cd[2]-all[i]->Ca[p].real_cd[2])
	   *(in->Ca[p].real_cd[2]-all[i]->Ca[p].real_cd[2]);

	 	if(d2>d2cut){
		 DisSim=true;
		 break;
		}
	}
	if (DisSim==false)
	 FindSimilar=true;
	
	//Compare and Swap
	if(FindSimilar==true && in->Score > all[i]->Score){
	 printf("Swap.. %d\n",i);
  	 all[i]->Start = in->Start;
  	 all[i]->End = in->End;
	 all[i]->Score = in->Score;
	 all[i]->RawScore = in->RawScore;
	 all[i]->Rmsd = in->Rmsd;
	 all[i]->Len = in->Len;
	 	for(int p=0;p<in->Len;p++){
	 	 all[i]->Ca[p].real_cd[0]=in->Ca[p].real_cd[0];
	 	 all[i]->Ca[p].real_cd[1]=in->Ca[p].real_cd[1];
	 	 all[i]->Ca[p].real_cd[2]=in->Ca[p].real_cd[2];
	 	 all[i]->Ca[p].AAtype=in->Ca[p].AAtype;
		}
	}
  	if(FindSimilar==true)
   	 break;
 }
 int i=Num;
 if(FindSimilar==false){//Add New Data
  //puts("Add New Model");
  Num++;

  all[i]->Start = in->Start;
  all[i]->End = in->End;
  all[i]->Score = in->Score;
  all[i]->RawScore = in->RawScore;
  all[i]->Rmsd = in->Rmsd;
  all[i]->Len = in->Len;

	for(int p=0;p<in->Len;p++){
	  all[i]->Ca[p].real_cd[0]=in->Ca[p].real_cd[0];
	  all[i]->Ca[p].real_cd[1]=in->Ca[p].real_cd[1];
	  all[i]->Ca[p].real_cd[2]=in->Ca[p].real_cd[2];
	  all[i]->Ca[p].AAtype=in->Ca[p].AAtype;
	  all[i]->Ca[p].node_id=in->Ca[p].node_id;
	}
 }
 //printf("END Comp... %d\n",Num);
 return Num;

}


//Compare with AF2 model
bool CompareFragAmodelMore(SEQ_NODE *result,SEQFG *sfg,int Nsfg,NODE *mod,bool mode){
 int MaxInter=20;

 //Malloc MODEL
 RESULT_MODEL **model,tmp_model;
 int Nmodel=0;
 int MaxModel=2000;

 if((model=(RESULT_MODEL **)malloc(sizeof(RESULT_MODEL *)*MaxModel))==NULL)
  return false;
 for(int i=0;i<MaxModel;i++)
  if((model[i]=(RESULT_MODEL *)malloc(sizeof(RESULT_MODEL)))==NULL)
   return false;





 //select two positions
 for(int s1=0;s1<Nsfg;s1++){
  int posi1=sfg[s1].pos;
  int posi1ed=sfg[s1].pos+sfg[s1].l-1;
  int N1=result[s1].Nali;
  SEQ_NODE *res1=&result[s1];
  double cd1[100][3],cd2[100][3];
  double mtx[3][3],mov_com[3],mov_ref[3];
  double rmsd;
	//Fragment
 	for(int s2=s1+1;s2<Nsfg;s2++){
	 int posi2=sfg[s2].pos;
  	 int posi2ed=sfg[s2].pos+sfg[s2].l-1;
  	 int N2=result[s2].Nali;
  	 SEQ_NODE *res2=&result[s2];
	 
	 if(posi1<posi2 && posi1ed > posi2)
	  continue;
	 if(posi2<posi1 && posi2ed > posi1)
	  continue;
	 if(posi1ed<posi2 && posi2-posi1ed > MaxInter)
	  continue;
	 if(posi2ed<posi1 && posi1-posi2ed > MaxInter)
	  continue;
	 printf("Comp %d-%d %d-%d %d vs %d\n",posi1,posi1ed,posi2,posi2ed,N1,N2);
	 //init
	 Nmodel=0;
	 	//Copy coords
		int RefSt,RefEd;
		RefSt=posi1;
		RefEd=posi1;
		if(posi2<posi1)
		 RefSt=posi2;
		for(int i1 =0;i1<N1;i1++){
		 //Zscore cutoff
		 if(res1->score[i1] < cmd.zcut)
		  continue;
		for(int i2 =0;i2<N2;i2++){
		 if(res2->score[i2] < cmd.zcut)
		  continue;
		 int Ncd=0;
		 for(int j=0;j<sfg[s1].l;j++){
	  	  if(res1->ali[i1][j]==NULL)
	   	   continue;
		  if(posi1+j>RefEd)
		   RefEd=posi1+j+1;

	   	  cd1[Ncd][0]=res1->ali[i1][j]->real_cd[0];
	   	  cd1[Ncd][1]=res1->ali[i1][j]->real_cd[1];
	   	  cd1[Ncd][2]=res1->ali[i1][j]->real_cd[2];

	   	  cd2[Ncd][0]=mod[posi1+j].real_cd[0];
	   	  cd2[Ncd][1]=mod[posi1+j].real_cd[1];
	   	  cd2[Ncd][2]=mod[posi1+j].real_cd[2];

	  	  Ncd++;
		 }
		 for(int j=0;j<sfg[s2].l;j++){
	  	  if(res2->ali[i2][j]==NULL)
	   	   continue;
		  if(posi2+j>RefEd)
		   RefEd=posi2+j+1;

	   	  cd1[Ncd][0]=res2->ali[i2][j]->real_cd[0];
	   	  cd1[Ncd][1]=res2->ali[i2][j]->real_cd[1];
	   	  cd1[Ncd][2]=res2->ali[i2][j]->real_cd[2];

	   	  cd2[Ncd][0]=mod[posi2+j].real_cd[0];
	   	  cd2[Ncd][1]=mod[posi2+j].real_cd[1];
	   	  cd2[Ncd][2]=mod[posi2+j].real_cd[2];

	  	  Ncd++;
		 }
		 //Check RMSD
		 fast_rmsd(cd1,cd2,Ncd,&rmsd);
		/*
		 res->rmsd[i]=rmsd;
	 	 if(!isfinite(rmsd))
	  	  res->rmsd[i]=100;
		*/
		 if(rmsd>=cmd.RMSD)
		  continue;
		
	 	  printf("RMSD %d-%d(%d) vs %d-%d(%d) RMSD= %.3f %d-%d\n",posi1,posi1ed,i1,posi2,posi2ed,i2,rmsd,RefSt,RefEd);
		 //Copy Coords again...then compute r&t
		 Ncd=0;
		 for(int j=0;j<sfg[s1].l;j++){
	  	  if(res1->ali[i1][j]==NULL)
	   	   continue;
	   	  cd1[Ncd][0]=res1->ali[i1][j]->real_cd[0];
	   	  cd1[Ncd][1]=res1->ali[i1][j]->real_cd[1];
	   	  cd1[Ncd][2]=res1->ali[i1][j]->real_cd[2];
	   	  cd2[Ncd][0]=mod[posi1+j].real_cd[0];
	   	  cd2[Ncd][1]=mod[posi1+j].real_cd[1];
	   	  cd2[Ncd][2]=mod[posi1+j].real_cd[2];
	  	  Ncd++;
		 }
		 for(int j=0;j<sfg[s2].l;j++){
	  	  if(res2->ali[i2][j]==NULL)
	   	   continue;
	   	  cd1[Ncd][0]=res2->ali[i2][j]->real_cd[0];
	   	  cd1[Ncd][1]=res2->ali[i2][j]->real_cd[1];
	   	  cd1[Ncd][2]=res2->ali[i2][j]->real_cd[2];
	   	  cd2[Ncd][0]=mod[posi2+j].real_cd[0];
	   	  cd2[Ncd][1]=mod[posi2+j].real_cd[1];
	   	  cd2[Ncd][2]=mod[posi2+j].real_cd[2];
	  	  Ncd++;
		 }
		 calculate_rotation_rmsd(cd1,cd2,Ncd,mov_com,mov_ref,mtx,&rmsd);
	 	 //printf("##POS %d RAWSCO= %f ZSCORE= %f RMSD= %f Amodel\n",
		 //RefSt,res1->raw_score[i1]+res2->raw_score[i2],res1->score[i1]+res2->score[i2],rmsd);


		 tmp_model.Start=RefSt;
		 tmp_model.End=RefEd;
		 tmp_model.Score=res1->score[i1]+res2->score[i2];
		 tmp_model.RawScore=res1->raw_score[i1]+res2->raw_score[i2];
		 tmp_model.Rmsd=rmsd;
		 tmp_model.Len=RefEd-RefSt;
		 //move rotate and extend
		 float rotated[3],tmp[3];
		 for(int j=RefSt;j<RefEd;j++){
		  //mov -> o
		  tmp[0]=mod[j].real_cd[0]-mov_com[0];
		  tmp[1]=mod[j].real_cd[1]-mov_com[1];
		  tmp[2]=mod[j].real_cd[2]-mov_com[2];
			for(int i1=0;i1<3;i1++){
				rotated[i1]=0.0;
				for(int i2=0;i2<3;i2++)
				 rotated[i1] +=mtx[i1][i2]*tmp[i2];
			}
		  //o -> mov -> ref
		  tmp[0]=rotated[0]+mov_com[0]+mov_ref[0];
		  tmp[1]=rotated[1]+mov_com[1]+mov_ref[1];
		  tmp[2]=rotated[2]+mov_com[2]+mov_ref[2];
		  //printf("ATOM  %5d  CA  %3s%6d    ",j+1,RES_NAMES[mod[j].AAtype],j+1);
  		  //printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,res1->score[i1]+res2->score[i2]);
		  tmp_model.Ca[j-RefSt].real_cd[0]=tmp[0];
		  tmp_model.Ca[j-RefSt].real_cd[1]=tmp[1];
		  tmp_model.Ca[j-RefSt].real_cd[2]=tmp[2];
		  tmp_model.Ca[j-RefSt].AAtype=mod[j].AAtype;
		 }
 	 	 //puts("ENDMDL");
		 //ADD tmp_model
		 Nmodel=AddModel(&tmp_model,model,Nmodel);
		 printf("#Nmodel= %d\n",Nmodel);
		}}
	 
		//ShowModel

		for(int m=0;m<Nmodel;m++){
	 	 printf("##POS %d RAWSCO= %f ZSCORE= %f RMSD= %f Amodel\n",
		 model[m]->Start,model[m]->RawScore,model[m]->Score,model[m]->Rmsd);
	 	 puts("MODEL");
			for(int p=0;p<model[m]->Len;p++){
		 printf("ATOM  %5d  CA  %3s%6d    ",model[m]->Start+1+p,RES_NAMES[model[m]->Ca[p].AAtype],model[m]->Start+1+p);
  		 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",model[m]->Ca[p].real_cd[0],model[m]->Ca[p].real_cd[1],model[m]->Ca[p].real_cd[2],1.0,model[m]->Score);
			}
 	 	 puts("ENDMDL");
		}

	}
 }
}





bool MakeDummyFrag(SEQ_NODE *new,SEQ_NODE *result,SEQFG *sfg,int Nsfg,bool mode){

 int Ncd=0;
 //Remove redundant fragments
 for(int s1=0;s1<Nsfg;s1++){
  SEQ_NODE *res1=&result[s1];
  int N1=res1->Nali;
  //double cd1[100][3],cd2[100][3];
  if(N1==0)
   continue;
	for(int i1=0;i1<N1;i1++){
	 if(res1->rmsd[i1]==100)
	  continue;
		for(int s2=s1+1;s2<Nsfg;s2++){
		 SEQ_NODE *res2=&result[s2];
		 int N2=res2->Nali;
		 if(N2==0)
		  continue;
			for(int i2=0;i2<N2;i2++){
		 	 if(res2->rmsd[i2]==100)
		 	  continue;
				//Comp two fragment
			 bool flag=false;
	 			for(int j=0;j<sfg[s1].l && j<sfg[s2].l;j++){
				 if(res1->ali[i1][j]!=res2->ali[i2][j])
				  break;
	  			 if(res1->ali[i1][j]==NULL)
	  			  break;
	  			 if(res2->ali[i2][j]==NULL)
	  			  break;
				 flag=true;
				}
			 if(flag==true){
			  //printf("Same %d %d\n",s1,s2);
			  res2->rmsd[i2]=100;
			 }
			}
		}

	}
 }
 int Nnr=0;
 for(int s1=0;s1<Nsfg;s1++){
  SEQ_NODE *res1=&result[s1];
  int N1=res1->Nali;
  //double cd1[100][3],cd2[100][3];
  if(N1==0)
   continue;
	for(int i1=0;i1<N1;i1++){
	 if(res1->rmsd[i1]==100)
	  continue;
	 Ncd=0;
		for(int j=0;j<sfg[s1].l;j++){
		 if(res1->ali[i1][j]==NULL)
	  	  break;
		 Ncd++;
		}
	 if(Ncd<5){
	  res1->rmsd[i1]=100;
	  continue;
	 }
	 Nnr++;
	}
 }
 printf("#Uniq fragment= %d * %d = %d\n",Nnr,Nsfg, Nnr*Nsfg);
 //Input to a new results
 //malloc res
 for(int i=0;i<Nsfg;i++){
  if((new[i].score=(float*)malloc(sizeof(float)*Nnr))==NULL)
   return true;
  if((new[i].raw_score=(float*)malloc(sizeof(float)*Nnr))==NULL)
   return true;
  if((new[i].rmsd=(float*)malloc(sizeof(float)*Nnr))==NULL)
   return true;
  if((new[i].ali=(NODE ***)malloc(sizeof(NODE **)*Nnr))==NULL)
   return true;
  for(int j=0;j<Nnr;j++)
   if((new[i].ali[j]=(NODE **)malloc(sizeof(NODE *)*100))==NULL)
    return true;

  SEQ_NODE *new_res=&new[i];
  new[i].Nali = Nnr;
  int N=0;
  	for(int s1=0;s1<Nsfg;s1++){
  	 SEQ_NODE *res1=&result[s1];
  	 int N1=res1->Nali;
  	 if(N1==0)
  	  continue;
		for(int i1=0;i1<N1;i1++){
	 	 if(res1->rmsd[i1]==100)
	 	  continue;
		 new_res->raw_score[N]=0.0001;
		 new_res->score[N]=0.0001;
		 new_res->rmsd[N]=0.00;
		 for(int j=0;j<sfg[s1].l;j++){
		  new_res->ali[N][j]=res1->ali[i1][j];
		  if(res1->ali[i1][j]==NULL)
	  	   break;
		 }
	 	 N++;
	 	}
	}
  //printf("Nali= %d\n",new[i].Nali); 
 }
 //Input Data
 
}

float CheckConnect(NODE *Ca1,NODE *Ca2, MRC *mrc,float Pcut){

 //vec Ca1->Ca2
 float vec[3],pos[3];
 int Ipos[3];
 int xydim=mrc->xdim*mrc->ydim;
 int xdim = mrc->xdim;
 int x,y,z,ind;
 vec[0]=Ca2->real_cd[0]-Ca1->real_cd[0];
 vec[1]=Ca2->real_cd[1]-Ca1->real_cd[1];
 vec[2]=Ca2->real_cd[2]-Ca1->real_cd[2];

 if(Ca1 == Ca2){ //Single point
  pos[0]=Ca1->real_cd[0];
  pos[1]=Ca1->real_cd[1];
  pos[2]=Ca1->real_cd[2];

  x=(int)round((pos[0]-mrc->orgxyz[0])/mrc->widthx);
  y=(int)round((pos[1]-mrc->orgxyz[1])/mrc->widthx);
  z=(int)round((pos[2]-mrc->orgxyz[2])/mrc->widthx);
  ind=xydim*z+xdim*y+x;

  printf("%d %d %d dens= %f\n",x,y,z,mrc->dens[ind]);
  if(mrc->dens[ind]< Pcut)
   return true;//No-connection

  return false;
 }
 //Sample 10 points
 float MinD = 1000000.000;
 for(int i=1;i<10;i++){ //1-9
  pos[0]=Ca1->real_cd[0] + 0.10*(i)*vec[0];
  pos[1]=Ca1->real_cd[1] + 0.10*(i)*vec[1];
  pos[2]=Ca1->real_cd[2] + 0.10*(i)*vec[2];

  //printf("ATOM  %5d  CA  %3s%6d    ",i,"ALA",i);
  //printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",pos[0],pos[1],pos[2],1.0,1.0);

  //check density value
  //g->node[i].real_cd[0]=pt->cd[i][0]*mrc->widthx+mrc->orgxyz[0];
  x=(int)round((pos[0]-mrc->orgxyz[0])/mrc->widthx);
  y=(int)round((pos[1]-mrc->orgxyz[1])/mrc->widthx);
  z=(int)round((pos[2]-mrc->orgxyz[2])/mrc->widthx);
  
  float SumD=0.00;
  ind=xydim*(z)+xdim*(y)+(x);
  SumD=mrc->dens[ind];
  //SumD = Best27;
  //printf("SumD= %f %d %d %d\n",SumD,x,y,z);
  if(SumD < MinD)
   MinD = SumD;
 }
 //printf("MinD= %f\n",MinD);
 //return MinD/27.00;
 return MinD;
}

int CountModels(PDB *p,SEQ *seq, int flen){
 int cnt=0;
 int tmpStart=-1;
 for(int i=0;i<p->NumOfAtom - flen ;i++){
  int ResNum = p->ResNum[i];//Start from 1
  int order = p->AtomOnRes[i];
 //Check Sequence
  if(p->TypeResId[order]!=seq->seq_code[ResNum-1])
   continue;
	//Check fragment
	int Ncd=0;
	for(int j=0;j<flen && i+j < p->NumOfAtom ;j++){
	 //Check Residue Number
	 if(ResNum + j != p->ResNum[i+j])
	  continue;
	 Ncd++;
	}
	if(Ncd != flen)
	 continue;
	tmpStart=ResNum-1;
	cnt++;
	//Check Homo-oligomer fragments
	for(int pos=0;pos<seq->len-flen;pos++){
		if(pos==tmpStart)
		 continue;
		bool homo_check=true;
		for(int k=0;k<flen;k++){
		 if(seq->seq_code[pos+k]!=seq->seq_code[tmpStart+k]){
		  homo_check=false;
		  break;
		 }
		}
		if(homo_check==true){
		 //printf("Found Homo %d %d\n",tmpStart,pos);
		 cnt++;
		}
	}
 }

 printf("##Total Models= %d\n",cnt);
 return cnt;
}

bool CheckContact(RESULT_MODEL *f1,RESULT_MODEL *f2,short int *tbl,int Ntbl){
  //int mode = CheckOverlap(f1->Start,f1->End,f2->Start,f2->End);
 float dist2=0.0;
 float FarDist = 8.0*8.0;

 //puts("Overlap");
 for(int a=0;a<f1->Len;a++){
  int pos1 = f1->Start+a;
  for(int b=0;b<f2->Len;b++){
   int pos2 = f2->Start+b;
   dist2=Dist2(&(f1->Ca[a]),&(f2->Ca[b]));
   if(dist2>FarDist){
	printf("pos0 %d %d = %d \n",pos1,pos2,Ntbl);
	if(tbl[pos1*Ntbl+pos2]==1)
	 return true;
	if(tbl[pos2*Ntbl+pos1]==1)
	 return true;

   }
  }
 }
 puts("DONE tbl");
 return false;
}
