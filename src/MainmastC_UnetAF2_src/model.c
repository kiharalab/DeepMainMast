#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "struct.h"
#include "mrc.h"
#include "thread.h"
#include "model.h"
#include "rmsd.h"

extern CMD cmd;

float Dist2(NODE *Ca1,NODE *Ca2){
	float d2=(Ca1->real_cd[0]-Ca2->real_cd[0])
		*(Ca1->real_cd[0]-Ca2->real_cd[0])
		+(Ca1->real_cd[1]-Ca2->real_cd[1])
		*(Ca1->real_cd[1]-Ca2->real_cd[1])
		+(Ca1->real_cd[2]-Ca2->real_cd[2])
		*(Ca1->real_cd[2]-Ca2->real_cd[2]);

 return d2;
}



bool ReadNodeCsv(char *filename,NODE_CSV *out){
 int Nnode=0;
 int Npath=0;
 int MaxChar=10000*5;
 printf("#Loading..%s\n",filename);
 FILE *fp = fopen(filename,"r");
 if (!fp)
  return true;
 
 char buf[MaxChar];

 int MaxPath=0;
 int Np;
 int Lseq=0;
 while(fgets(buf,MaxChar,fp)){
	char *value = strtok(buf,", ");
	if(!strcmp(value,"NODE"))
	 Nnode++;
	if(!strcmp(value,"SCO"))
	 Lseq++;
	if(!strcmp(value,"PATH")){
	 Npath++;
	 Np=0;
		while(value){
        	 //printf("%s\n",value);
        	 value = strtok(NULL,", ");
		 Np++;
        	}
	 if(Np>MaxPath)
	  MaxPath=Np;
	}
 }
 printf("##Nnode=%d Npath=%d Max=%d Lseq=%d\n",Nnode,Npath,MaxPath,Lseq);
 fclose(fp);

 out->Nnode=Nnode;
 out->Npath=Npath;
 out->LongestPath=MaxPath;
 //Malloc
 if((out->Ca=(NODE *)malloc(sizeof(NODE)*Nnode))==NULL)
  return true;
 if((out->Path=(int **)malloc(sizeof(int*)*Npath*2))==NULL)
  return true;
 for(int i=0;i<Npath*2;i++)
  if((out->Path[i]=(int *)malloc(sizeof(int)*MaxPath+1))==NULL)
   return true;
 if((out->Smtx=(int **)malloc(sizeof(int *)*Lseq))==NULL)
  return true;
 for(int i=0;i<Lseq;i++)
  if((out->Smtx[i]=(int *)malloc(sizeof(int)*(Nnode+1)))==NULL)
   return true;
 puts("##Fin Malloc");


 fp = fopen(filename,"r");
 if (!fp)
  return true;

 while(fgets(buf,MaxChar,fp)){
        char *value = strtok(buf,", ");
	//NODE
        if(!strncmp(buf,"NODE",4)){
	 //puts("NODE");
	 value = strtok(NULL,", ");
	 int id = atoi(value);
	 value = strtok(NULL,", ");
	 float x = atof(value);
	 value = strtok(NULL,", ");
	 float y = atof(value);
	 value = strtok(NULL,", ");
	 float z = atof(value);
	 out->Ca[id].real_cd[0]=x;
	 out->Ca[id].real_cd[1]=y;
	 out->Ca[id].real_cd[2]=z;
	}
	//PATH
        if(!strncmp(buf,"PATH",4)){
	 //puts("PATH");
         Np=0;
	 value = strtok(NULL,", ");
	 int id = atoi(value);
	 value = strtok(NULL,", ");
                while(value){
		 out->Path[id][Np]=atoi(value);
                 value = strtok(NULL,", ");
                 Np++;
		 //printf("%d %d/%d\n",Np,out->Path[id][Np-1],MaxPath);
                }
		out->Path[id][Np]=-1;//Terminal
		out->Lpath[id]=Np;
        }

	//SCO
	if(!strncmp(buf,"SCO",3)){
	 //puts("SCO");
         Np=0;
	 value = strtok(NULL,", ");
	 int id = atoi(value);
	 value = strtok(NULL,", ");

                while(value){
		 out->Smtx[id][Np]=atoi(value);
                 value = strtok(NULL,", ");
                 Np++;
                }
        }
 }
 fclose(fp);
/*
 //Check splitted path d>10
 for(int i=0;i<out->Npath;i++){
  printf("##Check %d\n",i);
	for(int j=0;j+1<out->Lpath[i];j++){
	 int i1=out->Path[i][j];
	 int i2=out->Path[i][j+1];
	 float d2=Dist2(&(out->Ca[i1]),&(out->Ca[i2]));
	 if(d2>100.0){//d>10A
	  printf("Split Path %d %d %f\n",i,out->Npath,sqrt(d2));
		int p=0;
		for(p=0;j+p<out->Lpath[i];p++){
		 //printf("p=%d /%d\n",p,MaxPath);
		 out->Path[out->Npath][p]=out->Path[i][j+p];
		 //out->Path[out->Npath][p]=0;
		}
		//out->Lpath[out->Npath]=p;
	  out->Npath++;
	  break;
	 }
	}
 }
*/
 return false;
}

bool GetFragFromNODE(NODE_CSV *p,int pos,int L,double out[100][3]){
 //Check sequence position
 int res0=p->Ca[pos].seq_id;
 for(int i=0;i<L;i++){
  if(pos+i >=p->Nnode)
   return true;
  int res1=p->Ca[pos+i].seq_id;
  if(res1-res0 != i)//seq_id check
   return true;
  if(p->Ca[pos].cid != p->Ca[pos+i].cid) //Same chain
   return true;
  out[i][0]=p->Ca[pos+i].real_cd[0];
  out[i][1]=p->Ca[pos+i].real_cd[1];
  out[i][2]=p->Ca[pos+i].real_cd[2];
 }

 return false;
}

int GetChainFragFromNODE(NODE_CSV *p,int pos,double out[2000][3],NODE *node[2000]){
 int cid=p->Ca[pos].cid;
 //printf("pos %d cid %d\n",pos,cid);
 //init
 int Nchain=0;
 for(int i=0;i<p->Nnode;i++){
  if(cid == p->Ca[i].cid){

   out[Nchain][0]=p->Ca[i].real_cd[0];
   out[Nchain][1]=p->Ca[i].real_cd[1];
   out[Nchain][2]=p->Ca[i].real_cd[2];
   node[Nchain]=&(p->Ca[i]);
   Nchain++;
  }

 }
 return Nchain;
}


bool GetFragFromPDB(PDB *p,int pos,int L,double out[100][3]){
 //Check sequence position
 int res0=p->ResNum[pos];
 for(int i=0;i<L;i++){
  if(pos+i >=p->NumOfAtom)
   return true;
  int res1=p->ResNum[pos+i];
  if(res1-res0 != i)//seq_id
   return true;
  out[i][0]=p->xyz[pos+i][0];
  out[i][1]=p->xyz[pos+i][1];
  out[i][2]=p->xyz[pos+i][2];
 }

 return false;
}

bool GetFragFromPATH(NODE_CSV *in,int path,int pos,int L,double out[100][3]){

 int L1=0;
 for(int pos1=pos;pos1<pos+L;pos1++){
  int nid=in->Path[path][pos1];
  if(in->Path[path][pos1]==-1)
   return true;
  out[L1][0]=in->Ca[nid].real_cd[0];
  out[L1][1]=in->Ca[nid].real_cd[1];
  out[L1][2]=in->Ca[nid].real_cd[2];
  L1++;
 }
 if(L1!=L)
  return true;
 return false;
}

bool GetFragFromTBL(NODE_CSV *in,int **tbl,int pos,int L,double out[100][3]){

 for(int i=0;i<L;i++){
  int nid=tbl[pos][i];
  out[i][0]=in->Ca[nid].real_cd[0];
  out[i][1]=in->Ca[nid].real_cd[1];
  out[i][2]=in->Ca[nid].real_cd[2];
 }
 return false;
}

bool OvCheck(int **tbl,int pos1,int pos2,int L){
 for(int i=0;i<L;i++){
  int nid1=tbl[pos1][i];
  for(int j=0;j<L;j++){
   int nid2=tbl[pos2][j];
   if(nid1==nid2)
    return true;
  }
 }
 return false;
}

bool ChkContPairs(int *tbl,int N,int id1,int id2){
 for(int i=0;i<N;i++){
  if(tbl[2*i]==id1 && tbl[2*i+1]==id2)
   return true;
  if(tbl[2*i]==id2 && tbl[2*i+1]==id1)
   return true;
 }
 return false;
}

float lb3d(double cd11[100][3],double cd12[100][3],double cd21[100][3],double cd22[100][3],int L){
 double g11[3]={0.0,0.0,0.0};
 double g12[3]={0.0,0.0,0.0};
 double g21[3]={0.0,0.0,0.0};
 double g22[3]={0.0,0.0,0.0};
 for(int i=0;i<L;i++){
  g11[0]+=cd11[i][0]; g11[1]+=cd11[i][1]; g11[2]+=cd11[i][2];
  g12[0]+=cd12[i][0]; g12[1]+=cd12[i][1]; g12[2]+=cd12[i][2];
  g21[0]+=cd21[i][0]; g21[1]+=cd21[i][1]; g21[2]+=cd21[i][2];
  g22[0]+=cd22[i][0]; g22[1]+=cd22[i][1]; g22[2]+=cd22[i][2];
 }
 for(int i=0;i<3;i++){
  g11[i]/=(double)(L);
  g12[i]/=(double)(L);
  g21[i]/=(double)(L);
  g22[i]/=(double)(L);
 }
 double tmp[3];
 tmp[0]=g11[0]-g12[0]; tmp[1]=g11[1]-g12[1]; tmp[2]=g11[2]-g12[2];
 double Fs = sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
 tmp[0]=g21[0]-g22[0]; tmp[1]=g21[1]-g22[1]; tmp[2]=g21[2]-g22[2];
 double Ft = sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
 //Too far?
/*
 if(Fs > 20.0)
  return 100.0;
 if(Ft > 20.0)
  return 100.0;
*/
 if(Fs-Ft > 0)
  return Fs-Ft;
 else
  return Ft-Fs;
}

//impose cd2 ->cd2 and then rotate&shift cd3
bool move_by_mtx(double cd1[2000][3],double cd2[2000][3],int Ncd,double cd3[2000][3],
		 double out[2000][3],int Ncd_out){ 
 float rotated[3],tmp[3];
 double mtx[3][3],mov_com[3],mov_ref[3];
 double rmsd;
 double CD1[2000][3],CD2[2000][3];
 //Before rotation, coodinates need to be copied
 
 for(int i=0;i<Ncd;i++){
  CD1[i][0]=cd1[i][0]; CD1[i][1]=cd1[i][1]; CD1[i][2]=cd1[i][2];
  CD2[i][0]=cd2[i][0]; CD2[i][1]=cd2[i][1]; CD2[i][2]=cd2[i][2];
 }


 calculate_rotation_rmsd(CD1,CD2,Ncd,mov_com,mov_ref,mtx,&rmsd);

 //printf("MODEL\n");
 for(int i=0;i<Ncd_out;i++){
 	//mov -> o
	//tmp[0]=cd2[i][0]-mov_com[0];
	//tmp[1]=cd2[i][1]-mov_com[1];
	//tmp[2]=cd2[i][2]-mov_com[2];
	tmp[0]=cd3[i][0]-mov_com[0];
	tmp[1]=cd3[i][1]-mov_com[1];
	tmp[2]=cd3[i][2]-mov_com[2];
	for(int i1=0;i1<3;i1++){
	 rotated[i1]=0.0;
	 for(int i2=0;i2<3;i2++)
	  rotated[i1] +=mtx[i1][i2]*tmp[i2];
	}
 	//o -> mov -> ref
 	tmp[0]=rotated[0]+mov_com[0]+mov_ref[0];
 	tmp[1]=rotated[1]+mov_com[1]+mov_ref[1];
 	tmp[2]=rotated[2]+mov_com[2]+mov_ref[2];
 	//printf("ATOM  %5d  CA  %3s%6d    ",i,"ALA",i);
  	//printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,1.0);
	out[i][0]=tmp[0];
	out[i][1]=tmp[1];
	out[i][2]=tmp[2];
 }
 //printf("TER\nENDMDL\n");
 return false;
}


bool Homo_check(SEQ *seq,int pos,int N,RESULT_MODEL *in,RESULT_MODEL *out){
 int St=in->Start;
 int Ed=in->End;
 int Len=Ed-St;
 bool homo=true;
 for(int s=0;s<Len;s++){
  if(pos+s>=seq->len)
   return false;
  if(seq->seq_code[pos+s]!=seq->seq_code[St+s])//Start from 0
   return false;

 }
 
 //printf("Found Homo! %d -> %d\n",St,pos);
 int seq_dev;
 for(int j=0;j<in->Len;j++){
  out->Ca[j].real_cd[0]=in->Ca[j].real_cd[0];
  out->Ca[j].real_cd[1]=in->Ca[j].real_cd[1];
  out->Ca[j].real_cd[2]=in->Ca[j].real_cd[2];
  out->Ca[j].node_id=in->Ca[j].node_id;
  seq_dev=in->Ca[j].seq_id-in->Ca[0].seq_id;
  out->Ca[j].seq_id = (pos+1)+seq_dev;//Start from 1
 }
 out->Start=pos+1;
 out->End=pos+Len;
 out->Score=in->Score;//Not yet
 out->RawScore=0;//not yet
 out->Rmsd=in->Rmsd;
 out->Len=in->Len;//Number of Coords;
 return true;
}

int GetNrPath(NODE_CSV *in,int Fl,int **tbl){
 int Ntbl=0;
 int tmp_path[100];
 for(int pid=0;pid<in->Npath;pid++){
	for(int pos=0;pos<in->LongestPath;pos++){
	 if(in->Path[pid][pos]==-1)
	  break;
		bool chk=false;
	 	for(int j=0;j<Fl;j++){
		 if(in->Path[pid][j+pos]==-1){
		  chk=true;
		  break;
		 }
		 tmp_path[j]=in->Path[pid][j+pos];
		}
		if(chk==true)
		 break;
		//Add data to tbl
		chk=false;
		for(int i=0;i<Ntbl;i++){
		 bool same=true;
		 for(int j=0;j<Fl;j++){
		  if(tbl[i][j]!=tmp_path[j]){
		   same=false;
		   break;
		  }
		 }
		 if(same==true){
		  chk=true;
		  break;
		 }
		}
		//printf("TEST %d\n",Ntbl);
		if(chk==true)//Already used
		 continue;
		//printf("##NR add %d\n",Ntbl);
		if((tbl[Ntbl]=(int *)malloc(sizeof(int)*Fl))==NULL) return 0;
		for(int j=0;j<Fl;j++)
                 tbl[Ntbl][j]=tmp_path[j];
		Ntbl++;
		//Reverse
		if((tbl[Ntbl]=(int *)malloc(sizeof(int)*Fl))==NULL) return 0;
		for(int j=0;j<Fl;j++)
                 tbl[Ntbl][j]=tmp_path[Fl-1-j];
		Ntbl++;
	}
 }
 printf("##Ntbl= %d\n",Ntbl);
 return Ntbl;

}

float QScoring(RESULT_MODEL *mod,NODE_CSV *in){
 float sco=0.0;

 for(int i=0;i<mod->Len;i++){
  int node_id=mod->Ca[i].node_id;
  int seq_id=mod->Ca[i].seq_id;
	//printf("%d %d\n",seq_id,node_id);
  sco += in->Smtx[seq_id][node_id];
 }

 return sco;
}

int FindAlignment(int ali[1000][2],NODE_CSV *in,double cd[1000][3],int N,float r){
 int Nali=0;
 float r2=r*r;
 float d2;
 //int ali[1000][2];
 float dtbl[2000];
 for(int i=0;i<N;i++){
  float MinD=r2;
  int MinID=-1;
	for(int j=0;j<in->Nnode;j++){

	 d2=(cd[i][0]-in->Ca[j].real_cd[0])*(cd[i][0]-in->Ca[j].real_cd[0])
	   +(cd[i][1]-in->Ca[j].real_cd[1])*(cd[i][1]-in->Ca[j].real_cd[1])
	   +(cd[i][2]-in->Ca[j].real_cd[2])*(cd[i][2]-in->Ca[j].real_cd[2]);
	 if(d2<MinD){
	  MinD=d2;
	  MinID=j;
	 }
	}
	if(MinID>0){
	 dtbl[Nali]=MinD;
	 ali[Nali][0]=MinID;//node
	 ali[Nali][1]=i;//PDB position
	 Nali++;
	}
 }
 //No overlap
 for(int i=0;i<Nali-1;i++){
  for(int j=i+1;j<Nali;j++){
   if(ali[i][0]==ali[j][0]){
    if(dtbl[i]<dtbl[j]){
     ali[j][0]=-1;
    }else{
     ali[i][0]=-1;
     break;
    }
   }
 }}

/*
 printf("MODEL\n");
 for(int i=0;i<Nali;i++){
  if(ali[i][0]==-1)
   continue;
  int pos=ali[i][1];
  printf("ATOM  %5d  CA  %3s%6d    ",i+1,RES_NAMES[0],i+1);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",cd[pos][0],cd[pos][1],cd[pos][2],1.0,1.0);
 }
 printf("ENDMDL\n");
*/


 //Remove Isolated position
 int gp[2000];
 int group;
 //init
 for(int i=0;i<Nali;i++)
  gp[i]=i;

 for(int i=0;i<Nali;i++){
  if(ali[i][0]==-1)
   continue;
  group=gp[i];
  for(int j=i+1;j<Nali;j++){
   if(ali[j][0]==-1)
    continue;
   int p1=ali[i][1];//PDB cd
   int p2=ali[j][1];//PDB cd
  d2=(cd[p1][0]-cd[p2][0])*(cd[p1][0]-cd[p2][0])
    +(cd[p1][1]-cd[p2][1])*(cd[p1][1]-cd[p2][1])
    +(cd[p1][2]-cd[p2][2])*(cd[p1][2]-cd[p2][2]);
  //close <8A
  if(d2<8.0*8.0){
   int renew=gp[j];
   for(int k=0;k<Nali;k++)
    if(gp[k]==renew)
     gp[k]=group;
  }
 }}
 int cnt[2000]={0};
 for(int i=0;i<Nali;i++){
  //printf("%d,",gp[i]);
  cnt[gp[i]]++;
 }
 //puts("");
 int MaxGp=-1;
 int MaxCnt=0;
 for(int i=0;i<Nali;i++){
  //printf("%d %d\n",i,cnt[i]);
  if(cnt[i]>MaxCnt){
   MaxGp=i;
   MaxCnt=cnt[i];
  }
 }
 //printf("MaxGp %d\n",MaxGp);
 //Filter
 int Nali2=0;
 for(int i=0;i<Nali;i++){
  if(gp[i]==MaxGp){
	if(i!=Nali2){
  	 ali[Nali2][0]=ali[i][0];
  	 ali[Nali2][1]=ali[i][1];
	}
   Nali2++;
  }
 }

 //printf("Nali= %d\n",Nali);
 return Nali2;

 //return Nali;
}

int SetUpCd(int ali[2000][2],int Nali,NODE_CSV *in,PDB *p,double cd1[2000][3],double cd2[2000][3]){
 int L=0;
 for(int i=0;i<Nali;i++){
  int node_id=ali[i][0];
  int res_pos=ali[i][1];
  if(node_id==-1 || res_pos==-1)
   continue;
  cd1[L][0]=in->Ca[node_id].real_cd[0];
  cd1[L][1]=in->Ca[node_id].real_cd[1];
  cd1[L][2]=in->Ca[node_id].real_cd[2];

  cd2[L][0]=p->xyz[res_pos][0];
  cd2[L][1]=p->xyz[res_pos][1];
  cd2[L][2]=p->xyz[res_pos][2];

  L++;
 }
 return L;
}

int SetUpCd2(int ali[2000][2],int Nali,NODE_CSV *in,NODE *nodes[2000],double cd1[2000][3],double cd2[2000][3]){
 int L=0;
 for(int i=0;i<Nali;i++){
  int node_id=ali[i][0];
  int pos=ali[i][1];
  if(node_id==-1 || pos==-1)
   continue;
  cd1[L][0]=in->Ca[node_id].real_cd[0];
  cd1[L][1]=in->Ca[node_id].real_cd[1];
  cd1[L][2]=in->Ca[node_id].real_cd[2];
/*
  cd2[L][0]=p->xyz[res_pos][0];
  cd2[L][1]=p->xyz[res_pos][1];
  cd2[L][2]=p->xyz[res_pos][2];
*/
  cd2[L][0]=nodes[pos]->real_cd[0];
  cd2[L][1]=nodes[pos]->real_cd[1];
  cd2[L][2]=nodes[pos]->real_cd[2];
  
  L++;
 }
 return L;
}

//Extend loop and NC-terminal
int DetectFittedRegions(PDB *p,double cd[2000][3],int ali[2000][2],int Nali,RESULT_MODEL **model,int Nmodel){

 RESULT_MODEL mod;
 //printf("#Nali= %d\n",Nali);
 int Filled[2000],Nfill=0;
 int Pass[2000],Npass=0;

 //Filter
 for(int i=0;i<Nali;i++){
  int pos1=ali[i][1];
  int Lfrag=0;
	//To N-ter
	for(int j=i-1;j>=0;j--){
	 if((i-j)==(pos1-ali[j][1]))
	  Lfrag++;
	 else
	  break;
	}
	//To C-ter
	for(int j=i+1;j<Nali;j++){
	 if((j-i)==(ali[j][1]-pos1))
	  Lfrag++;
	 else
	  break;
	}
  if(Lfrag>=5){
   Pass[Npass++]=pos1;
  }
 }

 //How Pass model
/*
 printf("MODEL\n");
 for(int i=0;i<Npass;i++){
  //printf("%d,",Pass[i]);
  int pos=Pass[i];
  int res = p->ResNum[pos];
  int aa_type = p->TypeResId[pos];
  printf("ATOM  %5d  CA  %3s%6d    ",i+1,RES_NAMES[aa_type],res);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",cd[pos][0],cd[pos][1],cd[pos][2],
1.0,1.0);
 }
 puts("TER\nENDMDL");
*/


 //Fill small parts
 for(int i=0;i<Npass;i++){
  int pos1=Pass[i];
  int pos2=Pass[i+1];
  int gap = pos2 - pos1;
  Filled[Nfill++]=pos1;
  if( gap >1 && gap < 11 ){//Up to 10AAs gap
   	for(int j=pos1+1;j<pos2;j++){
	 Filled[Nfill++]=j;
	}
  }
 }
 Filled[Nfill++]=Pass[Npass-1];
/*
 printf("MODEL\n");
 for(int i=0;i<Nfill;i++){
  //printf("%d,",Pass[i]);
  int pos=Filled[i];
  int res = p->ResNum[pos];
  int aa_type = p->TypeResId[pos];
  printf("ATOM  %5d  CA  %3s%6d    ",i+1,RES_NAMES[aa_type],res);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",cd[pos][0],cd[pos][1],cd[pos][2],
1.0,1.0);
 }
 puts("TER\nENDMDL");
*/
 //Remove Isolated parts

 //Input Pass Model
 for(int j=0;j<Npass;j++){
  int pos=Pass[j];
  mod.Ca[j].real_cd[0]=cd[pos][0];
  mod.Ca[j].real_cd[1]=cd[pos][1];
  mod.Ca[j].real_cd[2]=cd[pos][2];
  mod.Ca[j].node_id = -1;//Not use
  mod.Ca[j].seq_id = p->ResNum[pos];//Start from 1
 }
 mod.Len=Npass;
 mod.Start=p->ResNum[Pass[0]];
 mod.End=Pass[Npass-1]+1;
 mod.Score=0;//Not yet
 mod.RawScore=0.0;
 mod.Rmsd=0.0;
 if(mod.Len > 10)
  Nmodel=AddModel(&mod,model,Nmodel,false);

 if(Npass==Nfill)
  return Nmodel;
 //Input Pass Model
 for(int j=0;j<Nfill;j++){
  int pos=Filled[j];
  mod.Ca[j].real_cd[0]=cd[pos][0];
  mod.Ca[j].real_cd[1]=cd[pos][1];
  mod.Ca[j].real_cd[2]=cd[pos][2];
  mod.Ca[j].node_id = -1;//Not use
  mod.Ca[j].seq_id = p->ResNum[pos];//Start from 1
 }
 mod.Len=Nfill;
 mod.Start=p->ResNum[Filled[0]];
 mod.End=Filled[Nfill-1]+1;
 mod.Score=0;//Not yet
 mod.RawScore=0.0;
 mod.Rmsd=0.0;
 if(mod.Len > 10)
  Nmodel=AddModel(&mod,model,Nmodel,false);

 return Nmodel;
}

//Add loops and NC-term
int DetectFittedRegions2(NODE *nodes[2000],double cd[2000][3],int Ncd,int ali[2000][2],int Nali,RESULT_MODEL **model,int Nmodel, int MinNode){

 RESULT_MODEL mod;
 //printf("#Nali= %d\n",Nali);
 int Filled[2000],Nfill=0;
 int Pass[2000],Npass=0;

 //Fill Small Gaps <=3AA
 for(int i=0;i<Nali;i++){
  int pos1=ali[i][1];
  int pos2=ali[i+1][1];
  int gap = pos2 - pos1;
  Filled[Nfill++]=pos1;
  if( gap >1 && gap < 2 ){//Up to 3AAs gap
        //printf("Gap:%d\n",gap);
        for(int j=pos1+1;j<pos2;j++){
         Filled[Nfill++]=j;
        }
  }
 }
 //puts("#Done check gap");
 //Filter
 for(int i=0;i<Nfill;i++){
  int pos1=Filled[i];
  int Lfrag=0;
	//To N-ter
	for(int j=i-1;j>=0;j--){
	 if((i-j)==(pos1-Filled[j]))
	  Lfrag++;
	 else
	  break;
	}
	//To C-ter
	for(int j=i+1;j<Nfill;j++){
	 if((j-i)==(Filled[j]-pos1))
	  Lfrag++;
	 else
	  break;
	}
  if(Lfrag>=4){//4+1=5AAs
   Pass[Npass++]=pos1;
  }
 }

 //puts("#Done filter gap");
 //Input Pass Model at least 5AAs
 for(int j=0;j<Npass;j++){
  int pos=Pass[j];
  mod.Ca[j].real_cd[0]=cd[pos][0];
  mod.Ca[j].real_cd[1]=cd[pos][1];
  mod.Ca[j].real_cd[2]=cd[pos][2];
  mod.Ca[j].node_id = -1;//Not use
  mod.Ca[j].seq_id = nodes[pos]->seq_id;//Start from 1
  mod.Ca[j].AAtype = nodes[pos]->AAtype;//Start from 1
  mod.Ca[j].cid = nodes[pos]->cid;//Start from 1
  //printf("%d cid= %d\n",mod.Ca[j].seq_id,mod.Ca[j].cid);
 }
 mod.Len=Npass;
 mod.Start=mod.Ca[0].seq_id;
 mod.End=mod.Ca[Npass-1].seq_id+1;
 mod.Score=0;//Not yet
 mod.RawScore=0.0;
 mod.Rmsd=0.0;
 mod.Weight=cmd.Weight;
 //printf("Start: %d\n",mod.Start);
 if(mod.Len < MinNode)
  return Nmodel;
 //printf("#Adding...L=%d Nmodel=%d\n",mod.Len,Nmodel);
 Nmodel=AddModel(&mod,model,Nmodel,false);
 //printf("#Done Adding...L=%d Nmodel=%d\n",mod.Len,Nmodel);

/*
 printf("MODEL\n");
 for(int i=0;i<Npass;i++){
  int pos=Pass[i];
  int res = mod.Ca[i].seq_id;
  int aa_type = mod.Ca[i].AAtype;
  printf("ATOM  %5d  CA  %3s%6d    ",i+1,RES_NAMES[aa_type],res);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",cd[pos][0],cd[pos][1],cd[pos][2],
1.0,1.0);
 }
 puts("TER\nENDMDL");
*/

 //Show tbl
 //for(int i=0;i<Npass;i++)
 // printf("##%d %d\n",i,Pass[i]);

 //Fill small parts
 for(int i=0;i<Npass;i++){
  int pos1=Pass[i];
  int pos2=Pass[i+1];
  int gap = pos2 - pos1;
  int StartPos=pos1-1;
  int EndPos=pos2+1;
  if(StartPos<0)
   StartPos=0;
  if(EndPos>=Ncd)
   EndPos=Ncd-1;
  
  
  if( gap >1 && gap < 30){//Up to 20AAs gap
        printf("#Adding..Gap:%d-%d %d\n",pos1,pos2,gap);
	int L=0;
        //for(int j=pos1-1;j<=pos2+1;j++){ //++GAP++
        for(int j=StartPos;j<=EndPos;j++){ //++GAP++
	 mod.Ca[L].real_cd[0]=cd[j][0];
	 mod.Ca[L].real_cd[1]=cd[j][1];
	 mod.Ca[L].real_cd[2]=cd[j][2];
	 mod.Ca[L].seq_id=nodes[j]->seq_id;
	 mod.Ca[L].AAtype=nodes[j]->AAtype;
	 mod.Ca[L].cid=nodes[j]->cid;
	 L++;
        }
	mod.Len=L;
 	mod.Start=mod.Ca[0].seq_id;
 	mod.End=mod.Ca[L-1].seq_id+1;
	mod.Score=0;//Not yet
 	mod.RawScore=0.0;
 	mod.Rmsd=0.0;
 	//mod.Weight=cmd.Weight;
 	mod.Weight=1.00;//Small parts
        //printf("#Adding..GapL:%d\n",L);
	Nmodel=AddModel(&mod,model,Nmodel,false);
        //printf("#DONE Adding..GapL:%d\n",L);
  }
 }
 //puts("#Done fill small-gap");
 //Fill N-ter
 if(Pass[0]<10){
  printf("#Fill N-ter 0-%d\n",Pass[0]);
  int L=0;
  int ExpLen=Pass[0];
  if(ExpLen<5)
   ExpLen=5;
        for(int j=0;j<=ExpLen;j++){ //++GAP++
         mod.Ca[L].real_cd[0]=cd[j][0];
         mod.Ca[L].real_cd[1]=cd[j][1];
         mod.Ca[L].real_cd[2]=cd[j][2];
         mod.Ca[L].seq_id=nodes[j]->seq_id;
         mod.Ca[L].AAtype=nodes[j]->AAtype;
         mod.Ca[L].cid=nodes[j]->cid;
         L++;
        }
        mod.Len=L;
        mod.Start=mod.Ca[0].seq_id;
        mod.End=mod.Ca[L-1].seq_id+1;
        mod.Score=0;//Not yet
        mod.RawScore=0.0;
        mod.Rmsd=0.0;
        //mod.Weight=cmd.Weight;
        mod.Weight=1.00;//Small parts
        Nmodel=AddModel(&mod,model,Nmodel,false);
 }
 //Fill C-ter
 if(Ncd-Pass[Npass-1]<10){
  printf("#Fill C-ter %d-%d\n",Pass[Npass-1],Ncd);
  int L=0;
  int ExpStart=Pass[Npass-1];
  if(Ncd-ExpStart<5)
   ExpStart=Ncd-5;
        for(int j=ExpStart;j<Ncd;j++){ //++GAP++
         mod.Ca[L].real_cd[0]=cd[j][0];
         mod.Ca[L].real_cd[1]=cd[j][1];
         mod.Ca[L].real_cd[2]=cd[j][2];
         mod.Ca[L].seq_id=nodes[j]->seq_id;
         mod.Ca[L].AAtype=nodes[j]->AAtype;
         mod.Ca[L].cid=nodes[j]->cid;
         L++;
        }
        mod.Len=L;
        mod.Start=mod.Ca[0].seq_id;
        mod.End=mod.Ca[L-1].seq_id+1;
        mod.Score=0;//Not yet
        mod.RawScore=0.0;
        mod.Rmsd=0.0;   
        //mod.Weight=cmd.Weight;
        mod.Weight=1.00;//Small parts
        Nmodel=AddModel(&mod,model,Nmodel,false);
 }
 return Nmodel;
}




int GenerateFragments(SEQ *seq,NODE_CSV *in,PDB *p,NODE_CSV *af2,int Fl,float Rcut,RESULT_MODEL **model){
 int **stbl;
 int Ntbl=0;
 int l;
 int Nmodel=0;
 if((stbl=(int **)malloc(sizeof(int *)*in->Npath*in->LongestPath*2))==NULL)
  return 0;

 //Find Unique path with Fl
 int Nstbl=GetNrPath(in,Fl,stbl);
 printf("##Nstbl=%d / %d\n",Nstbl,in->Npath*in->LongestPath*2);
/*
 for(int i=0;i<Nstbl;i++){
  double cd11[2000][3];
  GetFragFromTBL(in,stbl,i,Fl,cd11);
	  printf("MODEL\n");
	for(int j=0;j<Fl;j++){
 	 if(stbl[i][j]==-1)
	  break;
	 int pos=stbl[i][j];
  	printf("ATOM  %5d  CA  %3s%6d    ",i+1,RES_NAMES[0],j);
  	//printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",in->Ca[pos].real_cd[0],in->Ca[pos].real_cd[1],in->Ca[pos].real_cd[2],1.0,1.0);
  	printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",cd11[j][0],cd11[j][1],cd11[j][2],1.0,1.0);
 	}
 puts("TER\nENDMDL");
	
 }
*/
 //return 0;

 RESULT_MODEL **model_th[100];//Up to 100 threads
 int Nth=omp_get_max_threads();
 int Nmodel_th[100];
 for(int i=0;i<Nth;i++){
  Nmodel_th[i]=0;
  //if((model_th[i]=(RESULT_MODEL **)malloc(sizeof(RESULT_MODEL *)*af2->Nnode*in->Nnode))==NULL)
  if((model_th[i]=(RESULT_MODEL **)malloc(sizeof(RESULT_MODEL *)*af2->Nnode*Nstbl))==NULL)
   return 0;
 }

 #pragma omp parallel for schedule(dynamic,5)
 //for(int pdb_pos1=0;pdb_pos1<p->NumOfAtom;pdb_pos1++){
 for(int pdb_pos1=0;pdb_pos1<af2->Nnode-Fl;pdb_pos1++){
  //puts("Start");
  int th=omp_get_thread_num();

  //int ResNum = p->ResNum[pdb_pos1];//Start from 1
  int ResNum = af2->Ca[pdb_pos1].seq_id;//Start from 1
  //int order = p->AtomOnRes[pdb_pos1];

  double cd11[2000][3],cd21[2000][3],imposed[2000][3],pdb_cd[2000][3];
  int Nali,ali[2000][2];
  double rmsd;
  NODE *nodes[2000];

	//Use each path
        for(int pos1=0;pos1<Nstbl;pos1++){
         //if(GetFragFromPDB(p,pdb_pos1,Fl,cd21))
	 //puts("GetFrag");
         if(GetFragFromNODE(af2,pdb_pos1,Fl,cd21))
          continue;
	 //puts("GetTBL1");
         if(GetFragFromTBL(in,stbl,pos1,Fl,cd11))
          continue;
	 //puts("GetRMSD");
         fast_rmsd(cd11,cd21,Fl,&rmsd);
	 if(rmsd>Rcut)
	  continue;
	 //Initial impose
	 //if(GetFragFromPDB(p,pdb_pos1,Fl,cd21))//PDB coords
	 //puts("GetNODE");
         if(GetFragFromNODE(af2,pdb_pos1,Fl,cd21))
          continue;
	 //puts("GetTBL");
         if(GetFragFromTBL(in,stbl,pos1,Fl,cd11))
          continue;
	 Nali=Fl;
	 //printf("#RMSD= %f\n",rmsd);
	 int Npre=9;
	 int Nchain=0;
		for(int iter=0;iter < 10;iter++){
	 	 //if(GetFragFromPDB(p,0,p->NumOfAtom,pdb_cd))
	 	 Nchain=GetChainFragFromNODE(af2,pdb_pos1,pdb_cd,nodes);
		 //printf("Nchain= %d\n",Nchain);
	 	 //move_by_mtx(cd11,cd21,Nali,pdb_cd,imposed,p->NumOfAtom);
	 	 move_by_mtx(cd11,cd21,Nali,pdb_cd,imposed,Nchain);
	 	 //Find alignment
	 	 Nali=FindAlignment(ali,in,imposed,Nchain,Rcut);
	 	 //printf("#Lali%d= %d\n",iter,Nali);
	 	 //SetUp
	 	 //Nali=SetUpCd(ali,Nali,in,p,cd11,cd21);
	 	 Nali=SetUpCd2(ali,Nali,in,nodes,cd11,cd21);
		 if(Nali<=Npre){
		  //printf("##Done\n");
		  break;
		 }
		 Npre=Nali;
		}

	 //Add model
	 if(Nali>Fl){
	  //puts("#Adding");
	  Nmodel_th[th] = DetectFittedRegions2(nodes,imposed,Nchain,ali,Nali,model_th[th],Nmodel_th[th],cmd.MinReg);
	  //puts("#Done Adding");
	 }
	}
        printf("##Th %d pos:%d Nmodel= %d / %d\n",th,pdb_pos1,Nmodel_th[th],af2->Nnode*in->Nnode);
	//puts("END....");
 }

 //Merge model_th
 for(int i=0;i<Nth;i++){
	printf("#AddingThread %d Nmodel= %d\n",i,Nmodel_th[i]);
	for(int j=0;j<Nmodel_th[i];j++){
	 //printf("%d Start %d Nmodel= %d\n",j,model_th[i][j]->Start,Nmodel);
	 Nmodel=AddModel(model_th[i][j],model,Nmodel,false);
	}
	printf("##Thread%d Nmodel= %d Total= %d\n",i,Nmodel_th[i],Nmodel);
 }

 return Nmodel;
}
//Compare path & AF model
int GenerateFragmentsLocal(SEQ *seq,NODE_CSV *in,PDB *p,int Fl,float Rcut,RESULT_MODEL **model){

 //Malloc NR path tbl
 //Make NR path
 int tmp_path[100];
 int **tbl,**stbl;
 int Ntbl=0;
 int l;
 int Nmodel=0;
 if((tbl=(int **)malloc(sizeof(int *)*in->Npath*in->LongestPath*2))==NULL)
  return 0;
 if((stbl=(int **)malloc(sizeof(int *)*in->Npath*in->LongestPath*2))==NULL)
  return 0;

 //single
 int Nstbl=GetNrPath(in,Fl*2,stbl);
 printf("##Nstbl=%d / %d\n",Nstbl,in->Npath*in->LongestPath*2);
 for(int pdb_pos1=0;pdb_pos1<p->NumOfAtom;pdb_pos1++){
  int ResNum = p->ResNum[pdb_pos1];//Start from 1
  int order = p->AtomOnRes[pdb_pos1];
  double cd11[100][3],cd21[100][3],imposed[100][3];
  RESULT_MODEL tmp_model,homo_model;
  	for(int pos1=0;pos1<Nstbl;pos1++){
  	 if(GetFragFromPDB(p,pdb_pos1,Fl*2,cd21))
  	  continue;
	 double rmsd;
	 int Fl_single=2*Fl;
	 if(GetFragFromTBL(in,stbl,pos1,2*Fl,cd11))
	  continue;
	 rmsd=100;
		while(1){
	 	 fast_rmsd(cd11,cd21,Fl_single,&rmsd);
		 //printf("##%d %f Fl_single: %d\n",pos1,rmsd,Fl_single);
		 if(Fl_single<=Fl+1)
		  break;
		 if(rmsd<=Rcut)
		  break;
		 Fl_single--;
		}
		if(rmsd>Rcut)
		 continue;
		//printf("%d %f Fl_single: %d\n",pos1,rmsd,Fl_single);
		//Impose
		if(GetFragFromPDB(p,pdb_pos1,Fl_single,cd21))
  	 	 continue;
	 	if(GetFragFromTBL(in,stbl,pos1,Fl_single,cd11))
	 	 continue;
	 	move_by_mtx(cd11,cd21,Fl_single,cd21,imposed,Fl_single);

		int ResNum1 = p->ResNum[pdb_pos1];//Start from 1
		//float sco=0;
	  	for(int j=0;j<Fl_single;j++){
	 	 tmp_model.Ca[j].real_cd[0]=imposed[j][0];
	 	 tmp_model.Ca[j].real_cd[1]=imposed[j][1];
	 	 tmp_model.Ca[j].real_cd[2]=imposed[j][2];
	 	 tmp_model.Ca[j].node_id = stbl[pos1][j];
		 tmp_model.Ca[j].seq_id = ResNum1+j;//Start from 1
		 //printf("IN %d %d\n",tmp_model.Ca[j].seq_id ,tmp_model.Ca[j].node_id );
	 	}
	 	tmp_model.Len=Fl_single;
	 	tmp_model.Start=ResNum1;
	 	tmp_model.End=tmp_model.Start + Fl_single;
	 	tmp_model.Score=0;//Not yet
	 	tmp_model.RawScore=QScoring(&tmp_model,in);
	 	tmp_model.Rmsd=rmsd;
		//printf("sco= %f\n",tmp_model.RawScore/(float)Fl_single);
		if(tmp_model.RawScore>0){
		 Nmodel=AddModel(&tmp_model,model,Nmodel,true);
		 	//Check Homo-oligomer
			for(int spos=0;spos<seq->len;spos++){
			 if(spos==tmp_model.Start-1)
			  continue;
			  if(Homo_check(seq,spos,Fl,&tmp_model,&homo_model)==true){
			   //printf("HOMO %d\n",spos);
			   if(tmp_model.RawScore>0)
			    Nmodel=AddModel(&homo_model,model,Nmodel,true);
	 		 }
			}
		}
	}

	printf("##Nmodel= %d\n",Nmodel);
 }


 //two regions
 Ntbl=GetNrPath(in,Fl,tbl);

 int Npair=0;
 int *Pair;
 if((Pair=(int *)malloc(sizeof(int)*p->NumOfAtom*2*Ntbl))==NULL)
  return 0;
 //single position
 //take a fragment from PDB
 for(int pdb_pos1=0;pdb_pos1<p->NumOfAtom;pdb_pos1++){
  int ResNum = p->ResNum[pdb_pos1];//Start from 1
  int order = p->AtomOnRes[pdb_pos1];
  double cd11[100][3],cd21[100][3],imposed[100][3];
  double cd12[100][3],cd22[100][3];
  if(GetFragFromPDB(p,pdb_pos1,Fl,cd21))
   continue;
  //printf("POS %d\n",pdb_pos1);
	for(int pos1=0;pos1<Ntbl;pos1++){
	 double rmsd;
	 if(GetFragFromTBL(in,tbl,pos1,Fl,cd11))
	  continue;

	  //rms
	  fast_rmsd(cd11,cd21,Fl,&rmsd);
	  if(rmsd>Rcut)
	   continue;
	  if(!isfinite(rmsd))
	   continue;
	 //Keep similar fragment
	 Pair[Npair*2]=pdb_pos1;
	 Pair[Npair*2+1]=pos1;
	 Npair++;
	  //printf("##PDB %d  vs %d RMSD= %f %d\n",pdb_pos1,pos1,rmsd,Npair);
	}
 }
 printf("##Npair= %d\n",Npair);
 //Check Closest Dist between PDB pairs
 int *ContPairs;
 int Ncont=0;
 if((ContPairs=(int *)calloc(sizeof(int),p->NumOfAtom*p->NumOfAtom))==NULL)
  return 0;
 for(int pdb_pos1=0;pdb_pos1<p->NumOfAtom;pdb_pos1++){
  int ResNum = p->ResNum[pdb_pos1];//Start from 1
  int order = p->AtomOnRes[pdb_pos1];
  double cd1[100][3],cd2[100][3];
  if(GetFragFromPDB(p,pdb_pos1,Fl,cd1))
   continue;
  for(int pdb_pos2=pdb_pos1+Fl;pdb_pos2<p->NumOfAtom;pdb_pos2++){
   if(GetFragFromPDB(p,pdb_pos2,Fl,cd2))
    continue;
	float MinD=100.0;
	float d2;
	for(int i=0;i<Fl;i++){
	 for(int j=0;j<Fl;j++){
	  d2=(cd1[i][0]-cd2[j][0])*(cd1[i][0]-cd2[j][0])
	    +(cd1[i][1]-cd2[j][1])*(cd1[i][1]-cd2[j][1])
	    +(cd1[i][2]-cd2[j][2])*(cd1[i][2]-cd2[j][2]);
	  if(d2<MinD)
	   MinD=d2;
	 }
	}

   if(MinD>=64.0)//8A*8A
    continue;
   //printf("##PDBfrag %d %d MinD= %f\n",pdb_pos1,pdb_pos2,sqrt(MinD));
   ContPairs[p->NumOfAtom*pdb_pos1+pdb_pos2]=1;
   ContPairs[p->NumOfAtom*pdb_pos2+pdb_pos1]=1;
  }
 }
 //Pair

 RESULT_MODEL **model_th[100];//Up to 100 threads
 int Nth=omp_get_max_threads();
 int Nmodel_th[100];
 for(int i=0;i<Nth;i++){
  Nmodel_th[i]=0;
  if((model_th[i]=(RESULT_MODEL **)malloc(sizeof(RESULT_MODEL *)*p->NumOfAtom*in->Nnode))==NULL)
   return 0;
 }

 #pragma omp parallel for schedule(dynamic,5)
 for(int i=0;i<Npair;i++){
  RESULT_MODEL tmp_model,homo_model;
  int th=omp_get_thread_num();
  int pdb_pos1=Pair[2*i];
  int tbl_pos1=Pair[2*i+1];
  double cd11[100][3],cd21[100][3],imposed[100][3];
  double cd12[100][3],cd22[100][3];
  double cd1[100][3],cd2[100][3];
  int SeqGap = Fl;
  	for(int j=i+1;j<Npair;j++){
  	 int pdb_pos2=Pair[2*j];
  	 int tbl_pos2=Pair[2*j+1];
	 double rmsd;
	 if(pdb_pos1+Fl>pdb_pos2)//Overlap
	  continue;
	 //Different chain ID
	 if(p->Chain[pdb_pos1]!=p->Chain[pdb_pos2])
	  continue;

	 //Remove No-contact
	 if(ContPairs[pdb_pos1*p->NumOfRes+pdb_pos2]==0)
	  continue;

	 //Overlap Check
	 if(OvCheck(tbl,tbl_pos1,tbl_pos2,Fl)==true)
	  continue;

	 if(tbl_pos1==tbl_pos2)
	  continue;
	 if(GetFragFromTBL(in,tbl,tbl_pos1,Fl,cd11))
	  continue;
	 if(GetFragFromTBL(in,tbl,tbl_pos2,Fl,cd12))
	  continue;
  	 if(GetFragFromPDB(p,pdb_pos1,Fl,cd21))
  	  continue;
  	 if(GetFragFromPDB(p,pdb_pos2,Fl,cd22))
  	  continue;

	 //LB3D filter
	 if(lb3d(cd11,cd12,cd21,cd22,Fl)>Rcut)
	  continue;


	 //continue;
	 //Merge
	 for(int k=0;k<Fl;k++){
	  cd11[Fl+k][0]=cd12[k][0]; cd11[Fl+k][1]=cd12[k][1]; cd11[Fl+k][2]=cd12[k][2];
	  cd21[Fl+k][0]=cd22[k][0]; cd21[Fl+k][1]=cd22[k][1]; cd21[Fl+k][2]=cd22[k][2];
	 }
	 fast_rmsd(cd11,cd21,Fl*2,&rmsd);
	 if(rmsd>Rcut)
	  continue;
	 if(!isfinite(rmsd))
	  continue;
	
	//Impose
	 if(GetFragFromTBL(in,tbl,tbl_pos1,Fl,cd11))
	  continue;
	 if(GetFragFromTBL(in,tbl,tbl_pos2,Fl,cd12))
	  continue;
  	 if(GetFragFromPDB(p,pdb_pos1,Fl,cd21))
  	  continue;
  	 if(GetFragFromPDB(p,pdb_pos2,Fl,cd22))
  	  continue;
	 //Merge
	 for(int k=0;k<Fl;k++){
	  cd11[Fl+k][0]=cd12[k][0]; cd11[Fl+k][1]=cd12[k][1]; cd11[Fl+k][2]=cd12[k][2];
	  cd21[Fl+k][0]=cd22[k][0]; cd21[Fl+k][1]=cd22[k][1]; cd21[Fl+k][2]=cd22[k][2];
	 }
	 move_by_mtx(cd11,cd21,Fl*2,cd21,imposed,Fl*2);

	  //input data 1st+2nd region
	 int ResNum1 = p->ResNum[pdb_pos1];//Start from 1
	 int ResNum2 = p->ResNum[pdb_pos2];//Start from 1
	 for(int j=0;j<Fl*2;j++){
	  tmp_model.Ca[j].real_cd[0]=imposed[j][0];
	  tmp_model.Ca[j].real_cd[1]=imposed[j][1];
	  tmp_model.Ca[j].real_cd[2]=imposed[j][2];
	  //tmp_model.Ca[j].node_id = -1;
	  //new
	  if(j<Fl){
	   tmp_model.Ca[j].seq_id = ResNum1+j;//Start from 1
	   tmp_model.Ca[j].node_id = tbl[tbl_pos1][j];//Start from 1
	  }
 	  else{
	   tmp_model.Ca[j].seq_id = ResNum2+j-Fl;//Start from 1
	   tmp_model.Ca[j].node_id = tbl[tbl_pos2][j-Fl];//Start from 1
	  }
	 }
	 tmp_model.Len=Fl*2;//Number of Coords;
	 tmp_model.Start=ResNum1;
	 tmp_model.End=ResNum2 + Fl;
	 tmp_model.Score=QScoring(&tmp_model,in);//Not yet
	 tmp_model.RawScore=0;//not yet
	 tmp_model.Rmsd=rmsd;
	 if(tmp_model.Score>0){
	  Nmodel_th[th]=AddModel(&tmp_model,model_th[th],Nmodel_th[th],true);
	 
	  int PosGap=ResNum2-ResNum1;
	  bool homo_check;


	 	for(int spos=0;spos<seq->len;spos++){
	 	 if(spos==tmp_model.Start-1)
	 	  continue;
	 	 if(Homo_check(seq,spos,Fl,&tmp_model,&homo_model)==true){
	 	  //printf("OK\n");
	 	  Nmodel_th[th]=AddModel(&homo_model,model_th[th],Nmodel_th[th],true);
	 	 }
	 	}
	 }


	 //input extended region
	 //pos1 -- Gap -- pos2
	 if (pdb_pos1 + Fl + SeqGap > pdb_pos2 && tmp_model.Score>0){
	  //printf("pdbpos1 %d-%d pos2 %d-%d\n",pdb_pos1,pdb_pos1+Fl, pdb_pos2,pdb_pos2+Fl);
	  int Fl2=(pdb_pos2+Fl)-pdb_pos1;
	  double ext[100][3];
	  //Impose
	  if(GetFragFromTBL(in,tbl,tbl_pos1,Fl,cd11))
	  continue;
	  if(GetFragFromTBL(in,tbl,tbl_pos2,Fl,cd12))
	   continue;
  	  if(GetFragFromPDB(p,pdb_pos1,Fl,cd21))
  	   continue;
  	  if(GetFragFromPDB(p,pdb_pos2,Fl,cd22))
  	   continue;
	  //Merge
	  for(int k=0;k<Fl;k++){
	   cd11[Fl+k][0]=cd12[k][0]; cd11[Fl+k][1]=cd12[k][1]; cd11[Fl+k][2]=cd12[k][2];
	   cd21[Fl+k][0]=cd22[k][0]; cd21[Fl+k][1]=cd22[k][1]; cd21[Fl+k][2]=cd22[k][2];
	  }
	  move_by_mtx(cd11,cd21,2*Fl,ext,imposed,Fl2);
	  //printf("##Adding %d, %d\n",pdb_pos1,pdb_pos2);
	  	//input data 1st-gap-2nd region
	 	int ResNum2 = p->ResNum[pdb_pos1];//Start from 1
	  	for(int j=0;j<Fl2;j++){
	 	 tmp_model.Ca[j].real_cd[0]=imposed[j][0];
	 	 tmp_model.Ca[j].real_cd[1]=imposed[j][1];
	 	 tmp_model.Ca[j].real_cd[2]=imposed[j][2];
	 	 tmp_model.Ca[j].node_id = -1;
		 tmp_model.Ca[j].seq_id = ResNum2+j;
	 	}
	 	tmp_model.Start=ResNum2;
	 	tmp_model.End=tmp_model.Start + Fl2;
	 	tmp_model.Score=0;//Not yet
	 	tmp_model.RawScore=0;//not yet
	 	tmp_model.Rmsd=rmsd;
	 	tmp_model.Len=Fl2;
		Nmodel_th[th]=AddModel(&tmp_model,model_th[th],Nmodel_th[th],true);

		//Check Homo-oligomer
		for(int spos=0;spos<seq->len;spos++){
	 	 if(spos==tmp_model.Start-1)
	 	  continue;
	  	 if(Homo_check(seq,spos,Fl,&tmp_model,&homo_model)==true){
	   	  Nmodel_th[th]=AddModel(&homo_model,model_th[th],Nmodel_th[th],true);
	  	 }
		}
	 }
	}
	if(i%100==0)
   	 printf("##%d/%d thd= %d Nmodel= %d\n",i,Npair,th,Nmodel_th[th]);
 }

 //Merge model_th
 //Nmodel=0;
 for(int i=0;i<Nth;i++){
	for(int j=0;j<Nmodel_th[i];j++){
	 Nmodel=AddModel(model_th[i][j],model,Nmodel,true);
	}
	printf("##Thread%d Nmodel= %d Total= %d\n",i,Nmodel_th[i],Nmodel);
 }

 return Nmodel;
}

int AddModel_ori(RESULT_MODEL *in,RESULT_MODEL **all, int N){
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
		//Different Sequence
		if(in->Ca[p].seq_id!=all[i]->Ca[p].seq_id){
		 DisSim=true;
		 break;
		}
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
	//if(FindSimilar==true && in->Score > all[i]->Score){
	if(FindSimilar==true && in->Rmsd < all[i]->Rmsd){
	 //printf("Swap.. %d\n",i);
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
	 	 all[i]->Ca[p].seq_id=in->Ca[p].seq_id;
	 	 all[i]->Ca[p].cid=in->Ca[p].cid;
		}
	}
  	if(FindSimilar==true)
   	 break;
 }
 int i=Num;
 if(FindSimilar==false){//Add New Data
  //puts("Add New Model");
  if((all[i]=(RESULT_MODEL *)malloc(sizeof(RESULT_MODEL)))==NULL)
   return 0;
  //puts("DONE malloc");

  Num++;

  all[i]->Start = in->Start;
  all[i]->End = in->End;
  all[i]->Score = in->Score;
  all[i]->RawScore = in->RawScore;
  all[i]->Rmsd = in->Rmsd;
  all[i]->Len = in->Len;
  //puts("ADD Data");
  //printf("%d\n",in->Len);

	for(int p=0;p<in->Len;p++){
	  all[i]->Ca[p].real_cd[0]=in->Ca[p].real_cd[0];
	  all[i]->Ca[p].real_cd[1]=in->Ca[p].real_cd[1];
	  all[i]->Ca[p].real_cd[2]=in->Ca[p].real_cd[2];
	  all[i]->Ca[p].AAtype=in->Ca[p].AAtype;
	  all[i]->Ca[p].node_id=in->Ca[p].node_id;
	  all[i]->Ca[p].seq_id=in->Ca[p].seq_id;
	  all[i]->Ca[p].cid=in->Ca[p].cid;
	}
 }
 //printf("END Comp... %d\n",Num);
 return Num;
}

float dist2(float a[3],float b[3]){
 float d2=(a[0]-b[0])*(a[0]-b[0])+
	 (a[1]-b[1])*(a[1]-b[1])+
	 (a[2]-b[2])*(a[2]-b[2]);
 return d2;
}


void CopyResultModel(RESULT_MODEL *in,RESULT_MODEL *out){
  out->Start = in->Start;
  out->End = in->End;
  out->Score = in->Score;
  out->RawScore = in->RawScore;
  out->Rmsd = in->Rmsd;
  out->Len = in->Len;
  out->Weight = in->Weight;
        for(int p=0;p<in->Len;p++){
          out->Ca[p].real_cd[0]=in->Ca[p].real_cd[0];
          out->Ca[p].real_cd[1]=in->Ca[p].real_cd[1];
          out->Ca[p].real_cd[2]=in->Ca[p].real_cd[2];
          out->Ca[p].AAtype=in->Ca[p].AAtype;
          out->Ca[p].node_id=in->Ca[p].node_id;
          out->Ca[p].seq_id=in->Ca[p].seq_id;
          out->Ca[p].cid=in->Ca[p].cid;
        }
}



//*****NEED TO FIX********************
int AddModel(RESULT_MODEL *in,RESULT_MODEL **all, int N,bool mode){
 //puts("Start COmpare....");
 int Num=N;
 bool DisSim=false;
 bool FindSimilar=false;
 bool AddFlag=true;
 bool NewFlag=true;
 bool SkipFlag=false;
 float d2=0;
 float d2cut=1.00;
 int RemoveList[10000];
 int Nrm=0;
 for(int i=0;i<Num;i++){
  FindSimilar=false;
  //printf("Comp..%d\n",i);
  if(all[i]->Len==0){
   RemoveList[Nrm++]=i;
   continue;
  }
//---in---
//  -all----
  if(in->Start < all[i]->Start && in->End < all[i]->End )
   continue;
//  ---in---
//---all--
  if(in->Start > all[i]->Start && in->End > all[i]->End )
   continue;

  	//Check positions
	DisSim=false;
	if(in->Len < all[i]->Len){
	//  ++in++
	//++++all+++
		for(int p=0;p<in->Len;p++){//shorter
		 bool find=false;
			for(int q=0;q<all[i]->Len;q++){
			 if(in->Ca[p].seq_id==all[i]->Ca[q].seq_id){
			  if(dist2(in->Ca[p].real_cd,all[i]->Ca[q].real_cd)<1.00)
			   find=true;
			  break;
			 }
			}
		 if(find==false){
		  DisSim=true;
		   break;
		 }
		}
	}else{
		for(int q=0;q<all[i]->Len;q++){//shorter
		 bool find=false;
			for(int p=0;p<in->Len;p++){
			 if(in->Ca[p].seq_id==all[i]->Ca[q].seq_id){
			  if(dist2(in->Ca[p].real_cd,all[i]->Ca[q].real_cd)<1.00)
			   find=true;
			  break;
			 }
			}
		 if(find==false){
		  DisSim=true;
		   break;
		 }
		}

	}
	if (DisSim==false){
	 FindSimilar=true;
	 NewFlag=false;
	}
	//printf("Step1\n");
	//Compare and Swap
	if(FindSimilar==true){
	 if(mode==true){
	  if(in->Len > all[i]->Len){//Keep longer
	   //all[i]->Len = 0;
	   RemoveList[Nrm++]=i;
	  }else{
	   SkipFlag=true;//Do not add
	   break;
	  }
	 }
	 //always filter
	 if(in->Len == all[i]->Len && in->RawScore >=all[i]->RawScore ){
	  //printf("##Swap.. %d\n",i);
	  //all[i]->Len = 0;
	  RemoveList[Nrm++]=i;
	 }else{
	  SkipFlag=true;
	  break;
	 }
	}
 }
 //Do not update
 if(SkipFlag==true){
  return Num;
 }else{
  //printf("Nrm= %d\n",Nrm);
  	if(Nrm==0){//No remove
	 printf("###Assign Memory..%d\n",Num);
	 if((all[Num]=(RESULT_MODEL *)malloc(sizeof(RESULT_MODEL)))==NULL)
          return 0;
 	 CopyResultModel(in,all[Num]);
	 Num++;
	}else{
	 int FirstRm=RemoveList[0];
	 //printf("###Swappin..%d\n",FirstRm);
 	 CopyResultModel(in,all[FirstRm]);
	 //Update Others
		for(int i=1;i<Nrm;i++){
		 all[RemoveList[i]]->Len=0;
		}
	 //printf("###Done..%d\n",FirstRm);
	}

 }

/*
 int Num2=0;
 for(int i=0;i<Num;i++){
  if(all[i]->Len!=0){
   all[Num2]=all[i];
   all[i]=NULL;
   Num2++;
  }else{
   printf("Skip...%d\n",i);
   all[i]=NULL;
  }
 }
*/
/*
 for(int i=0;i<Num;i++){
  if(all[i]->Len == 0){//Need to remove
   printf("Swap.. %d\n",i);
   for(int j=i;j+1<Num;j++){
    all[j]=all[j+1];
    all[j+1]=NULL;
   }
   Num--;
   i--;
   AddFlag=true;
  }
 }
*/
/*
 printf("Step2 %d\n",Num2);
 //Add
 //int i=Num;
 int i=Num2;
 if(AddFlag==true||NewFlag==true){//Add New Data
  printf("#Add New Model to %d\n",Num2);
  if(all[i]==NULL||i==0||NewFlag==true){
   puts("Malloc");
   if((all[i]=(RESULT_MODEL *)malloc(sizeof(RESULT_MODEL)))==NULL)
    return 0;
  }
  Num2++;

  all[i]->Start = in->Start;
  all[i]->End = in->End;
  all[i]->Score = in->Score;
  all[i]->RawScore = in->RawScore;
  all[i]->Rmsd = in->Rmsd;
  all[i]->Len = in->Len;
  all[i]->Weight = in->Weight;
  puts("ADD Data");
  printf("%d\n",in->Len);

	for(int p=0;p<in->Len;p++){
	 
	  all[i]->Ca[p].real_cd[0]=in->Ca[p].real_cd[0];
	  all[i]->Ca[p].real_cd[1]=in->Ca[p].real_cd[1];
	  all[i]->Ca[p].real_cd[2]=in->Ca[p].real_cd[2];
	  all[i]->Ca[p].AAtype=in->Ca[p].AAtype;
	  all[i]->Ca[p].node_id=in->Ca[p].node_id;
	  all[i]->Ca[p].seq_id=in->Ca[p].seq_id;
	  all[i]->Ca[p].cid=in->Ca[p].cid;
	}
 }
 */
 return Num;
}

typedef struct{
  int Dia;
  int Sco;
 } MTX;

//mode1 reverse
//mode2 non-negative
int SimpleDP(int **OUT,int *p,int Lp,int Lseq,int **Smtx,MTX **dpmtx,int **smtx, int Niter,bool mode1,bool mode2){
 int pth[10000];
 int GapNodePen = -10000;
 int GapAAPen = -10000;
 //same order mode1
 if(mode1==false){
  for(int i=0;i<Lp;i++)
   pth[i]=p[i];
 }else{//reverse
  for(int i=0;i<Lp;i++)
   pth[i]=p[Lp-i-1];
 }
 //init smtx
 for(int i=0;i<Lseq+1;i++)
  for(int j=0;j<Lp+1;j++)
   smtx[i][j]=0;
 //Fill smtx
 for(int i=1;i<Lseq+1;i++){//Sequence starts from 1
  for(int j=1;j<Lp+1;j++){//Path starts from 1
   smtx[i][j]=Smtx[i][pth[j-1]];//smtx from 1,1 seq[1:] and path[1:]
   if(mode2==true && smtx[i][j] < 0) //ignore negative
    smtx[i][j]=0;
   //printf("seq%d Node%d %d %d\n",i,j,pth[j-1],smtx[i][j]);
 }}

 for(int iter=0;iter<Niter;iter++){
  	int Hsco;
  	//init
  	for(int i=0;i<Lseq+1;i++){
  	for(int j=0;j<Lp+1;j++){
  	 dpmtx[i][j].Dia=0;
  	 dpmtx[i][j].Sco=0;
  	}}

	//Fill dpmtx
	for(int x=1;x<Lseq+1;x++){
		for(int y=1;y<Lp+1;y++){
		 int UpSco   = dpmtx[x][y-1].Sco + GapNodePen;
		 int LeftSco = dpmtx[x-1][y].Sco + GapAAPen;
		 int DiaSco = dpmtx[x-1][y-1].Sco + smtx[x][y];


		 if(DiaSco <=0 && LeftSco <=0 && UpSco <=0 ){
		  dpmtx[x][y].Dia = 0;
		 }else if(DiaSco >= UpSco && DiaSco >= LeftSco){
		  dpmtx[x][y].Sco = DiaSco;
		  dpmtx[x][y].Dia = 1;
		 }else if(UpSco >= DiaSco && UpSco >= LeftSco){
		  dpmtx[x][y].Sco = UpSco;
		  dpmtx[x][y].Dia = 2;
		 }else if(LeftSco >= DiaSco && LeftSco >= UpSco){
		  dpmtx[x][y].Sco = LeftSco;
		  dpmtx[x][y].Dia = 3;
		 }

		}
	}
	//Find highest
	Hsco = 0;
	int Hpos[2];
	Hpos[0]=Hpos[1]=0;
	for(int x=1;x<Lseq+1;x++){
		for(int y=1;y<Lp+1;y++){
		 if(dpmtx[x][y].Dia != 1)
		  continue;
		 if(dpmtx[x][y].Sco >= Hsco){
		  Hsco = dpmtx[x][y].Sco;
		  Hpos[0]=x;
		  Hpos[1]=y;
		 }
		}
	}
	if(Hsco == 0){//No alignments
		for(int iter2=iter;iter2<Niter;iter2++){
	 		OUT[iter2][0]=-1;
			OUT[iter2][1]=-1;
		}
	 	break;
	}
	//printf("#TraceBack:%d Hsco= %d Hpos= %d %d\n",iter,Hsco,Hpos[0],Hpos[1]);
	int ALI[10000][2];
	int Lali=0;

	//Trace Back
	while(1){
		if(dpmtx[Hpos[0]][Hpos[1]].Dia == 1){//Dia
		 ALI[Lali][0]=Hpos[0];
		 ALI[Lali][1]=Hpos[1];
		 Hpos[0]--;
		 Hpos[1]--;
		}
		else if(dpmtx[Hpos[0]][Hpos[1]].Dia == 2){//Up
		 ALI[Lali][0]=-1;
		 ALI[Lali][1]=Hpos[1];
		 Hpos[1]--;
		}
		else if(dpmtx[Hpos[0]][Hpos[1]].Dia == 3){//Left
		 ALI[Lali][0]=Hpos[0];
		 ALI[Lali][1]=-1;
		 Hpos[0]--;
		}
		else if(dpmtx[Hpos[0]][Hpos[1]].Dia == 0){//Stop
		 break;
		}
		Lali++;
	}
	//printf("Lali=%d\n",Lali);
	//Convert
	int Lout=0;

	for(int i=0;i<Lali;i++){
	 int seq_pos=ALI[Lali-1-i][0];//start from 1
	 int path_pos=ALI[Lali-1-i][1];//start from 1
	 int node_id = pth[ALI[Lali-1-i][1]-1];//start from 0
	 //printf("%d %d seq:%d path:%d sco= %d\n",iter,Lout,seq_pos,node_id,Smtx[seq_pos][node_id]);
	 OUT[iter][2*Lout]=seq_pos;
	 OUT[iter][2*Lout+1]=node_id;
	 Lout++;

		//masking
		if(seq_pos > 0 && path_pos > 0 && smtx[seq_pos][path_pos] > 0)
		 smtx[seq_pos][path_pos] = 0;



	}
	OUT[iter][2*Lout]=-1;
	OUT[iter][2*Lout+1]=-1;
 }

 return 0;
}

int ConvertToModelFrg(SEQ *seq,int **OUT,int id,NODE_CSV *in, RESULT_MODEL *mod,RESULT_MODEL **model,int Nmodel,bool mode){
 //printf("#Converting...%d\n",id);
 int Flen=9;
 bool stop=false;

	for(int pos=0;pos<10000;pos+=3){
	 int L=0;
	 float sco=0.0;
	 for(int i=0;i<Flen;i++){
	  int seq_id=OUT[id][2*(i+pos)];
	  int node_id=OUT[id][2*(i+pos)+1];
	  //printf("%d %d %d-%d\n",id,i,seq_id,node_id);
	  if(seq_id == -1 ||node_id < 0 ){
	   stop=true;
	   break;
	  }

	  mod->Ca[i].real_cd[0]=in->Ca[node_id].real_cd[0];
	  mod->Ca[i].real_cd[1]=in->Ca[node_id].real_cd[1];
	  mod->Ca[i].real_cd[2]=in->Ca[node_id].real_cd[2];
	  mod->Ca[i].node_id = -1;
	  mod->Ca[i].seq_id = seq_id;
	  mod->Ca[i].AAtype = seq->seq_code[seq_id-1];
	  mod->Ca[i].cid = seq->chain_id[seq_id-1];
	  sco += in->Smtx[seq_id][node_id];
	
	  L++;
	 }
	 if(stop==true)
	  break;
 	 mod->Start=mod->Ca[0].seq_id;
 	 mod->End=mod->Ca[L-1].seq_id+1;
 	 mod->Score=0;//Not yet
 	 mod->RawScore=sco;//temp
 	 mod->Rmsd=0.0;
 	 mod->Len=L;

	 //printf("sco= %f L= %d\n",sco,L);
 	 Nmodel=AddModel(mod,model,Nmodel,false);
	}
 return Nmodel;

}


//mode true remove shorter models
//mode false no-filtering
int ConvertToModel(SEQ *seq,int **OUT,int id,NODE_CSV *in, RESULT_MODEL *mod,RESULT_MODEL **model,int Nmodel,int mask_tbl[],bool mode){
 printf("#Converting...%d\n",id);
 printf("#Npath: %d\n",in->Npath);
 RESULT_MODEL tmp_mod;
 float dis_d2=7.0*7.0;
	 int L=0;
	 float sco=0.0;
	 //for(int i=0;i<10000;i++){
	 for(int i=0;i<5000;i++){//Up to 5k-long alighnment
	  int seq_id=OUT[id][2*(i)];
	  int node_id=OUT[id][2*(i)+1];
	  printf("ALI %d %d %d-%d\n",id,i,seq_id,node_id);
	  if(seq_id == -1||node_id<0||node_id>=in->Nnode){ //ignore strange alignments
	   break;
	  }
	  //Masking
	  if(mask_tbl[seq_id]>5)
	   continue;

	  mod->Ca[L].real_cd[0]=in->Ca[node_id].real_cd[0];
	  mod->Ca[L].real_cd[1]=in->Ca[node_id].real_cd[1];
	  mod->Ca[L].real_cd[2]=in->Ca[node_id].real_cd[2];
	  mod->Ca[L].node_id = -1;
	  mod->Ca[L].seq_id = seq_id;
	  mod->Ca[L].AAtype = seq->seq_code[seq_id-1];
	  mod->Ca[L].cid = seq->chain_id[seq_id-1];
	  sco += in->Smtx[seq_id][node_id];
	
	  L++;
	 }
	 if(L<5)//Remove small frag
	  return Nmodel;
 	 mod->Start=mod->Ca[0].seq_id;
 	 mod->End=mod->Ca[L-1].seq_id+1;
 	 mod->Score=0;//Not yet
 	 mod->RawScore=sco;//temp
 	 mod->Rmsd=0.0;
 	 mod->Len=L;
 	 mod->Weight=1.00;//Threading w=1.0

	 //printf("sco= %f L= %d\n",sco,L);
 	 //Nmodel=AddModel(mod,model,Nmodel,false);

	//Find strange bonds
	int PartL=0;
	sco=0;
	for(int i=0;i<L;i++){
	 tmp_mod.Ca[PartL].real_cd[0]=mod->Ca[i].real_cd[0];
	 tmp_mod.Ca[PartL].real_cd[1]=mod->Ca[i].real_cd[1];
	 tmp_mod.Ca[PartL].real_cd[2]=mod->Ca[i].real_cd[2];
         tmp_mod.Ca[PartL].node_id=mod->Ca[i].node_id;
         tmp_mod.Ca[PartL].seq_id=mod->Ca[i].seq_id;
         tmp_mod.Ca[PartL].cid=mod->Ca[i].cid;
	 tmp_mod.Ca[PartL].AAtype = seq->seq_code[mod->Ca[i].seq_id-1];
         sco += in->Smtx[tmp_mod.Ca[PartL].seq_id][tmp_mod.Ca[PartL].node_id];
	 PartL++;
	 float d2=0.0;
	 if(i+1<L)
	  d2=dist2(mod->Ca[i].real_cd,mod->Ca[i+1].real_cd);
	 if(d2>dis_d2||i==L-1){//disconnect or end
	  printf("#Disconnections %d-%d %d-%d/%d %f\n",i,i+1,mod->Ca[i].seq_id,mod->Ca[i+1].seq_id,L,sqrt(d2));
	  //input
	  tmp_mod.Start=tmp_mod.Ca[0].seq_id;
          tmp_mod.End=tmp_mod.Ca[PartL-1].seq_id+1;
          tmp_mod.Score=0;//Not yet
          tmp_mod.RawScore=sco;//temp
          tmp_mod.Rmsd=0.0;
          tmp_mod.Len=PartL;
	  tmp_mod.Weight=1.0;
	  //reset
	  PartL=0;
	  sco=0;
	  //printf("Adding.. %d-%d\n",tmp_mod.Start,tmp_mod.End);
	  if(tmp_mod.Len > 4){
	   Nmodel=AddModel(&tmp_mod,model,Nmodel,false);
	  }
	 }


	}
 printf("###Nmodel= %d\n",Nmodel);
 return Nmodel;

}




int Threading(SEQ *seq,NODE_CSV *in,RESULT_MODEL **model,int Nmodel,int Niter, bool mask){

 int **smtx;
 int **OUT;//Alignments Up to 50
 RESULT_MODEL mod;
 if((smtx=(int **)malloc(sizeof(int *)*(seq->len+1)))==NULL)
  return 0;
 for(int i=0;i<seq->len+1;i++)
  if((smtx[i]=(int *)malloc(sizeof(int)*(in->Nnode+1)))==NULL)
   return 0;
 
 if((OUT=(int **)malloc(sizeof(int*)*50))==NULL)
  return 0;
 for(int i=0;i<50;i++)
  if((OUT[i]=(int *)malloc(sizeof(int)*10000))==NULL)
   return 0;

 MTX **dpmtx;
 if((dpmtx=(MTX**)malloc(sizeof(MTX *)*(seq->len+1)))==NULL)
  return 0;
 for(int i=0;i<seq->len+1;i++)
  if((dpmtx[i]=(MTX *)malloc(sizeof(MTX)*(in->Nnode+1)))==NULL)
   return 0;

 int mask_tbl[10000];
 for(int i=0;i<seq->len+1;i++)
  mask_tbl[i]=0;

 	if(mask==true){ 

 		for(int m=0;m<Nmodel;m++){
	 	//if(model[m]->Weight==1.0)
	 	// continue;
	 	for(int pos=0;pos<model[m]->Len;pos++){
  		  mask_tbl[model[m]->Ca[pos].seq_id]++;
	 	}
 		}
		for(int i=0;i<seq->len+1;i++){
	 	printf("seq %d cnt=%d\n",i,mask_tbl[i]);
		}
 	}



 for(int p=0;p<in->Npath;p++){
  	printf("PATH%d L= %d/%d\n",p,in->Lpath[p],in->Npath);
	if(in->Lpath[p] < 5) //ignore too small path
	 continue;
	//Original
	SimpleDP(OUT,in->Path[p],in->Lpath[p],seq->len,in->Smtx,dpmtx,smtx,Niter,false,false);
	for(int iter=0;iter<Niter;iter++){
	 Nmodel=ConvertToModel(seq,OUT,iter,in,&mod,model,Nmodel,mask_tbl,false);//No length filtering
	 //Nmodel=ConvertToModelFrg(seq,OUT,iter,in,&mod,model,Nmodel,false);
	}

	//Reverse
	SimpleDP(OUT,in->Path[p],in->Lpath[p],seq->len,in->Smtx,dpmtx,smtx,Niter,true,false);
	for(int iter=0;iter<Niter;iter++){
	 Nmodel=ConvertToModel(seq,OUT,iter,in,&mod,model,Nmodel,mask_tbl,false);
	 //Nmodel=ConvertToModelFrg(seq,OUT,iter,in,&mod,model,Nmodel,false);
	}

	//Non negative
	SimpleDP(OUT,in->Path[p],in->Lpath[p],seq->len,in->Smtx,dpmtx,smtx,Niter,false,true);
	for(int iter=0;iter<Niter;iter++){
	 //Nmodel=ConvertToModel(seq,OUT,iter,in,&mod,model,Nmodel,false);
	 Nmodel=ConvertToModel(seq,OUT,iter,in,&mod,model,Nmodel,mask_tbl,false);
	 //Nmodel=ConvertToModelFrg(seq,OUT,iter,in,&mod,model,Nmodel,false);
	}

	//Reverse Non-negative
	SimpleDP(OUT,in->Path[p],in->Lpath[p],seq->len,in->Smtx,dpmtx,smtx,Niter,true,true);
	for(int iter=0;iter<Niter;iter++){
	 //Nmodel=ConvertToModel(seq,OUT,iter,in,&mod,model,Nmodel,false);
	 Nmodel=ConvertToModel(seq,OUT,iter,in,&mod,model,Nmodel,mask_tbl,false);
	 //Nmodel=ConvertToModelFrg(seq,OUT,iter,in,&mod,model,Nmodel,false);
	}
  printf("##Nmodel= %d\n",Nmodel);

 }




 return Nmodel;
}
