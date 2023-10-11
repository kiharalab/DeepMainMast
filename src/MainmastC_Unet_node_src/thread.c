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
#include "thread.h"

bool readseq(SEQ *seq,char *fname){
 int num=0;
 FILE *fp;
 char line[LIN],buf[LIN];
 int  ftype=-1;
 int i,len,Nch=0;
 char aa,ss;

 seq->len=0;
 //Detect File type
 if((fp=fopen(fname,"r"))==NULL)
  return true;
 while(fgets(line,LIN,fp)!=NULL){
  if(!strncmp(line,"#\tAA",4)){
   ftype=1;//*.spd3
   break;
  }
  if(!sscanf(line,">%s",seq->id)){
  //if(!strncmp(line,">",1)){
   ftype=0;//fasta format
   break;
  }
 }
 fclose(fp);
 //Read chain-break '/' or multiple chains

 Nch=0;
 if(ftype==0){
  puts("#Fasta file");
  if((fp=fopen(fname,"r"))==NULL)
   return true;
  while(fgets(line,LIN,fp)!=NULL){
   if(!strncmp(line,">",1)){
    //strncpy(seq->id,&line[1],len-1);
    //seq->id[len-2]='\0';
    printf("#ID= '%s'\n",seq->id);
    Nch++;
    if(Nch!=1){
     //Ignnore
     //add chain-break
     //seq->seq_txt[seq->len]='/';
     //seq->len++;
    }
   }else{
    len=strlen(line);
	for(int p=0;p<len;p++)
	 seq->chain_id[seq->len+p]=Nch-1;
    //strncpy(&seq->seq_txt[seq->len],line,len);
    strcpy(&seq->seq_txt[seq->len],line);
    seq->len+=len;
    if(seq->seq_txt[seq->len-1]=='\n')
     seq->len--;
   }
  }
  seq->seq_txt[seq->len]='\0';
  printf("#SEQ= %s\n",seq->seq_txt);
  fclose(fp);
  //convert to code
  int cid=0;
  for(i=0;i<seq->len;i++){
	 //seq->chain_id[i]=cid;
         seq->seq_code[i]=A2int(seq->seq_txt[i]);
         seq->ss[i]=3;//All C
	 seq->Pss[i][0]=seq->Pss[i][1]=seq->Pss[i][2]=0.00;//HEC
   //printf("%d %c %d\n",i,seq->seq_txt[i],seq->seq_code[i]);
  }
  printf("#CHAIN: ");
  for(int i=0;i<seq->len;i++)
   printf("%d",seq->chain_id[i]);
  puts("");
 }
if(ftype==1){
  puts("#spd3 file");
  if((fp=fopen(fname,"r"))==NULL)
   return true;
  while(fgets(line,LIN,fp)!=NULL){
   if(!strncmp(line,"#",1)){
     strncpy(seq->id,"Unknown",7);
     seq->id[8]='\0';
     printf("#ID= '%s'\n",seq->id);
     //chain-break
     Nch++;
     if(Nch!=1){
      seq->seq_code[seq->len]=-9;
      seq->seq_txt[seq->len]='/';
      seq->ss[seq->len]=-1;
      seq->len++;
     }
   }else{
    float Pc,Pe,Ph;
//#	AA	SS3	SS8	ASA 	HSEa-u	HSEa-d	CN13	theta	tau 	phi 	psi 	P(3-C)	P(3-E)	P(3-H)	P(8-C)	P(8-S)	P(8-T)	P(8-H)	P(8-G)	P(8-I)	P(8-E)	P(8-B)
    sscanf(line,"%d\t%c\t%c\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f",&i,&aa,&ss,buf,buf,buf,buf,buf,buf,buf,buf,buf,&Pc,&Pe,&Ph);
    seq->seq_code[seq->len]=A2int(aa);
    seq->seq_txt[seq->len]=aa;
    	seq->Pss[seq->len][0]=Ph*0.01;
    	seq->Pss[seq->len][1]=Pe*0.01;
    	seq->Pss[seq->len][2]=Pc*0.01;
	//printf("%f %f %f\n",Ph,Pe,Pc);
        //ss
        switch(ss){
                case 'C':seq->ss[seq->len]=3;
                 break;
                case 'H':seq->ss[seq->len]=1;
                 break;
                case 'E':seq->ss[seq->len]= 2;
                 break;
                default: seq->ss[seq->len]= -1;
        }

    //printf("%d %c %c %d %d\n",seq->len,seq->seq_txt[seq->len],ss,seq->seq_code[seq->len],seq->ss[seq->len]);
    seq->len++;
    
   }
  }
 fclose(fp);
  seq->seq_txt[seq->len]='\0';
  printf("#SEQ    '%s'\n",seq->seq_txt);
  printf("#CHAIN: ");
  for(int i=0;i<seq->len;i++)
   printf("%d",seq->chain_id[i]);
  puts("");
  for(int i=seq->len;i<seq->len+5;i++)
   seq->ss[seq->len]= -1;
 }
 if(ftype==-1){
  printf("Cannot detect file format..%s\n",fname);
  return true;
 }
 return false;
}

int splitseq(SEQ *seq,SEQFG *fg,int len){
 int n=0;
 int i,j,pos;
 int Nc=0;
 bool flag;
 //int flen=len+4;//???
 int flen=len;//???
 for(i=0;i<seq->len-flen+1;i++){
  if(seq->seq_code[i]==-9){
   Nc++;
   continue;
  }
  //check
        flag=false;
        for(j=0;j<flen;j++){
         if(seq->seq_code[i+j]==-9){//chain break
          flag=true;
          break;
         }
        }
        if(flag==true) continue;
  //input
  fg[n].ss=(int *)malloc(sizeof(int)*flen);
  fg[n].seq=(int *)malloc(sizeof(int)*flen);
  for(j=0;j<flen;j++){
   fg[n].seq[j]=seq->seq_code[i+j];
   fg[n].ss[j]=seq->ss[i+j];
   fg[n].Pss[j][0]=seq->Pss[i+j][0];
   fg[n].Pss[j][1]=seq->Pss[i+j][1];
   fg[n].Pss[j][2]=seq->Pss[i+j][2];
   printf("%3d",fg[n].seq[j]);
   printf("(%1d,[%2.0f,%2.0f,%2.0f])",fg[n].ss[j],fg[n].Pss[j][0]*10,fg[n].Pss[j][1]*10,fg[n].Pss[j][2]*10);
  }
  printf(" chain=%d\n",Nc);
  fg[n].cid=Nc;
  fg[n].l=flen;
  fg[n].pos=i;
  n++;
 }
 return n;
}


