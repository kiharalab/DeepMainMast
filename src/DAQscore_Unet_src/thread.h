

typedef struct{
 char id[100];
 int len;
 int  seq_code[10000];
 int  ss[10000];
 char seq_txt[10000];
 float Pss[10000][3];
 int chain_id[10000];
}  SEQ;

typedef struct{
 double **cd;//xyz
 int seq_pos;
 int N;//length
 double score;
} SEQ_FRAG;


bool readseq(SEQ *,char *);
int splitseq(SEQ *,SEQFG *,int);
