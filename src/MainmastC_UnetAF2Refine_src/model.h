
typedef struct{
 NODE *Ca;
 int Nnode;
 int Npath,LongestPath,Lseq;
 int **Path;
 int Lpath[2000];//Max 2k
 int **Smtx;
} NODE_CSV;

//For each model
typedef struct{
 NODE Ca[5000];//Max 5k length
 int Start,End,Len;
 float RawScore,Score,Rmsd,Weight;
} RESULT_MODEL;

float Dist2(NODE *,NODE *);
bool ReadNodeCsv(char *,NODE_CSV *);
int GenerateFragmentsLocal(SEQ *,NODE_CSV *,PDB *,int,float,RESULT_MODEL **);
int GenerateFragments(SEQ *,NODE_CSV *,PDB *,NODE_CSV *,int,float,RESULT_MODEL **);

int AddModel(RESULT_MODEL *,RESULT_MODEL **, int,bool);

int Threading(SEQ *,NODE_CSV *,RESULT_MODEL **,int,int,bool);
