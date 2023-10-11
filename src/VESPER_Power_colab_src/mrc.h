typedef struct{
	char filename[LIN];
	float map_t;
	int Nthreads;
	int xdim,ydim,zdim;
	int ncstart,nrstart,nsstart;
	int mx,my,mz;
	float xlen,ylen,zlen;
	float alpha,beta,gamma;
	int mapc,mapr,maps;
	float dmin,dmax,dmean;
	int ispg;
	int nsymbt;
	float orgxyz[3];
	int NumVoxels;
	float *dens,*sco;
	float **vec;
	float widthx,widthy,widthz;
	unsigned int Nact;
	float cent[3];
	float dmax2,dsum,std,ave;
	float std_norm_ave;
} MRC;


typedef struct{
 float **cd,*dens;
 int **origrid;
 int Ncd,Nori;
 int *member;
 float *mask;
} POINTS;

typedef struct{
 float d;
 float dens;
 int id1,id2,eid;
 bool mst,local,keep;
} EDGE;

typedef struct{
 int N;
 EDGE *e[20];//!!!!! 10 is not enough
 int cid;//chain
} NODE;

typedef struct{
 int cut_id,add_id;
 float score;
 //float score;
} MOVE;



typedef struct{
 float len,bf_len;
 int Nnode,Ne;
 NODE *node;
 int Ntotal,Etotal;
 int St,Ed;
 bool *ActE;//Active edge
 bool *ActN;//Active node
 //bool *MstE://MST edge
 int *stock,Nstock;
 int *nextv;
 float *cost;
 int *cid;
 int *CutTbl,Ncut;//Edge table for cut;
 int *AddTbl,Nadd;//Edge table for add;

 MOVE *mv; //Movement
 int Nmv;
 int *Path,Lpath;
 float score;
} TREE;

typedef struct{
 EDGE *edge;//All edge data
 int Ne;
 bool **adj;
 NODE *node;//All connections
 int Nnode;
 int *cid;

 int Nt;
 EDGE **tree;

 //TREE mst;

} GRAPH;


bool readmrc(MRC *,char *);
bool ToCubic(MRC *);
bool upsampling(MRC *,float);
void out_situs(MRC *);
bool meanshift(MRC *,POINTS *);
bool fastLDP(MRC *,POINTS *);

bool fastVEC(MRC *,MRC *);
void SetUpVoxSize(MRC *,MRC *,float,float);
void SetUpComVoxSize(MRC *,MRC *,MRC *,MRC *,float,float,float);
void ShowVec(MRC *);
void ShowVec2(MRC *,MRC *,int [3]);
void ShowVec3(MRC *,MRC *,int [3]);
void ShowMRC(MRC *);

float meanshift_pos(MRC *,float [3]);
bool MergePoints(MRC *,POINTS *);
bool MergeLDP(MRC *,POINTS *);
void ShowModel(MRC *,POINTS *);
void ShowLDP(MRC *,POINTS *);
void ShowOri(MRC *,POINTS *);
void ShowGraph(GRAPH *);
void ShowTree(GRAPH *,TREE *);
void ShowPath(MRC *,POINTS *,GRAPH *,TREE *,int n);
void ShowPath2(MRC *,POINTS *,GRAPH *,TREE *,int n);

bool SetUpGraph(POINTS *, GRAPH *,MRC *,TREE *);
bool Tabu(GRAPH *,TREE *,TREE *);
bool InitTree(TREE *,bool);
bool ConstTree(TREE *,GRAPH *);
float QualityTree(GRAPH *, TREE *);
bool CopyTree(TREE *,TREE *,bool);
bool SetCutTbl(GRAPH *,TREE *,int *,int,float);
bool MoveTree(GRAPH *, TREE *, int,int);

int CutTree(GRAPH *, TREE *,float,int *,int);
int AddEdge(GRAPH *, TREE *);
int SplitChain(GRAPH *, TREE *,int);

float RandDouble();
int RandInt(int);
void ShuffleTbl(int *, int);
void ShuffleMv(MOVE *, int);

float Kendall(int *,int,int *,int);

bool PairExhaust(GRAPH *,TREE *,TREE *);
int ListStEd(GRAPH *,TREE *,int *);

float OptPath2points(GRAPH *,TREE *,TREE *,int,int,int);
bool get_path(GRAPH *,TREE *, int, int);
