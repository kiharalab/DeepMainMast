
#define GRID(N,a,b,c) c+N*(b+N*a)
#define GRID3D(N1,N2,a,b,c) c+N2*(b+N1*a)
#define GRID2D(N1,N2,a,b) (b)+(N2)*(a)


typedef struct{ 
	float sco; 
	short int poi;
	short int pre_align1,pre_align2;
} DPMTX;


float dp(DPMTX *,float *,float ,float ,int *,int,int *,int ,int *,int *);
float dp_fast(DPMTX *,float *,float *,float ,float ,int *,int,int *,int ,int *,int *,bool);
float dp_local(DPMTX *,float *,float *,float ,float ,int *,int,int *,int ,int *,int *,bool);

