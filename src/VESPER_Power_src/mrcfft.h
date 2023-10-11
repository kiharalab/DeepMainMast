

typedef struct{
 int t[3];
 float r[3];
 float q[4];
 int code;
 float sco;
} TBL;

bool SearchMAPfft(MRC *,MRC *,float);
bool SearchMAPfftMT(MRC *,MRC *,MRC *,int,bool,MRC *,MRC *);
bool SimpleScoring(MRC *,MRC *);
bool SearchMAPfftMT_OVCC(MRC *,MRC *,int,int,bool); //Overlap or CCC
float GetScore(MRC *, MRC *, int [3]);
float GetScore2(MRC *, MRC *, int [3], float [5]);

