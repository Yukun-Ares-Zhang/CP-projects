#ifndef __MATRIX_H__
#define __MATRIX_H__

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifndef EPS
#define EPS 1e-7
#endif

#ifndef TOL
#define TOL 1e-10
#endif

#ifndef PI
#define PI 3.1415926535897932
#endif

typedef struct SqMatrix{
    int n;
    double** el;
}SqMatrix;

int sgn(double x){//符号函数
    if(x>=0)return 1;
    return -1;
}
void printMatrix(SqMatrix* A){
    int i,j;
    for ( i = 0; i < A->n; i++)
    {
        for ( j = 0; j < A->n; j++)
        {
            printf("%lf",A->el[i][j]);
            if(j!=A->n-1)printf("\t");
            else printf("\n");
        }
    }
}
double** createNew2DArray(int m,int n){
    int i,j;
    double **a=(double **)malloc(m*sizeof(double *));
    for ( i = 0; i < m; i++)
    {
        a[i]=(double *)malloc(n*sizeof(double));
    }
    for ( i = 0; i < m; i++)
    {
        for ( j = 0; j < n; j++)
        {
            a[i][j]=0;
        }
    }
    return a;
}
SqMatrix* createZeroSqMatrix(int n){
    SqMatrix* p=(SqMatrix*)malloc(sizeof(SqMatrix));
    p->n=n;
    p->el=createNew2DArray(n,n);
    return p;
}
int is_sym_sq(SqMatrix* A){//判断方阵是否对称
    if(A->el==NULL)return 1;
    int i,j;
    for ( i = 0; i < A->n; i++)
    {
        for ( j = i+1; j < A->n; j++)
        {
            if(A->el[i][j]!=A->el[j][i])return 0;
        }
    }
    return 1;
}

double* eigenvalue_sym_sq(SqMatrix* A){//求出实对称方阵的特征值组
    if(!is_sym_sq(A)||A->n<=0)return NULL;
    //Householder变换化为三对角阵
    double *al=(double *)malloc((A->n+1)*sizeof(double));
    double *b=(double *)malloc((A->n+1)*sizeof(double));
    double *au=(double *)malloc((A->n+1)*sizeof(double));
    double *d=(double *)malloc(A->n*sizeof(double));
    if(A->n==1){
        b[1]=A->el[0][0];return b;
    }
    if(A->n==2){
        double b1=A->el[0][0],b2=A->el[1][1],a=A->el[0][1];
        b[1]=(b1+b2+sqrt(pow((b1-b2),2)+4*a*a))/2;
        b[2]=(b1+b2-sqrt(pow((b1-b2),2)+4*a*a))/2;
        return b;
    }
    int i,j,k;
    //printf("entered Householder\n");
    for ( i = 0; i < A->n-2; i++)
    {
        for ( j = i+2; j < A->n; j++)
        {
            if(fabs(A->el[j][i])>EPS) {
                //printf("A(%d,%d)=%.15lf   ",j,i,A->el[j][i]);
                goto Householder;
            }
        }
        continue;

        Householder:
        {
        //printf("start a householder transformation : %d\n",i);
        double sigma=0,beta,uAu;
        for ( j = i+1; j < A->n; j++)
        {
            sigma+=A->el[j][i]*A->el[j][i];
        }
        sigma=sqrt(sigma);
        sigma*=sgn(A->el[i+1][i]);    
        beta=sigma*(sigma+A->el[i+1][i]);
        A->el[i+1][i]+=sigma;//更新u
        //计算Au,储存在uT中
        for ( j = 1; j < A->n-i; j++)
        {
            A->el[i][i+j]=0;
            for ( k = 1; k < A->n-i; k++)
            {
                A->el[i][i+j]+=A->el[i+j][i+k]*A->el[i+k][i];
            }
        }
        //计算uTAu
        uAu=0;
        for ( j = 1; j < A->n-i; j++)
        {
            uAu+=A->el[i+j][i]*A->el[i][i+j];
        }
        //更新A
        for ( j = 1; j < A->n-i; j++)
        {
            for ( k = 1; k < A->n-i; k++)
            {
                A->el[i+j][i+k]-=1/beta*(A->el[i+j][i]*A->el[i][i+k]+A->el[i][i+j]*A->el[i+k][i]);
                A->el[i+j][i+k]+=uAu/(beta*beta)*A->el[i+j][i]*A->el[i+k][i];
            }
        }
        A->el[i+1][i]=A->el[i][i+1]=-sigma;
        for ( j = 2; j < A->n-i; j++)
        {
            A->el[i+j][i]=A->el[i][i+j]=0;
        }
        //printf("%.20lf %.20lf\n",sigma,beta);     
        }
    }
    b[1]=A->el[0][0];
    for ( i = 2; i <= A->n; i++)
    {
        al[i]=A->el[i-1][i-2];
        b[i]=A->el[i-1][i-1];
        au[i-1]=A->el[i-2][i-1];
    }
    //QR迭代求特征值
    int n=A->n;
    double* cosine=(double *)malloc(n*sizeof(double));
    double* sine=(double *)malloc(n*sizeof(double));
    while(n>=2){
        double s=b[n];
        double temp1,temp2;
        b[1]-=s;
        for ( k = 1; k < n; k++)
        {
            b[k+1]-=s;
            double r=sqrt(b[k]*b[k]+al[k+1]*al[k+1]);
            cosine[k]=b[k]/r;sine[k]=al[k+1]/r;
            //左变换
            b[k]=r;al[k+1]=0;
            temp1=au[k],temp2=b[k+1];
            au[k]=temp1*cosine[k]+temp2*sine[k];
            b[k+1]=-temp1*sine[k]+temp2*cosine[k];
            if(k!=n-1){d[k]=au[k+1]*sine[k];au[k+1]*=cosine[k];}
        }
        for ( k = 1; k < n; k++)
        {
            if(k!=1){
                temp1=au[k-1];temp2=d[k-1];
                au[k-1]=temp1*cosine[k]+temp2*sine[k];
                d[k-1]=-temp1*sine[k]+temp2*cosine[k];
            }
            temp1=b[k];temp2=au[k];
            b[k]=temp1*cosine[k]+temp2*sine[k];
            au[k]=temp1*(-sine[k])+temp2*cosine[k];
            temp1=al[k+1];temp2=b[k+1];
            al[k+1]=temp1*cosine[k]+temp2*sine[k];
            b[k+1]=temp1*(-sine[k])+temp2*cosine[k];
            if(k!=n-1)al[k+2]*=cosine[k];

            b[k]+=s;
        }
        b[n]+=s;
        if(fabs(al[n])<EPS){
            al[n]=0;
            //printf("%d\n",n);
            n--;
        }
    }
    return b;
}


typedef struct vector{
    int n;
    double* el;
}vector;

typedef struct pair{
    int original;
    int current;
}pair;


pair* initialize_pair(int n){
    pair* p;
    p=(pair*)malloc(sizeof(pair)*n);
    if(p==NULL)return NULL;
    int i;
    for(i=0;i<n;i++){
        p[i].original=p[i].current=i;
    }
    return p;
}

vector* createNewVector(int n){
    vector* pv;
    pv=(vector*)malloc(sizeof(vector));
    if(pv==NULL)return NULL;
    pv->n=n;
    pv->el=(double*)malloc(sizeof(double)*n);
    if(pv->el==NULL)return NULL;
    int i;
    for ( i = 0; i < n; i++)
    {
        pv->el[i]=0;
    }
    return pv;
}

int copy_vector(vector* src,vector* des){
    if(src->n!=des->n)return 0;
    int i;
    for(i=0;i<src->n;i++){
        des->el[i]=src->el[i];
    }
    return 1;
}
int copy_SqMatrix(SqMatrix* src,SqMatrix* des){
    if(src->n!=des->n)return 0;
    int i,j;
    for(i=0;i<src->n;i++){
        for(j=0;j<src->n;j++){
            des->el[i][j]=src->el[i][j];
        }
    }
    return 1;
}
int average_vector(vector* src1,vector* src2,vector* des){
    if(src1->n!=src2->n||des->n!=src2->n)return 0;
    int i;
    for(i=0;i<src1->n;i++){
        des->el[i]=(src1->el[i]+src2->el[i])/2;
    }
    return 1;
}

void print_vector(vector* b){
    int i;
    for ( i = 0; i < b->n; i++)
    {
        printf("%lf",b->el[i]);
        if(i!=b->n-1)printf(",");
    }
    printf("\n");
}

int cmp_pair_original(const void* a,const void* b){
    return (*(pair*)a).original>(*(pair*)b).original?1:-1;
}

int Solve_RankN_LinearEq(SqMatrix* A, vector* b){
    if(A->n<1||b->n<1||A->n!=b->n)return 0;
    int N=A->n;
    pair* p;
    p=initialize_pair(N);
    //LU分解 
    int r,i,k;
    for ( r = 0; r < N; r++)
    {
        for ( i = r; i < N; i++)
        {
            double s=0;
            for ( k = 0; k< r; k++)
            {
                s+=A->el[i][k]*A->el[k][r];
            }
            A->el[i][r]-=s;
        }
        
        int ir_max=r;
        for ( i = r+1; i < N; i++)
        {
            if(fabs(A->el[i][r])>fabs(A->el[ir_max][r]))ir_max=i;
        }

        if(fabs(A->el[r][r])<EPS){return 0;}
        if (ir_max!=r)
        {
            double temp;
            for ( i = 0; i < N; i++)
            {
                temp=A->el[r][i];
                A->el[r][i]=A->el[ir_max][i];
                A->el[ir_max][i]=temp;
            }
            temp=p[r].original;
            p[r].original=p[ir_max].original;
            p[ir_max].original=temp;    
        }
        
        for ( i = r+1; i < N; i++)
        {
            A->el[i][r]/=A->el[r][r];
            double s=0;
            for ( k = 0; k<r; k++)
            {
                s+=A->el[r][k]*A->el[k][i];
            }
            A->el[r][i]-=s;
        }
    }
    //解方程
    vector* Pb=createNewVector(N);
    copy_vector(b,Pb);
    for ( i = 0; i < N; i++)
    {
        b->el[i]=Pb->el[p[i].original];
    }
    double s;
    for ( i = 0; i < N; i++)
    {
        s=0;
       for ( k = 0; k < i; k++)
        {
            s+=A->el[i][k]*b->el[k];
        }
        b->el[i]-=s;
    }
    for ( i = N-1; i >=0; i--)
    {
        s=0;
        for ( k = i+1; k < N; k++)
        {
            s+=A->el[i][k]*b->el[k];
        }
        b->el[i]=(b->el[i]-s)/A->el[i][i];
    }
    
    return 1;
}


#endif