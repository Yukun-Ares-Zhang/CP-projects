#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS 1e-8
typedef struct SqMetrix{
    int n;
    double** el;
}SqMetrix;

int sgn(double x){//符号函数
    if(x==0)return 0;
    if(x>0)return 1;
    return -1;
}

double** createNewMatrix(int m,int n){
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
void printMetrix(SqMetrix* A){
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
void printArray(double *a, int ini,int fin){
    int i;
    for ( i = ini; i <= fin; i++)
    {
        printf("%lf",a[i]);
        if(i!=fin)printf("\t");
        else printf("\n");
    }
}
int is_sym_sq(SqMetrix* A){//判断方阵是否对称
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
int main(){//求出实对称方阵的特征值组
    int i,j,k;
    SqMetrix A;
    A.n=3;
    A.el=createNewMatrix(A.n,A.n);
    for ( i = 0; i < A.n; i++)
    {
        A.el[i][i]=i+2;
        if(i!=A.n-1){
            A.el[i+1][i]=A.el[i][i+1]=1;
        }
    }
    
    //Householder变换化为三对角阵
    double *al=(double *)malloc((A.n+1)*sizeof(double));
    double *b=(double *)malloc((A.n+1)*sizeof(double));
    double *au=(double *)malloc((A.n+1)*sizeof(double));
    double *d=(double *)malloc(A.n*sizeof(double));

    for ( i = 0; i < A.n-2; i++)
    {
        for ( j = i+2; j < A.n; j++)
        {
            if(A.el[j][i]!=0) goto Householder;
        }
        continue;

        Householder:
        {
        double sigma=0,beta;
        for ( j = i+1; j < A.n; j++)
        {
            sigma+=A.el[j][i]*A.el[j][i];
        }
        sigma=sqrt(sigma);
        sigma*=sgn(A.el[i+1][i]);    
        beta=sigma*(sigma+A.el[i+1][i]);
        A.el[i+1][i]+=sigma;

        SqMetrix R;//生成反射矩阵
        R.n=A.n-1-i;
        R.el=createNewMatrix(R.n,R.n);
        for ( j = 0; j < R.n; j++)
        {
            R.el[j][j]++;
        }
        for ( j = 0; j < R.n; j++)
        {
            for ( k = 0; k < R.n; k++)
            {
                R.el[j][k]-=(A.el[i+1+j][i]*A.el[i+1+k][i])/beta;
            }
        }

        A.el[i+1][i]=-sigma;A.el[i][i+1]=-sigma;
        for ( j = i+2; j < A.n; j++)
        {
            A.el[j][i]=0;A.el[i][j]=0;
        }

        //计算RAR，并替代A
        double** C=createNewMatrix(R.n,R.n);//中间矩阵储存RA
        for ( j = 0; j < R.n; j++)
        {
            for ( k= 0; k < R.n; k++)
            {
                int m;
                for ( m = 0; m < R.n; m++)
                {
                    C[j][k]+=R.el[j][m]*A.el[i+1+m][i+1+k];
                }
            }
        }
        for ( j = 0; j < R.n; j++)
        {
            for ( k= 0; k < R.n; k++)
            {
                A.el[i+1+j][i+1+k]=0;
                int m;
                for ( m = 0; m < R.n; m++)
                {
                    A.el[i+1+j][i+1+k]+=C[j][m]*R.el[m][k];
                }
            }
        }
        }
    }
    b[1]=A.el[0][0];
    for ( i = 2; i <= A.n; i++)
    {
        al[i]=A.el[i-1][i-2];
        b[i]=A.el[i-1][i-1];
        au[i-1]=A.el[i-2][i-1];
    }

    
    //QR迭代求特征值
    int n=A.n;
    printArray(au,1,A.n-1);
    printArray(b,1,A.n);
    printArray(al,2,A.n);
    printf("\n");
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
        printArray(d,1,A.n-2);
        printArray(au,1,A.n-1);
        printArray(b,1,A.n);
        printArray(al,2,A.n);
        printf("\n");
        for ( k = 1; k < n; k++)
        {
            if(k!=1){
                temp1=au[k-1];temp2=d[k-1];
                au[k-1]=temp1*cosine[k]+temp2*sine[k];
                d[k-1]=-temp1*sine[k]+temp2*cosine[k];
            }
            temp1=b[k];temp2=au[k];
            b[k]=temp1*cosine[k]+temp2*sine[k];
            au[k]=-temp1*sine[k]+temp2*cosine[k];
            temp1=al[k+1];temp2=b[k+1];
            al[k+1]=temp1*cosine[k]+temp2*sine[k];
            b[k+1]=-temp1*sine[k]+temp2*cosine[k];
            if(k!=n-1)al[k+2]*=cosine[k];

            b[k]+=s;
        }
        b[n]+=s;
        if(fabs(al[n])<EPS){
            al[n]=0;
            n--;
        }
        printArray(d,1,A.n-2);
        printArray(au,1,A.n-1);
        printArray(b,1,A.n);
        printArray(al,2,A.n);
        printf("\n");
    }
    for ( i = 1; i <= A.n; i++)
    {
        printf("%lf ",b[i]);
    }
    return 0;
}
