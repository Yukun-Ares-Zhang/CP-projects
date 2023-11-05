#include "matrix.h"

#define N 1024

int main(){
    FILE *fdata=fopen("C:/C_files/HW3/q1-N=1024.txt","w");
    SqMatrix* A;
    A=createZeroSqMatrix(N);
    //构建矩阵
    A->el[A->n-1][0]--;;A->el[0][A->n-1]--;A->el[0][0]=2;
    int i;
    for ( i = 1; i < A->n; i++)
    {
        A->el[i][i]=2;
        A->el[i-1][i]--;
        A->el[i][i-1]--;
    }
    //求特征值
    double* test=eigenvalue_sym_sq(A);
    if(test==NULL)return 0;
    for ( i = 1; i <= A->n; i++)
    {
        fprintf(fdata,"%lf %lf\n",test[i],2*(1-cos(2*PI*i/N)));
    }
    
    return 0;
}