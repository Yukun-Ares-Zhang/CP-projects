#include "matrix.h"

double a=1./2;
double r=2;
int n=720;//n>=2

int main(){
    FILE *fdata=fopen("C:/C_files/HW3/q2-N=720.txt","w");
    SqMatrix* A;
    A=createZeroSqMatrix(n);
    A->el[A->n-1][0]=-1;A->el[0][A->n-1]=-1;
    A->el[0][0]=2-r;
    int i;
    for ( i = 1; i < A->n ; i++)
    {
        A->el[i-1][i]=A->el[i][i-1]=-1;
        A->el[i][i]=2+r*cos(2*PI*(i+1)*a);
    }
    double *test=eigenvalue_sym_sq(A);
    if(test==NULL)return 0;
    for(i=1;i<=A->n;i++){
        fprintf(fdata,"%lf  %lf %lf\n",2-test[i],(2-test[i])*(2-test[i]),r*r+2*(1+cos(4*PI*i/n)));
    }
    return 0;
}