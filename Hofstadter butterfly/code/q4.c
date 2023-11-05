#include "matrix.h"

double a=1./4;
double r=-2;
int n=720;//n>=2

int main(){
    FILE *fdata=fopen("C:/C_files/HW3/q4-explore-r=-2.txt","w");
    SqMatrix* A;
    A=createZeroSqMatrix(n);
    A->el[A->n-1][0]=-1;A->el[0][A->n-1]=-1;
    A->el[0][0]=2+r*cos(2*PI*a);
    int i;
    for ( i = 1; i < A->n ; i++)
    {
        A->el[i-1][i]=A->el[i][i-1]=-1;
        A->el[i][i]=2+r*cos(2*PI*(i+1)*a);
    }
    double *test=eigenvalue_sym_sq(A);
    if(test==NULL)return 0;
    for(i=1;i<=A->n;i++){
        fprintf(fdata,"%lf\n",2-test[i]);
    }
    return 0;
}