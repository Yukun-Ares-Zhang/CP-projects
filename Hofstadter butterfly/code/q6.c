#include "matrix.h"
#define ITV 1e-3
#define N 720
double a=ITV;
double r=2;

int main(){
    FILE *fdata=fopen("C:/C_files/HW3/q6.txt","w");
    for(a=ITV;a<1;a+=ITV)
    {
    SqMatrix* A;
    A=createZeroSqMatrix(N);
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
        fprintf(fdata,"%lf  %lf\n",a,test[i]);
    }
     printf("a=%lf finished...\n",a);
    }
    fclose(fdata);
    return 0;
}