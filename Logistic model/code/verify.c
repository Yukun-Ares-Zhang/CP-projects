#include <stdio.h>
#include <math.h>
#define e 2.718281828459045
double r=1;
long long N1=10000000;
long long N2=2048;
int itv=512;

double f(double x){
    double y=r*sin(x);
    return y;
}
int main(){
    FILE* fdes=fopen("C://C_files//HW1//periodNum.txt","w");
    if(fdes==NULL){
        return 0;
    }
    for ( r = 2.7192965; r < 2.7192975; r+=1e-7)
    {
        fprintf(fdes,"r=%.7lf\n",r);
       double x=0.3;
        long long i;
        for ( i = 0; i < N1; i++)
        {
            x=f(x);
        }
        for(i=0;i<N2;i++){
            x=f(x);
            if(i%itv==0){
                fprintf(fdes,"%.10lf ",x);
            }
        }  
        fprintf(fdes,"\n");
    }
    
    return 0;
}

