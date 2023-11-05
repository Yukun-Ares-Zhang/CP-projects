#include <stdio.h>
#include <math.h>
#define EPS 1e-10//收敛精度
#define MAX 10000
#define itv 1e-2
#define RMIN 3.4
#define RMAX 3.7
#define POWER_MAX 10

double r;
int TWO[POWER_MAX+1];

double f(double x){
    double y=r*x*(1-x);
    return y;
}
void creatPower2(void){
    int i;
    TWO[0]=1;
    for ( i = 0; i < POWER_MAX; i++)
    {
        TWO[i+1]=2*TWO[i];
    }
}

int main(){
    FILE* fdes=fopen("C://C_files//HW1//periodNum.txt","w");
    if(fdes==NULL){
        return 0;
    }
    creatPower2();
    for ( r = RMIN; r < RMAX; r+=itv)
    {
        int count=0;
        double x[POWER_MAX+2]={0.9};
        int i,j;
        double y=x[0];
        double x00=-9;
        for ( i = 0; i < TWO[POWER_MAX]; i++)
        {
            y=f(y);
            for ( j = 0; j < POWER_MAX+1; j++){
                if(j==POWER_MAX){
                    if(i==TWO[POWER_MAX]-1)x[j+1]=y;
                }
                else if(i==TWO[POWER_MAX]-TWO[POWER_MAX-j-1]-1)x[j+1]=y;
            }
        }
        

        while (fabs(x00-x[0])>EPS||fabs(x[0]-x[POWER_MAX+1])>EPS){
            y=x[0];
            x00=x[0];
            x[0]=x[POWER_MAX+1];
            for ( i = 0; i < TWO[POWER_MAX]; i++)
            {
                y=f(y);
                for ( j = 0; j < POWER_MAX+1; j++){
                    if(j==POWER_MAX){
                        if(i==TWO[POWER_MAX]-1)x[j+1]=y;
                    }
                    else if(i==TWO[POWER_MAX]-TWO[POWER_MAX-j-1]-1)x[j+1]=y;
                }
            }
            
        }  
        fprintf(fdes,"%lf ",r);
        int noOutput=1;
        if(fabs(x[0]-x[POWER_MAX+1])>EPS){
            fprintf(fdes,">%d",TWO[POWER_MAX]);noOutput=0;
        }
        else{
            for ( j = 0; j < POWER_MAX+1; j++){
                if(fabs(x[POWER_MAX+1]-x[POWER_MAX-j])<EPS){
                    fprintf(fdes,"%d\n",TWO[j]);
                    noOutput=0;
                    break;
                }
            }
        }
        if (noOutput)
        {
           fprintf(fdes,"\n");
        }
        for ( j = 0; j < POWER_MAX+2; j++)
        {
            fprintf(fdes,"%.10lf ",x[j]);
        }
        fprintf(fdes,"\n");
    }
    return 0;
}