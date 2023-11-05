//本程序求一个r值下迭代序列第100001-100100项
#include <stdio.h>
#include <math.h>

double r=3.55;
long long N=100000;

double f(double x){
    double y=r*x*(1-x);
    return y;
}
int main(){
    double x=0.9;
    long long i;
    for ( i = 0; i < N; i++)
    {
        x=f(x);
    }
    for(i=0;i<100;i++){
        x=f(x);
        printf("%lf\n",x);
    } 
    return 0;
}