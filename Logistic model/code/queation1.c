#include <stdio.h>
double r=0.5;//可调参数r
int N=10;//迭代次数
double x=0.1;//迭代初值
double f(double x){
    double y=r*x*(1-x);
    return y;
}
int main(){
    int i;
    for ( i = 0; i < N+1; i++)
    {
        printf("%lf\n",x);
        x=f(x);
    }
    return 0;    
}