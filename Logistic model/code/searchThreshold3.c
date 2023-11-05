#include <stdio.h>
#include <math.h>
#define EPS1 1e-10//周期分点的精度
#define EPS2 1e-16//判断数列收敛的精度
//以下是需要手动设置的三参数
#define A 2.7192969
#define B 2.7192971
#define N 512

double f(double x, double r){
    double y=r*sin(x);
    return y;
}
int periodIsM(double r,int m){//对于参数r，判断序列是不是周期为m
    double x0=-9,x=0.3;//序列初值0.3，x0用于储存前一个x值
    int i;
    while ( fabs(x-x0)>EPS2){//迭代使序列收敛
        x0=x;
        for(i=0;i<2*m;i++){x=f(x,r);}
    }//由于序列周期只可能为m,m/2,2m，每隔2m取出的子序列一定收敛，那么每隔2m取一个值，当与前一个值做差小于EPS2时认为序列收敛
    double y=x;
    for ( i = 0; i < m; i++){x=f(x,r);}//再迭代m次看两个值是否相近，若在小于EPS1精度内相等，认为周期为m
    if(fabs(y-x)<EPS1){return 1;}
    else return 0;
}
double SearchThreshold(double a, double b,int m,int n, double tol){//寻找周期分点函数
    if(!(periodIsM(a,m)&&periodIsM(b,n))){//区间错误
        printf("FALSE PARAMETER!\n");
        return -404404;
    }
    double r;
    while(fabs(b-a)>tol){
        r=(a+b)/2;
        if(periodIsM(r,m)){a=r;}
        else{b=r;}
    }
    return (a+b)/2;
}
int main(){
    double rt;
    rt=SearchThreshold(A,B,N,2*N,EPS1);
    printf("%.10lf",rt);
    return 0;
}
