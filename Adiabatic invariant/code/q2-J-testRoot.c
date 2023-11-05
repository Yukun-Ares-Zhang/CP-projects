#include <stdio.h>
#include <math.h>
#define EPS_ITG 1e-10
#define EPS_ROOT 1e-16
#define PI 3.1415926535897932
double N=256;
double dt=1e-3;
double Y_0=0.1;
double PY_0=0;
double Y_o,Y_n,PY_o,PY_n;
double xt;
double x(double t){
    if(t<0)return -404404;
    return 2-t/N;
}

double pySquare(double y){
    double pySquare=Y_0*Y_0-y*y+2*(sqrt(xt*xt+y*y)-sqrt(xt*xt+Y_0*Y_0))+PY_0*PY_0;
    return pySquare;
}

int isInside(double a,double b,double x){//判断x是否在a,b之间
    if((a<x&&x<b)||(b<x&&x<a))return 1;
    else return 0;
}

double DekkerPySquare(double a, double b){//返回正根，若区间错误返回-1
    double temp,m,s,fa,fb;
    if(pySquare(a)*pySquare(b)>0)return -1;
    if(fabs(pySquare(b))>fabs(pySquare(a))){
        temp=a;
        a=b;
        b=temp;
    }
    fa=pySquare(a);
    fb=pySquare(b);

    while(fabs(fb)>EPS_ROOT){
        m=(a+b)/2;
        if(fa!=fb){
            s=(a*fb-b*fa)/(fb-fa);
            if(isInside(m,b,s)&&fabs(pySquare(s))<fabs(fb)){
                if(fb*pySquare(s)<0){
                    a=b;fa=fb;
                }
                b=s;fb=pySquare(s);
            }
            else{
                pos_1:
                if(pySquare(m)*fb<0){
                    a=m;
                    fa=pySquare(m);
                }
                else{
                    b=m;
                    fb=pySquare(m);
                }
            }
        }
        else goto pos_1;

        if(fabs(fb)>fabs(fa)){
            temp=b;b=a;a=temp;
            temp=fb;fb=fa;fa=temp;
        }
    }
    return b;
}
double py(double y){
    double py_pingfang=pySquare(y);
    double py;
    if(py_pingfang<0)return 0;
    else{
        py=sqrt(py_pingfang);
        return py;
    }
}
double integrate_phaseSpace(double y1,double y2){
    double h=fabs(y2-y1);
    double T=h/2*(py(y1)+py(y2));
    double T0=T-1;
    long long n=2;
    while (fabs(T-T0)>EPS_ITG){
        T0=T;
        double S=0;
        double dh=h/n;
        long long i;
        for ( i = 0; 2*i+1 < n ; i++)
        {
            S+=py(y1+(2*i+1)*dh);
        }
        T=S*dh+0.5*T0;
        n*=2;
    }
    return T;
}
double getJ(){
    double y0=fabs(Y_0),py0=fabs(PY_0);
    double yc=DekkerPySquare(y0,2.1);
    double J;
    if(yc<0)return -404404;
    if(y0*y0+2*(xt-sqrt(xt*xt+y0*y0))+py0*py0>0){
        J=2/PI*(integrate_phaseSpace(0,yc));
    }
    else{
        double ym=DekkerPySquare(0,y0);
        if(ym<0)return -404404;
        J=1/PI*(integrate_phaseSpace(ym,yc));
    }
    return J;
}

int main(){
    xt=0.975859;
    Y_0=0.268230479181774;
    PY_0=-0.000005846044208;
    double J;
    J=getJ();
    printf("%lf\n",J);
}