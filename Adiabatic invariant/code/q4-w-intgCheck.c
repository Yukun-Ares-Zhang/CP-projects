#include <stdio.h>
#include <math.h>
#include <limits.h>
#define MAX 8
#define EPS_ITG 1e-6
#define EPS_ROOT 1e-14
#define PI 3.1415926535897932

double N=256;
double dt=1e-2;
double Y_0=-2.1;
double PY_0=0;
double t,xt,gt;
double Ymin;

double g(double t){
    double g=2*cos(2*PI*t/N);
    return g;
}

double x(double t){
    double x=2*sin(2*PI*t/N);
    return x;
}

double eq1(double y){
    double eq1=y/sqrt(xt*xt+y*y)-y-gt;
    return eq1;
}

double pySquare(double y){
    double pySquare=PY_0*PY_0+Y_0*Y_0-y*y+2*(sqrt(xt*xt+y*y)-sqrt(xt*xt+Y_0*Y_0))+2*(Y_0-y)*gt;
    return pySquare;
}

double V(double y){//给定y返回体系势能
    double potential=0.5*pow(sqrt(xt*xt+y*y)-1,2)+y*gt;
    return potential;
}

int isInside(double a,double b,double x){//判断x是否在a,b之间
    if((a<x&&x<b)||(b<x&&x<a))return 1;
    else return 0;
}

double DekkerSolveEq1(double a, double b){//若区间错误返回-404404
    double temp,m,s,fa,fb;
    if(eq1(a)*eq1(b)>0)return -404404;
    if(fabs(eq1(b))>fabs(eq1(a))){
        temp=a;
        a=b;
        b=temp;
    }
    fa=eq1(a);
    fb=eq1(b);

    while(fb!=0&&fabs(a-b)>EPS_ROOT){
        m=(a+b)/2;
        if(fa!=fb){
            s=(a*fb-b*fa)/(fb-fa);
            if(isInside(m,b,s)&&fabs(eq1(s))<fabs(fb)){
                if(fb*eq1(s)<0){
                    a=b;fa=fb;
                    b=s;fb=eq1(s);
                }
                else{
                    b=s;fb=eq1(s);
                    if(fb*eq1(m)<0){
                        a=m;fa=eq1(m);
                    }
                }

            }
            else{
                pos_1:
                if(eq1(m)*fb<0){
                    a=m;
                    fa=eq1(m);
                }
                else{
                    b=m;
                    fb=eq1(m);
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

double DekkerPySquare(double a, double b){//若区间错误返回-404404
    double temp,m,s,fa,fb;
    if(pySquare(a)*pySquare(b)>0)return -404404;
    if(fabs(pySquare(b))>fabs(pySquare(a))){
        temp=a;
        a=b;
        b=temp;
    }
    fa=pySquare(a);
    fb=pySquare(b);

    while(fb!=0&&fabs(a-b)>EPS_ROOT){
        m=(a+b)/2;
        if(fa!=fb){
            s=(a*fb-b*fa)/(fb-fa);
            if(isInside(m,b,s)&&fabs(pySquare(s))<fabs(fb)){
                if(fb*pySquare(s)<0){
                    a=b;fa=fb;
                    b=s;fb=pySquare(s);
                }
                else{
                    b=s;fb=pySquare(s);
                    if(fb*pySquare(m)<0){
                        a=m;fa=pySquare(m);
                    }
                }
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
    if(fb>=0)return b;
    else return a;
    
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

double integratePeriod(double y1,double y2){
    double T=0,T0=-1;
    double yi;
    int n,i;
    for ( n = 1; fabs(T-T0)>EPS_ITG; n++)
    {
        double S=0;
        for ( i = 0; i < n; i++)
        {
            yi=(y1+y2)/2+(y2-y1)/2*cos((2*i+1)*PI/(2*n));
            S+=sqrt((y2-yi)*(yi-y1))/py(yi);
        }
        T0=T;
        T=S*PI/n;
    }
    return T;
}

double getW(void){
    double H=0.5*PY_0*PY_0+0.5*pow(sqrt(xt*xt+Y_0*Y_0)-1,2)+gt*Y_0;
    Ymin=DekkerSolveEq1(-3,3);
    if(Ymin<-100)return -404404;
    if(H<V(Ymin))return -404;

    double y1,y2;//相空间积分上下限
    y1=DekkerPySquare(-8,Ymin);
    y2=DekkerPySquare(Ymin,8);
    printf("%.15lf %.15lf\n",pySquare(y1),pySquare(y2));
    double T;//周期
    T=2*integratePeriod(y1,y2);
    if(T==0)return INT_MAX;
    double w=2*PI/T;
    return w;
}


int main(){
    double y1,y2;
    y1=-3.9021312790940375;
    y2=-2.1000000011249784;
    Y_0=-2.100163386565657;
    PY_0=0.156731404156210;
    t=5981*dt;
    xt=x(t);
    gt=g(t);

    double w=getW();
    printf("%lf\n",w);
    return 0;
}