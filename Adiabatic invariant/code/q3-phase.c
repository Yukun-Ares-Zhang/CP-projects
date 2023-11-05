#include <stdio.h>
#include <math.h>
#define MAX 8//求根的经验上界
#define EPS_ITG 1e-8//积分特征精度
#define EPS_ROOT 1e-14//求根特征精度
#define PI 3.1415926535897932
#define X 0.2
#define Y_THRESHOLD 0.277417932961286//Y临界，与X取值有关
#define G_THRESHOLD 0.53375//g临界，与X取值有关
double N=16;
double dt=5e-4;
double Y_0=-2;//储存初值条件y0
double PY_0=0;//储存初值条件py0
double t;
double Ymin1,Ymax,Ymin2;//储存三极值时的极值点
double Ymin;//储存单极值时的极值点

double g(double t){//给定t返回g
    double g=2*cos(2*PI*t/N);
    return g;
}

double eq1(double y){//求势函数极值点需要使用的方程
    double eq1=y/sqrt(X*X+y*y)-y-g(t);
    return eq1;
}

double pySquare(double y){//py的平方
    double pySquare=PY_0*PY_0+Y_0*Y_0-y*y+2*(sqrt(X*X+y*y)-sqrt(X*X+Y_0*Y_0))+2*(Y_0-y)*g(t);
    return pySquare;
}

double V(double y){//给定y返回体系势能
    double potential=0.5*pow(sqrt(X*X+y*y)-1,2)+y*g(t);
    return potential;
}
int isInside(double a,double b,double x){//判断x是否在a,b之间
    if((a<x&&x<b)||(b<x&&x<a))return 1;
    else return 0;
}

double DekkerExtreme(double a, double b){//若区间错误返回-404404
    double temp,m,s,fa,fb;
    if(eq1(a)*eq1(b)>0)return -404404;
    if(fabs(eq1(b))>fabs(eq1(a))){
        temp=a;
        a=b;
        b=temp;
    }
    fa=eq1(a);
    fb=eq1(b);

    while(fabs(a-b)>EPS_ROOT){
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

    while(fabs(a-b)>EPS_ROOT){
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
    double H=0.5*PY_0*PY_0+0.5*pow(sqrt(X*X+Y_0*Y_0)-1,2)+Y_0*g(t);
    double gt=g(t);
    double y1,y2;
    double J;
    if(fabs(gt)<0.53375){//势函数双极值
        Ymin1=DekkerExtreme(Y_THRESHOLD,1.6);
        Ymin2=DekkerExtreme(-1.6,-Y_THRESHOLD);
        if(gt>0)Ymax=DekkerExtreme(0,Y_THRESHOLD);
        else Ymax=DekkerExtreme(-Y_THRESHOLD,0);

        if(Ymin1<-100||Ymin2<-100||Ymax<-100)return -404404;

        if(H>V(Ymax)){
            y1=DekkerPySquare(-MAX,Ymin2);
            y2=DekkerPySquare(Ymin1,MAX);
        }
        else{
            if(Y_0>Ymax){
                y1=DekkerPySquare(Ymax,Ymin1);
                y2=DekkerPySquare(Ymin1,MAX);
            }
            else{
                y1=DekkerPySquare(-MAX,Ymin2);
                y2=DekkerPySquare(Ymin2,Ymax);
            }
        }
    }
    else{//势函数单极值
        if(gt>0){
            Ymin=DekkerExtreme(-3,0);
        }
        else Ymin=DekkerExtreme(0,3);

        if(Ymin<-100)return -404404;
        
        if(H>V(Ymin)){
            y1=DekkerPySquare(-MAX,Ymin);
            y2=DekkerPySquare(Ymin,MAX);
        }
        else return 0;
    }
    if(y1<-100||y2<-100)return -404404;

    J=1/PI*integrate_phaseSpace(y1,y2);

    return J;
}
int main(){
    FILE *fdata=fopen("C://C_files//HW2//q3-2.txt","w");
    if(fdata==NULL){
        return 0;
    }
    double y_o,y_n,py_o,py_n;
    y_n=Y_0;
    py_n=PY_0+dt/2*((Y_0+dt/4*PY_0)*(-1+1/sqrt(X*X+pow(Y_0+dt/4*PY_0,2)))-g(dt/4));
    t=0;
    double J;
    J=getJ();
    fprintf(fdata,"%.15lf  %.15lf %lf %lf\n",Y_0,PY_0,g(t),J);
    printf("%lf\n",g(t));//控制台输出一个g值，方便查看进度
    
    
    for(t=dt;t<=N/2;t+=dt){//蛙跳法
        y_o=y_n;
        py_o=py_n;
        y_n=y_o+dt*py_o;
        py_n=py_o+dt*(y_n*(-1+1/sqrt(X*X+y_n*y_n))-g(t));
        Y_0=y_n;
        PY_0=(py_n+py_o)/2;
        J=getJ();
        fprintf(fdata,"%.15lf  %.15lf %lf %lf\n",Y_0,PY_0,g(t),J);
        printf("%lf\n",g(t));//控制台输出一个g值，方便查看进度
    }
    fclose(fdata);
    return 0;
}