#include <stdio.h>
#include <math.h>
#include <limits.h>
#define MAX 8
#define EPS_ITG 1e-8
#define EPS_ROOT 1e-15
#define PI 3.1415926535897932

double N=256;
double dt=1e-2;
double Y_0=0.52;
double PY_0=0;
double t,xt,gt;
double Ymin1,Ymin2,Ymax;
double Y_threshold;

double g(double t){
    double g=0.3*cos(2*PI*t/N);
    return g;
}

double x(double t){
    double x=0.3*sin(2*PI*t/N);
    return x;
}

double f(double y){
    double f=pow(y*y+xt*xt,3)-pow(xt,4);
    return f;
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

double DekkerRoot(double a, double b,double(*f)(double y)){//若区间错误返回-404404
    double temp,m,s,fa,fb;
    if((*f)(a)*(*f)(b)>0)return -404404;
    if(fabs((*f)(b))>fabs((*f)(a))){
        temp=a;
        a=b;
        b=temp;
    }
    fa=(*f)(a);
    fb=(*f)(b);

    while(fb!=0&&fabs(a-b)>EPS_ROOT){
        m=(a+b)/2;
        if(fa!=fb){
            s=(a*fb-b*fa)/(fb-fa);
            if(isInside(m,b,s)&&fabs((*f)(s))<fabs(fb)){
                if(fb*(*f)(s)<0){
                    a=b;fa=fb;
                    b=s;fb=(*f)(s);
                }
                else{
                    b=s;fb=(*f)(s);
                    if(fb*(*f)(m)<0){
                        a=m;fa=(*f)(m);
                    }
                }
            }
            else{
                pos_1:
                if((*f)(m)*fb<0){
                    a=m;
                    fa=(*f)(m);
                }
                else{
                    b=m;
                    fb=(*f)(m);
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
//右端点瑕点的瑕积分
double integratePeriod_right_infty(double y1,double y2,double y){
    double T=0,T0=-1;
    double yi;
    if(!isInside(y1,y2,y))return -404404;
    double y_sum=(y1+y2)/2,y_minus=(y2-y1)/2;
    double h=acos((y-y_sum)/y_minus);
    int n,i;
    for ( n = 1; fabs(T-T0)>EPS_ITG; n++)
    {
        double S=0;
        for ( i = 0; i < n; i++)
        {
            yi=(y1+y2)/2+(y2-y1)/2*cos((2*i+1)*h/(2*n));
            S+=sqrt((y2-yi)*(yi-y1))/py(yi);
        }
        T0=T;
        T=S*h/n;
    }
    return T;
}

double getW(void){
    double H=0.5*PY_0*PY_0+0.5*pow(sqrt(xt*xt+Y_0*Y_0)-1,2)+gt*Y_0;
    
    if(xt==0){
        Ymax=0;
        Ymin1=DekkerRoot(0,3,eq1);
        Ymin2=DekkerRoot(-3,0,eq1);
    }
    else{
        Y_threshold=DekkerRoot(0,2,f);
        Ymin1=DekkerRoot(Y_threshold,3,eq1);
        Ymin2=DekkerRoot(-3,-Y_threshold,eq1);
        if(gt>0)Ymax=DekkerRoot(0,Y_threshold,eq1);
        else Ymax=DekkerRoot(-Y_threshold,0,eq1);
    }
    if(Ymin1<-100||Ymin2<-100||Ymax<-100||Y_threshold<-100)return -404404;
    
    double y1,y2;//相空间积分上下限
    if(H>V(Ymax)){
            y1=DekkerRoot(-MAX,Ymin2,pySquare);
            y2=DekkerRoot(Ymin1,MAX,pySquare);
        }
    else{
        if(Y_0>Ymax){
            y1=DekkerRoot(Ymax,Ymin1,pySquare);
            y2=DekkerRoot(Ymin1,MAX,pySquare);
        }
        else{
            y1=DekkerRoot(-MAX,Ymin2,pySquare);
            y2=DekkerRoot(Ymin2,Ymax,pySquare);
        }
    }

    printf("%.15lf %.15lf\n",pySquare(y1),pySquare(y2));
    double T;//周期
    T=2*integratePeriod(y1,y2);
    if(T==0)return INT_MAX;
    double w=2*PI/T;
    return w;
}

int main(){
    FILE *fdata=fopen("C://C_files//HW2//q5-berry-N=256.txt","w");
    if(fdata==NULL){
        return 0;
    }
    double y_o,y_n,py_o,py_n;
    double I1,I2=0;
    double wi,yi,pyi,yf,pyf;
    yi=Y_0;pyi=PY_0;
    double w;
    y_n=Y_0;
    py_n=PY_0+dt/2*((Y_0+dt/4*PY_0)*(-1+1/sqrt(x(dt/4)*x(dt/4)+pow(Y_0+dt/4*PY_0,2)))-g(dt/4));
    t=0;xt=x(t);gt=g(t);
    w=getW();
    wi=w;
    fprintf(fdata,"%.15lf  %.15lf %lf %lf\n",Y_0,PY_0,t,w);
    
    
    for(t=dt;t<N;t+=dt){//蛙跳法
        xt=x(t);
        gt=g(t);
        y_o=y_n;
        py_o=py_n;
        y_n=y_o+dt*py_o;
        py_n=py_o+dt*(y_n*(-1+1/sqrt(xt*xt+y_n*y_n))-gt);
        Y_0=y_n;
        PY_0=(py_n+py_o)/2;
        w=getW();
        I2+=w*dt;
        fprintf(fdata,"%.15lf  %.15lf %lf %lf\n",Y_0,PY_0,t,w);
    }
    yf=Y_0;pyf=PY_0;
    
    yf=Y_0;pyf=PY_0;
    
    t=0;xt=x(t);gt=g(t);
    double H=0.5*pyi*pyi+0.5*pow(sqrt(xt*xt+yi*yi)-1,2)+gt*yi;

    if(xt==0){
        Ymax=0;
        Ymin1=DekkerRoot(0,3,eq1);
        Ymin2=DekkerRoot(-3,0,eq1);
    }
    else{
        Y_threshold=DekkerRoot(0,2,f);
        Ymin1=DekkerRoot(Y_threshold,3,eq1);
        Ymin2=DekkerRoot(-3,-Y_threshold,eq1);
        if(gt>0)Ymax=DekkerRoot(0,Y_threshold,eq1);
        else Ymax=DekkerRoot(-Y_threshold,0,eq1);
    }
    if(Ymin1<-100||Ymin2<-100||Ymax<-100||Y_threshold<-100)return -404404;
    
    double y1,y2;//相空间积分上下限
    if(H>V(Ymax)){
            y1=DekkerRoot(-MAX,Ymin2,pySquare);
            y2=DekkerRoot(Ymin1,MAX,pySquare);
        }
    else{
        if(Y_0>Ymax){
            y1=DekkerRoot(Ymax,Ymin1,pySquare);
            y2=DekkerRoot(Ymin1,MAX,pySquare);
        }
        else{
            y1=DekkerRoot(-MAX,Ymin2,pySquare);
            y2=DekkerRoot(Ymin2,Ymax,pySquare);
        }
    }

    I1=wi*(integratePeriod(y1,y2)+integratePeriod_right_infty(y1,y2,yf));
    fprintf(fdata,"%lf %lf %lf\n",I1,I2,(I1-I2)-round((I1-I2)/(2*PI))*2*PI);

    fclose(fdata);
    return 0;
}