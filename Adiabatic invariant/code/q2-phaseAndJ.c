#include <stdio.h>
#include <math.h>
#define EPS_ITG 1e-8//积分特征精度
#define EPS_ROOT 1e-14//求根特征精度
#define PI 3.1415926535897932
double N=4;//v=1/N
double dt=1e-4;//迭代步长
double Y_0=0.1;//储存每个时刻的初值位置
double PY_0=0;//储存每个时刻的初值动量
double Y_o,Y_n,PY_o,PY_n;//迭代参数
double xt;//xt=x(t)
double x(double t);//给定t，返回x值
double pySquare(double y);//计算给定位置动量的平方
int isInside(double a,double b,double x);
double DekkerPySquare(double a, double b);//Dekker算法求pySquare函数零点
double py(double y);//计算给定位置的动量
double integrate_phaseSpace(double y1,double y2);//相空间定上下限积分
double getJ();//计算J

int main(){
    FILE *fdata=fopen("C://C_files//HW2//q2-1.txt","w");
    if(fdata==NULL){
        return 0;
    }
    Y_n=Y_0;
    PY_n=PY_0+dt/2*(Y_0+dt/4*PY_0)*(-1+1/sqrt(pow(x(1/4*dt),2)+pow(Y_0+dt/4*PY_0,2)));
    xt=x(0);
    double J;
    J=getJ();
    fprintf(fdata,"%lf  %lf %lf %lf\n",Y_0,PY_0,x(0),J);
    printf("%lf\n",x(0));//控制台查看运行进度
    double t=dt;
    xt=x(t);
    while(!(xt<0)){//蛙跳法
        Y_o=Y_n;
        PY_o=PY_n;
        Y_n=Y_o+dt*PY_o;
        PY_n=PY_o+dt*Y_n*(-1+1/sqrt(xt*xt+Y_n*Y_n));
        Y_0=Y_n;
        PY_0=(PY_n+PY_o)/2;
        J=getJ();
        fprintf(fdata,"%lf  %lf %lf %lf\n",Y_0,PY_0,xt,J);
        printf("%lf\n",xt);//控制台查看运行进度，无特别含义

        t+=dt;
        xt=x(t);
    }
    fclose(fdata);
    return 0;
}

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

    while(fabs(fb)>EPS_ROOT&&fabs(a-b)>EPS_ROOT){
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
    //由于求积分上下限精度有限，有可能在积分过程中遇到一个y值
    //动量平方是个很小的负数，这时不妨将它抹去，于是有下面置零操作
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
    while (fabs(T-T0)>EPS_ITG){//二分法求定积分
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
    double yc=DekkerPySquare(y0,2.1);//2.1是估计的上界
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

