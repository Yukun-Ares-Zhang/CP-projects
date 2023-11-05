#include <stdio.h>
#include <math.h>
#define MAX 8
#define EPS_ITG 1e-8
#define EPS_ROOT 1e-14
#define PI 3.1415926535897932
#define X 0.2
#define Y_THRESHOLD 0.277417932961286//与X取值有关
#define G_THRESHOLD 0.53375//与X取值有关
double N=256;
double dt=1e-2;
double Y_0=0.52;
double PY_0=0;
double t;
double Ymin1,Ymax,Ymin2;
double Ymin;

double g(double t){
    double g=0.3*cos(2*PI*t/N);
    return g;
}

double x(double t){
    double x=0.3*sin(2*PI*t/N);
}
int main(){
    FILE *fdata=fopen("C://C_files//HW2//q45-phase.txt","w");
    if(fdata==NULL){
        return 0;
    }
    double y_o,y_n,py_o,py_n;
    y_n=Y_0;
    py_n=PY_0+dt/2*((Y_0+dt/4*PY_0)*(-1+1/sqrt(x(dt/4)*x(dt/4)+pow(Y_0+dt/4*PY_0,2)))-g(dt/4));
    t=0;
    fprintf(fdata,"%.15lf  %.15lf\n",Y_0,PY_0);
    printf("%lf\n",g(t));
    
    
    for(t=dt;t<=N/2;t+=dt){
        y_o=y_n;
        py_o=py_n;
        y_n=y_o+dt*py_o;
        py_n=py_o+dt*(y_n*(-1+1/sqrt(x(t)*x(t)+y_n*y_n))-g(t));
        Y_0=y_n;
        PY_0=(py_n+py_o)/2;
        fprintf(fdata,"%.15lf  %.15lf\n",Y_0,PY_0);
        printf("%lf\n",g(t));
    }
    fclose(fdata);
    return 0;
}