#include <stdio.h>
#include <math.h>
double N=256;
double dt=1e-3;
double y_0=0.1;
double py_0=0;

double x(double t){
    if(t<0)return -404404;
    return 2-t/N;
}
int main(){
    FILE *fdata=fopen("C://C_files//HW2//q2-phase-4.txt","w");
    if(fdata==NULL){
        return 0;
    }
    double y_o,y_n,py_o,py_n;
    y_n=y_0;
    py_n=py_0+dt/2*(y_0+dt/4*py_0)*(-1+1/sqrt(pow(x(1/4*dt),2)+pow(y_0+dt/4*py_0,2)));
    fprintf(fdata,"%lf  %lf\n",y_0,py_0);
    double t=dt;
    double xt=x(t);
    while(!(xt<0)){
        y_o=y_n;
        py_o=py_n;
        y_n=y_o+dt*py_o;
        py_n=py_o+dt*y_n*(-1+1/sqrt(xt*xt+y_n*y_n));
        fprintf(fdata,"%lf  %lf\n",y_n,(py_n+py_o)/2);

        t+=dt;
        xt=x(t);
    }
    fclose(fdata);
    return 0;
}