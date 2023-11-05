#include <stdio.h>
#include <math.h>
#define MAX 2000000
#define EPS 1e-6
#define PI 3.1415926535897932
double N=4;
double dt1=1e-3;
double dt2=1e-3;
double y_0=0.1;
double py_0=0;

double x(double t){
    if(t<0)return -404404;
    return 2-t/N;
}
int main(){
    FILE *fdata=fopen("C://C_files//HW2//q2-J-1-try.txt","w");
    if(fdata==NULL){
        return 0;
    }
    double t=0;
    double xt=x(t);
    while (!(xt<0))
    {
        double y_oo=-1,y_o=0,y_n,py_o,py_n;
        y_n=y_0;
        py_n=py_0+dt2/2*(y_0+dt2/4*py_0)*(-1+1/sqrt(pow(xt,2)+pow(y_0+dt2/4*py_0,2)));
        int i;
        double J=0;
        while(fabs(y_o-y_0)>EPS||fabs(y_oo-y_0)>EPS||fabs(y_n-y_0)>EPS||(y_o-y_oo)*(y_n-y_o)>0){
            y_oo=y_o;
            y_o=y_n;
            py_o=py_n;
            y_n=y_o+dt2*py_o;
            py_n=py_o+dt2*y_n*(-1+1/sqrt(xt*xt+y_n*y_n));
            J+=py_o*(y_n-y_o);
        }
        J-=py_o*(y_n-y_o);
        fprintf(fdata,"%lf  %lf\n",xt,J/(2*PI));

        t+=dt1;
        xt=x(t);
    }
    
    fclose(fdata);
    return 0;
}

//for(i=0;i<MAX;i++)