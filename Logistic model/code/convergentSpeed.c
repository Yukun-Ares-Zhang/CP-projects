#include <stdio.h>
#include <math.h>
#define EPS1 1e-10//原序列收敛特征长度
#define EPS2 1e-7//收敛速度收敛特征长度

double r=0.0001;
double f(double x){
    double y=r*x*(1-x);
    return y;
}
double Nf(double x, int N){//计算序列迭代N次后的结果
    int i;
    double y=x;
    for ( i = 0; i < N; i++)
    {
        y=f(y);
    }
    return y;
}
 int main(){//数据将输入到指定路径下的txt文件中
    FILE *fdata=fopen("C://C_files//HW1//convergentSpeed_Data.txt","w");
    if (fdata==NULL)
    {
        return 0;
    }
    int i;
    while (r < 3.569892)
    {
        double x00=0.9;//$x_n$
        double x0;//$x_{n+N}$
        double x;//$x_{n+2N}$
        double V=-9,V0=-99,V00=-999,V000=-9999,V0000=-99999;
        //分别储存临近五次计算的收敛速度V，只有在相邻收敛速度之差都小于EPS2时，迭代才会收敛
        //|x-x0|<EPS1时被除数过小，迭代也会终止
        if(r<3){
           x0=Nf(x00,1);
           x=Nf(x0,1);
            while(((fabs(V-V0)>EPS2)||(fabs(V0-V00)>EPS2)||(fabs(V00-V000)>EPS2)||
            (fabs(V000-V0000)>EPS2))&&fabs(x-x0)>EPS1){
            V0000=V000;
            V000=V00;
            V00=V0;
            V0=V;
            V=fabs((x0-x00)/(x-x0));
            x00=x0;
            x0=x;
            x=Nf(x,1);
            }
            fprintf(fdata,"%lf %lf\n",r,V); 
        }
        else if(r<3.449490){
            x0=Nf(x00,2);
           x=Nf(x0,2);
            while(((fabs(V-V0)>EPS2)||(fabs(V0-V00)>EPS2)||(fabs(V00-V000)>EPS2)||
            (fabs(V000-V0000)>EPS2))&&fabs(x-x0)>EPS1){
            V0000=V000;
            V000=V00;
            V00=V0;
            V0=V;
            V=fabs((x0-x00)/(x-x0));
            x00=x0;
            x0=x;
            x=Nf(x,2);
            }
            fprintf(fdata,"%lf %lf\n",r,V); 
        }
        else if(r<3.544091){
            x0=Nf(x00,4);
           x=Nf(x0,4);
            while(((fabs(V-V0)>EPS2)||(fabs(V0-V00)>EPS2)||(fabs(V00-V000)>EPS2)||
            (fabs(V000-V0000)>EPS2))&&fabs(x-x0)>EPS1){
            V0000=V000;
            V000=V00;
            V00=V0;
            V0=V;
            V=fabs((x0-x00)/(x-x0));
            x00=x0;
            x0=x;
            x=Nf(x,4);
            }
            fprintf(fdata,"%lf %lf\n",r,V); 
        }
        else if(r<3.564408){
            x0=Nf(x00,8);
           x=Nf(x0,8);
            while(((fabs(V-V0)>EPS2)||(fabs(V0-V00)>EPS2)||(fabs(V00-V000)>EPS2)||
            (fabs(V000-V0000)>EPS2))&&fabs(x-x0)>EPS1){
            V0000=V000;
            V000=V00;
            V00=V0;
            V0=V;
            V=fabs((x0-x00)/(x-x0));
            x00=x0;
            x0=x;
            x=Nf(x,8);
            }
            fprintf(fdata,"%lf %lf\n",r,V); 
        }
        else if(r<3.568760){
            x0=Nf(x00,16);
           x=Nf(x0,16);
            while(((fabs(V-V0)>EPS2)||(fabs(V0-V00)>EPS2)||(fabs(V00-V000)>EPS2)||
            (fabs(V000-V0000)>EPS2))&&fabs(x-x0)>EPS1){
            V0000=V000;
            V000=V00;
            V00=V0;
            V0=V;
            V=fabs((x0-x00)/(x-x0));
            x00=x0;
            x0=x;
            x=Nf(x,16);
            }
            fprintf(fdata,"%lf %lf\n",r,V); 
        }
        else if(r<3.569692){
            x0=Nf(x00,32);
           x=Nf(x0,32);
            while(((fabs(V-V0)>EPS2)||(fabs(V0-V00)>EPS2)||(fabs(V00-V000)>EPS2)||
            (fabs(V000-V0000)>EPS2))&&fabs(x-x0)>EPS1){
            V0000=V000;
            V000=V00;
            V00=V0;
            V0=V;
            V=fabs((x0-x00)/(x-x0));
            x00=x0;
            x0=x;
            x=Nf(x,32);
            }
            fprintf(fdata,"%lf %lf\n",r,V); 
        }
        else if(r<3.569892){
            x0=Nf(x00,64);
           x=Nf(x0,64);
            while(((fabs(V-V0)>EPS2)||(fabs(V0-V00)>EPS2)||(fabs(V00-V000)>EPS2)||
            (fabs(V000-V0000)>EPS2))&&fabs(x-x0)>EPS1){
            V0000=V000;
            V000=V00;
            V00=V0;
            V0=V;
            V=fabs((x0-x00)/(x-x0));
            x00=x0;
            x0=x;
            x=Nf(x,64);
            }
            fprintf(fdata,"%lf %lf\n",r,V); 
        }
        
        if(r>3.15&&r<3.3)r+=1e-6;//部分区间取点加密
        else if(r>3.47&&x<3.52)r+=1e-6;
        else if(r>3.544)r+=1e-6;
        else r+=1e-5;
    }
    return 0;
 }
 
 /* int main(){
    FILE *fdata=fopen("C://C_files//HW1//data.txt","w");
    if (fdata==NULL)
    {
        return 0;
    }
    double x=0.9;
    int i;
    for ( r = 0.0001; r < 3; r+=0.0001)
    {
        x=0.9;
        double x0=-99,x_EXACT;
        double V=-9,V0=-99,V00=-999,V000=-9999,V0000=-99999;
        while(fabs(x-x0)>EPS1){
            x0=x;
            x=f(x);
        }
        x_EXACT=x;

        x0=0.9;
        x=f(x0);
        while(((fabs(V-V0)>EPS2)||(fabs(V0-V00)>EPS2)||(fabs(V00-V000)>EPS2)||(fabs(V000-V0000)>EPS2))&&fabs(x-x_EXACT)>EPS1){
            V0000=V000;
            V000=V00;
            V00=V0;
            V0=V;
            V=fabs((x0-x_EXACT)/(x-x_EXACT));
            x0=x;
            x=f(x);
        }
        fprintf(fdata,"%lf %lf\n",r,V); 
    }
    return 0;

    
 } */