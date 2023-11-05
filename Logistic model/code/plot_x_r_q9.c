//本程序需要手动添加周期分点和周期数
#include <stdio.h>
#include <math.h>
#define EPS 1e-10//收敛精度

double r;
double x;

double f(double x){
    double y=r*sin(x);
    return y;
}
double converge(double x,int N){//迭代直到数列收敛,N为序列周期数
    double y=x;
    double y0=-9;
    int i;
    while(fabs(y-y0)>EPS){//收敛条件
        y0=y;
        for(i=0;i<N;i++){//每次迭代N次
            y=f(y);
        }
    }
    return y;
}
int main(){//数据将输入到指定路径下的txt文件中
    FILE *fdata=fopen("C://C_files//HW1//x_r_data_q9.txt","w");
    if (fdata==NULL)
    {
        return 0;
    }
    int j;
    r = 0.0001;
    while (r < 2.71925)
    {
        x=0.9;
        int i=0,j;
        if(r<2.2618){
            x=converge(x,1);//收敛后输出1个数
           fprintf(fdata,"%lf   %lf\n",r,x);
        }
        else if(r<2.6178){
            x=converge(x,2);//收敛后输出2个数
            for ( j = 0; j < 2; j++)
            {
               fprintf(fdata,"%lf   %lf\n",r,x);
                x=f(x);
            }
        }
        else if(r<2.6974){
            x=converge(x,4);//收敛后输出4个数
            for ( j = 0; j < 4; j++)
            {
               fprintf(fdata,"%lf   %lf\n",r,x);
                x=f(x);
            }
        }
        else if(r<2.7146){
            x=converge(x,8);//收敛后输出8个数
            for ( j = 0; j < 8; j++)
            {
               fprintf(fdata,"%lf   %lf\n",r,x);
                x=f(x);
            }
        }
        else if(r<2.7183){
            x=converge(x,16);//收敛后输出16个数
            for ( j = 0; j < 16; j++)
            {
               fprintf(fdata,"%lf   %lf\n",r,x);
                x=f(x);
            }
        }
        else if(r<2.7191){
            x=converge(x,32);//收敛后输出32个数
            for ( j = 0; j < 32; j++)
            {
               fprintf(fdata,"%lf   %lf\n",r,x);
                x=f(x);
            }
        }
        else if(r<2.7193){
            x=converge(x,64);//收敛后输出64个数
            for ( j = 0; j < 64; j++)
            {
               fprintf(fdata,"%lf   %lf\n",r,x);
                x=f(x);
            }
        }
        r+=1e-4;
    }    
    fclose(fdata);
    return 0;
}

