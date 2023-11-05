#include <stdio.h>
#include <math.h>
#define EPS 10e-15
#define ITV 1e-4
double x;//参数x从0变到1

double f(double y){
    double f=pow(y*y+x*x,3)-pow(x,4);
    return f;
}
int isInside(double a,double b,double x){//判断x是否在a,b之间
    if((a<x&&x<b)||(b<x&&x<a))return 1;
    else return 0;
}
//Dekker method to find root of f(x) in given range [x1,x2]
double DekkerRoot (double x1,double x2){
    double a=x1,b=x2;
    double temp,m,s,fa,fb;
    if(f(a)*f(b)>0){
        return -404404;
    }
    if(fabs(f(b))>fabs(f(a))){
        temp=a;
        a=b;
        b=temp;
    }
    fa=f(a);
    fb=f(b);

    while(fb!=0&&fabs(a-b)>EPS){
        m=(a+b)/2;
        if(fa!=fb){
            s=(a*fb-b*fa)/(fb-fa);
            if(isInside(m,b,s)&&fabs(f(s))<fabs(fb)){
                if(fb*f(s)<0){
                    a=b;fa=fb;
                    b=s;fb=f(s);
                }
                else{
                    b=s;fb=f(s);
                    if(fb*f(m)<0){
                        a=m;fa=f(m);
                    }
                }
            }
            else{
                pos_1:
                if(f(m)*fb<0){
                    a=m;
                    fa=f(m);
                }
                else{
                    b=m;
                    fb=f(m);
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

int main(){
    FILE *fdata=fopen("C://C_files//HW2//plot_x-g_data.txt","w");
    if(fdata==NULL){
        return 0;
    }
    double y,g;
    for ( x = ITV; x < 1; x+=ITV)
    {
        y=DekkerRoot(0,1);//计算切点
        if(y==-404404){
            fprintf(fdata,"x=%lf, Wrong section!\n",x);
            continue;
        }
        g=y/(sqrt(x*x+y*y))-y;//将切点回代，得到临界g
        fprintf(fdata,"%lf %.15lf %lf\n",x,y,g);
        printf("%lf\n",x);
    }
    fclose(fdata);
    return 0;
}