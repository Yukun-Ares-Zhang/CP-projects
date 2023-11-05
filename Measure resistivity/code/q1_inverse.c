#include "matrix.h"
#define NX 300
#define NY 200
#define OMEGA 1.5
#define H 0.01

int xytol(int x,int y);
int ltox(int l);
int ltoy(int l);
double AmX(int m,vector* x);
double b(int l);
int GaussSeidel(vector* x);

int main(){
    FILE* fdata=fopen("C:/C_files/HW5/q1_NX=300_inverse.txt","w");
    vector* x=createNewVector((NX+1)*(NY+1));
    if(GaussSeidel(x)){
        int i,j;
        for(j=0;j<=NY;j++){
            for(i=0;i<=NX;i++){
                fprintf(fdata,"%10lf\t",x->el[xytol(i,j)]);
            }
            fprintf(fdata,"\n");
        }
    }
    printf("%lf %lf\n",x->el[xytol(NX,0)],x->el[xytol(NX,NY)]);
    return 0;
}

int xytol(int x,int y){
    if(x<0||x>NX||y<0||y>NY)return -1;
    return x+y*(NX+1);
}
int ltox(int l){
    if(l<0||l>(NX+1)*(NY+1))return -1;
    return l%(NX+1);
}
int ltoy(int l){
    if(l<0||l>(NX+1)*(NY+1))return -1;
    return (int)l/(NX+1);
}
double AmX(int m,vector* x){
    int i,j;
    i=ltox(m);j=ltoy(m);
    if(i==0&&j==0){
        return -4*x->el[m]+4*x->el[xytol(i+1,j+1)];
    }
    else if(i==0&&j==NY){
        return -4*x->el[m]+4*x->el[xytol(i+1,j-1)];
    }
    else{
        double s=0;
        int counted1,counted2,counted3,counted4;
        counted1=counted2=counted3=counted4=0;
        s+=-4*x->el[m];
        if(i==0){
            s+=2*x->el[xytol(1,j)];counted1=1;
        }
        if(i==NX){
            s+=2*x->el[xytol(NX-1,j)];counted3=1;
        }
        if(j==0){
            s+=2*x->el[xytol(i,1)];counted2=1;
        }
        if(j==NY){
            s+=2*x->el[xytol(i,NY-1)];counted4=1;
        }
        if(i!=NX&&!counted1){
            s+=x->el[xytol(i+1,j)];
        }
        if(j!=NY&&!counted2){
            s+=x->el[xytol(i,j+1)];
        }
        if(i!=0&&!counted3){
            s+=x->el[xytol(i-1,j)];
        }
        if(j!=0&&!counted4){
            s+=x->el[xytol(i,j-1)];
        }
        return s;
    }
}
double b(int l){
    int i,j;
    i=ltox(l);j=ltoy(l);
    if(i==0&&j==0)return -4;
    if(i==0&&j==NY)return 4;
    return 0;
}
int GaussSeidel(vector* x){
    if(x->n<1)return 0;
    double dx=0;
    double dx_max=99;
    int count=0;
    while(dx_max>10*EPS){
        dx=0;
        dx_max=dx;
        count++;
        int i,j;
        for(i=0;i<x->n;i++){
            dx=b(i);
            dx-=AmX(i,x);
            dx*=OMEGA;
            dx/=(-4.);
            x->el[i]+=dx;
            if(fabs(dx)>dx_max)dx_max=fabs(dx);
        }
        //printf("%d  %lf\n",count,dx_max);
    }
    printf("%d\n",count);
    return 1;
}
