#include "matrix.h"
#define H 0.01
#define LY 1.
#define LX 1.5
#define OMEGA 1.5
int** ijtol_store;
int *ltoi,*ltoj;
int nx1,n2,total,ny_mid,nx_mid;

int isinside(int i,int j);
int add_position(int i,int j,int l);
int ijtol(int i, int j);
double AmX(int m,vector* x);
double b(int l);
int GaussSeidel(vector* x);

int main(){
    FILE* fdata=fopen("C:/C_files/HW5/q2_H=0.01.txt","w");
    int count=0,n1=0;
    n2=0,nx1=0;
    ny_mid=n2,nx_mid=n1;
    int renew_nx_mid=0;
    int i=1,j=1;
    while(isinside(i,j)){
        count++;n1++;n2++;i++;
        if(fabs(LY/2-n2*H)<fabs(LY/2-ny_mid*H)){
            ny_mid=n2;renew_nx_mid=1;
        }
        while(isinside(i,j)){
            count++;n1++;i++;
        }
        if(j==1)nx1=n1;
        if(renew_nx_mid){
            nx_mid=n1;
            renew_nx_mid=0;
        }
        j++;i=1;n1=0;
    }
    printf("%d  %d  %d  %d   %d\n",nx1,n2,count,nx_mid,ny_mid);
    int total=4*count+2*(n2+nx1)+1;
    printf("%d\n",total);

    ijtol_store=createNew2DIntArray(2*nx1+3,2*n2+3);
    ltoi=createNew1DIntArray(total+5);
    ltoj=createNew1DIntArray(total+5);
    int p,q,l=0;
    for(q=0;q<2*n2+3;q++){
        for(p=0;p<2*nx1+3;p++){
            ijtol_store[p][q]=-1;
        }
    }
    i=j=0; 
    add_position(i,j,l++);
    for(i=1;i<=nx1;i++){
        add_position(i,j,l++);
        add_position(-i,j,l++);
    }
    i=0;
    for(j=1;j<=n2;j++){
        add_position(i,j,l++);
        add_position(i,-j,l++);
    }
    i=1,j=1;
    while(isinside(i,j)){
        while(isinside(i,j)){
            add_position(i,j,l++);
            add_position(-i,j,l++);
            add_position(-i,-j,l++);
            add_position(i,-j,l++);
            i++;
        }
        j++;i=1;
    }
    printf("position stored\n");
    if(l!=total){
        printf("storage error: l should equal to total\n");
        return 0;
    }
    vector* x=createNewVector(total);
    if(GaussSeidel(x)){
        for(q=2*n2+2;q>=0;q--){
            for(p=0;p<2*nx1+3;p++){
                if(ijtol_store[p][q]<0)fprintf(fdata,"%8d+I\t",ijtol_store[p][q]);
                else{
                    fprintf(fdata,"%10lf\t",x->el[ijtol_store[p][q]]);
                }
            }
            fprintf(fdata,"\n");
        }
    }
    printf("%lf %lf\n",x->el[ijtol(nx_mid,ny_mid)],x->el[ijtol(-nx_mid,ny_mid)]);
    return 0;
}

int isinside(int i,int j){
    double x,y;
    x=abs(i)*H;y=abs(j)*H;
    if(y+LY/LX*x>LY)return 0;
    return 1;
}

int add_position(int i,int j,int l){
    if(i<-(nx1+1)||i>(nx1+1)||j<-(n2+1)||j>(n2+1))return 0;
    int m,n;
    m=i+(nx1+1);
    n=(n2+1)-j;
    ijtol_store[m][n]=l;
    ltoi[l]=i;
    ltoj[l]=j;
    return 1;
}
int ijtol(int i, int j){
    if(i<-(nx1+1)||i>(nx1+1)||j<-(n2+1)||j>(n2+1))return -1;
    return ijtol_store[i+(nx1+1)][(n2+1)-j];
}
double AmX(int m,vector* x){
    int i,j;
    i=ltoi[m];j=ltoj[m];
    if(i==(-nx_mid)&&j==-ny_mid){
        return -4*x->el[m]+4*x->el[ijtol(i+1,j+1)];
    }
    else if(i==(nx_mid)&&j==-ny_mid){
        return -4*x->el[m]+4*x->el[ijtol(i-1,j+1)];
    }
    else{
        double s=0;
        int counted1,counted2,counted3,counted4;
        counted1=counted2=counted3=counted4=0;
        s+=-4*x->el[m];
        if(!isinside(i-1,j)){
            s+=2*x->el[ijtol(i+1,j)];counted1=1;
        }
        if(!isinside(i+1,j)){
            s+=2*x->el[ijtol(i-1,j)];counted3=1;
        }
        if(!isinside(i,j-1)){
            s+=2*x->el[ijtol(i,j+1)];counted2=1;
        }
        if(!isinside(i,j+1)){
            s+=2*x->el[ijtol(i,j-1)];counted4=1;
        }
        if(isinside(i+1,j)&&!counted1){
            s+=x->el[ijtol(i+1,j)];
        }
        if(isinside(i,j+1)&&!counted2){
            s+=x->el[ijtol(i,j+1)];
        }
        if(isinside(i-1,j)&&!counted3){
            s+=x->el[ijtol(i-1,j)];
        }
        if(isinside(i,j-1)&&!counted4){
            s+=x->el[ijtol(i,j-1)];
        }
        return s;
    }
}

double b(int l){
    int i,j;
    i=ltoi[l];j=ltoj[l];
    if(i==(-nx_mid)&&j==-ny_mid)return -4;
    if(i==(nx_mid)&&j==-ny_mid)return 4;
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