#include "matrix.h"
#define H 0.05
#define OMEGA 1.9
int** ijtol_store;
int *ltoi,*ltoj;
int n1max,n2max,n2mid,n1diag,n2diag;
int flag;

int isinside(int i,int j);
int add_position(int i,int j,int l);
int ijtol(int i, int j);
int set_flag();
double AmX(int m,vector* x);
double b(int l);
int GaussSeidel(vector* x);

int main(){
    FILE* fdata=fopen("C:/C_files/HW5/q4_inverse.txt","w");
    int count=0;
    n1max=n2max=n2mid=n1diag=n2diag=0;
    int i=1,j=1;
    while(isinside(i,j)){
        count++;n2max++;
        if(isinside(0,j))n2mid++;
        if(j==1)n1max++;
        if(j>n2diag&&abs(i-j)<=abs(n1diag-n2diag)){n1diag=i;n2diag=j;}
        i++;
        while(isinside(i,j)){
            count++;
            if(j==1)n1max++;
            if(j>n2diag&&abs(i-j)<=abs(n1diag-n2diag)){n1diag=i;n2diag=j;}
            i++;
        }
        j++;i=1;
    }
    while(1){
        while(i<=n1max&&!isinside(i,j)){
            i++;
        }
        if(i>n1max)break;
        count++;n2max++;
        if(j>n2diag&&abs(i-j)<=abs(n1diag-n2diag)){n1diag=i;n2diag=j;}
        i++;
        while(isinside(i,j)){
            count++;
            if(j>n2diag&&abs(i-j)<=abs(n1diag-n2diag)){n1diag=i;n2diag=j;}
            i++;
        }
        j++;i=1;
    }
    printf("%d  %d  %d  %d  %d\n",n1max,n2max,n2mid,n1diag,n2diag);
    int total=4*count+2*(n2mid+n1max)+1;
    printf("%d\n",total);

    ijtol_store=createNew2DIntArray(2*n1max+3,2*n2max+3);
    ltoi=createNew1DIntArray(total+5);
    ltoj=createNew1DIntArray(total+5);
    int p,q,l=0;
    for(q=0;q<2*n2max+3;q++){
        for(p=0;p<2*n1max+3;p++){
            ijtol_store[p][q]=-1;
        }
    }
    i=j=0; 
    add_position(i,j,l++);
    for(i=1;i<=n1max;i++){
        add_position(i,j,l++);
        add_position(-i,j,l++);
    }
    i=0;
    for(j=1;j<=n2mid;j++){
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
    while(1){
        while(i<=n1max&&!isinside(i,j)){
            i++;
        }
        if(i>n1max)break;
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
    printf("%d\n",l);
    if(l!=total){
        printf("storage error: l should equal to total\n");
        return 0;
    }
    vector* x=createNewVector(total);
    set_flag();
    if(GaussSeidel(x)){
        for(q=2*n2max+2;q>=0;q--){
            for(p=0;p<2*n1max+3;p++){
                if(ijtol_store[p][q]<0)fprintf(fdata,"%8d+I\t",ijtol_store[p][q]);
                else{
                    fprintf(fdata,"%10lf\t",x->el[ijtol_store[p][q]]);
                }
            }
            fprintf(fdata,"\n");
        }

        // for(l=0;l<total;l++){
        //     fprintf(fdata,"%lf  %lf %lf\n",ltoi[l]*H,ltoj[l]*H,x->el[l]);
        // }
    }
    
    printf("%lf %lf\n",x->el[ijtol(0,-n2mid)],x->el[ijtol(-n1diag,n2diag)]);
    return 0;
}

int isinside(int i,int j){
    double x,y;
    x=abs(i)*H;y=abs(j)*H;
    if(x*x+y*y-(8*x*x)/(x*x+y*y)-(16*pow(x,4))/pow(x*x+y*y,2)>1)return 0;
    return 1;
}

int add_position(int i,int j,int l){
    if(i<-(n1max+1)||i>(n1max+1)||j<-(n2max+1)||j>(n2max+1))return 0;
    int m,n;
    m=i+(n1max+1);
    n=(n2max+1)-j;
    ijtol_store[m][n]=l;
    ltoi[l]=i;
    ltoj[l]=j;
    return 1;
}

int ijtol(int i, int j){
    if(i<-(n1max+1)||i>(n1max+1)||j<-(n2max+1)||j>(n2max+1))return -1;
    return ijtol_store[i+(n1max+1)][(n2max+1)-j];
}
int set_flag(){
    if(isinside(n1diag-1,n2diag)&&isinside(n1diag+1,n2diag))flag=1;
    else if(!isinside(n1diag-1,n2diag)){
        flag=2;
    }
    else{
        flag=3;
    }
}
double AmX(int m,vector* x){
    int i,j;
    i=ltoi[m];j=ltoj[m];
    if(i==(n1max)&&j==0){
        return -4*x->el[m]+4*x->el[ijtol(i-1,j)];
    }
    else if(i==(n1diag)&&j==-n2diag){
        switch (flag)
        {
        case 1:
            return -4*x->el[m]+4*x->el[ijtol(i,j+1)];
            break;

        case 2:
            return -4*x->el[m]+4*x->el[ijtol(i+1,j+1)];
            break;

        case 3:
            return -4*x->el[m]+4*x->el[ijtol(i-1,j+1)];
            break;
        }
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
    if(i==(n1max)&&j==0)return -4;
    if(i==(n1diag)&&j==-n2diag)return 4;
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
    printf("OMEGA=%lf,N=%d\n",OMEGA,count);
    return 1;
}
