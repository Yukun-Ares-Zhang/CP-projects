//本程序计算给定G值基态能量和波函数
//本程序自动判断真收敛和伪收敛
//真收敛给出收敛本征值本征矢，伪收敛给出两个横跳本征值本征矢
//需要G<0，G值需要在G_CHAR中设定
#include "matrix.h"
#include <string.h>
#define NUM 100
#define LOW -20
#define UP 20
#define G_CHAR "-2.0"//手动设置G值
#define PDIS 0.
#define EPS_J 1e-7
double G;
double h;
double trailFunc(double x){
    return 1/pow(PI,0.25)*pow(E,-x*x/2);
}

int normalize(vector* pv){
    int i;
    double s=0;
    for(i=0;i<pv->n;i++){
        s+=pow(pv->el[i],2)*h;
    }
    if(s<EPS){
        printf("Normalization error\n");
        return 0;
    }
    double sq=sqrt(s);
    for(i=0;i<pv->n;i++){
        pv->el[i]/=sq;
    }
    return 1;
}
int resetCoeff(vector *phi,vector *pa,vector *pb,vector *pc){
    if(phi==NULL||pa==NULL||pb==NULL||pc==NULL)return 0;
    if(phi->n!=pa->n||pa->n!=pb->n||pb->n!=pc->n){
        printf("Coefficient reset error\n");
        return 0;
    }
    int i;
    for(i=1;i<pa->n;i++){pa->el[i]=-1;}
    for(i=0;i<pb->n;i++){pb->el[i]=2+pow((LOW+(i+1)*h)*h,2)-2*G*pow(phi->el[i]*h,2);}
    for(i=0;i<pc->n-1;i++){pc->el[i]=-1;}
    return 1;
}

int main(){
    FILE* fdata;
    G=atof(G_CHAR);
    char name[50]="../data/extra/cmpVector_G=";
    char fmt[5]=".txt";
    strncat(name,G_CHAR,strlen(G_CHAR));
    strncat(name,fmt,strlen(fmt));
	fdata=fopen(name,"w");
    h=1/(double)NUM;
    eigen* peigen=(eigen*)malloc(sizeof(eigen));
    peigen->value=-999999;
    peigen->pv=createNewVector((UP-LOW)*NUM-1);
    int i;
    for(i=0;i<peigen->pv->n;i++){
        peigen->pv->el[i]=trailFunc(LOW+(i+1)*h);
    } 
    if(normalize(peigen->pv)==0)return 0;
    vector *pa,*pb,*pc;
    pa=createNewVector(peigen->pv->n);
    pb=createNewVector(peigen->pv->n);
    pc=createNewVector(peigen->pv->n);
    if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
    
    eigen* peigen_ex;
    peigen_ex=peigen;
    peigen=inversePowerMethod(pa,pb,pc,PDIS);
    if(normalize(peigen->pv)==0)return 0;

    int count=1;
    int nojump;
    double epsValue=99,epsVector=99;
    double epsValue_ex,epsVector_ex;
   do{
        if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
        freeVector(peigen_ex->pv);
        peigen_ex=peigen;
        peigen=inversePowerMethod(pa,pb,pc,PDIS);
        if(normalize(peigen->pv)==0)return 0;
        count++;
        epsValue_ex=epsValue;epsVector_ex=epsVector;
        epsValue=fabs((peigen->value-peigen_ex->value)/peigen->value);
        epsVector=distance(peigen->pv,peigen_ex->pv);
        printf("%e	%e\n",epsValue,epsVector);
        if(epsValue<EPS_J&&epsVector<EPS_J){
            nojump=1;break;
        }
        if(fabs((epsValue-epsValue_ex)/epsValue)<EPS_J&&fabs((epsVector-epsVector_ex)/epsVector)<EPS_J){
            nojump=0;break;
        }
    } while(1);
    if(nojump){
        fprintf(fdata,"eigenvalue:\n");
        fprintf(fdata,"%e\n",peigen->value/(2*h*h));
        fprintf(fdata,"eigenvector:\n");
        for(i=0;i<peigen->pv->n;i++){
            fprintf(fdata,"%lf  %e\n",LOW+(i+1)*h,peigen->pv->el[i]);
        }
        printf("%d\n",count);
    }
    else{
        fprintf(fdata,"eigenvalue 1:\n");
        fprintf(fdata,"%e\n",peigen->value/(2*h*h));
        fprintf(fdata,"eigenvector 1:\n");
        for(i=0;i<peigen->pv->n;i++){
            fprintf(fdata,"%lf  %e\n",LOW+(i+1)*h,peigen->pv->el[i]);
        }
        printf("%d\n",count);

        if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
        freeVector(peigen_ex->pv);
        peigen_ex=peigen;
        peigen=inversePowerMethod(pa,pb,pc,PDIS);
        if(normalize(peigen->pv)==0)return 0;

        fprintf(fdata,"eigenvalue 2:\n");
        fprintf(fdata,"%e\n",peigen->value/(2*h*h));
        fprintf(fdata,"eigenvector 2:\n");
        for(i=0;i<peigen->pv->n;i++){
            fprintf(fdata,"%lf  %e\n",LOW+(i+1)*h,peigen->pv->el[i]);
        }
        printf("%d\n",count);
    }
    return 0;
}
