//本程序计算相邻两次迭代的特征值相对偏差和波函数的距离(epsValue,epsVector)在迭代过程中的演化
//已经设定G=-2.0
#include "matrix.h"
#include <string.h> 
#define NUM 100
#define LOW -20
#define UP 20
#define MAX 200
#define EPS_J 1e-12

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
    G=-2;
    h=1/(double)NUM;
    eigen* peigen=(eigen*)malloc(sizeof(eigen));
    eigen* peigen_ex;
    peigen->value=-999999;
    peigen->pv=createNewVector((UP-LOW)*NUM-1);
    int i;
    
    vector *pa,*pb,*pc;
    pa=createNewVector(peigen->pv->n);
    pb=createNewVector(peigen->pv->n);
    pc=createNewVector(peigen->pv->n);

    FILE* fdata=fopen("../data/extra/plot_iteration_G=-2.txt","w");
    fprintf(fdata,"epsValue    epsVector\n");
    
    for(i=0;i<peigen->pv->n;i++){
        peigen->pv->el[i]=trailFunc(LOW+(i+1)*h);
    } 
    if(normalize(peigen->pv)==0)return 0;
    
    int pdis=0;
    if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
    
    peigen_ex=peigen;
    peigen=inversePowerMethod(pa,pb,pc,pdis);
    if(normalize(peigen->pv)==0)return 0;

    int count=1;
    double epsValue=99,epsVector=99;
    double epsValue_ex,epsVector_ex;

    do{
        if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
        freeVector(peigen_ex->pv);
        peigen_ex=peigen;
        peigen=inversePowerMethod(pa,pb,pc,pdis);
        if(normalize(peigen->pv)==0)return 0;
        count++;
        epsValue_ex=epsValue;epsVector_ex=epsVector;
        epsValue=fabs((peigen->value-peigen_ex->value)/peigen->value);
        epsVector=distance(peigen->pv,peigen_ex->pv);
        fprintf(fdata,"%e	%e\n",epsValue,epsVector);
        if(epsValue<EPS_J&&epsVector<EPS_J){
        break;
        }
        if(count>MAX){
            break;
        }
    } while(1);
       
    fclose(fdata);
    return 0;
}
