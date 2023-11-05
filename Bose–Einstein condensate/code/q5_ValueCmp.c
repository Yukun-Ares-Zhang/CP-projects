//本程序计算提升精度后G取0至-2区间的特征值
//本程序自动判断真收敛或伪收敛
#include "matrix.h"
#include <string.h> 
#define NUM 100
#define LOW -20
#define UP 20
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

    FILE* fdata=fopen("../data/extra/q5_values_cmp.txt","w");
    fprintf(fdata,"G    epsilon\n");
    for(G=0;G>=-2.001;G-=0.1){
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
        int nojump;
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
            printf("%.14e	%.14e\n",epsValue,epsVector);
            if(epsValue<EPS_J&&epsVector<EPS_J){//真收敛
            nojump=1;break;
            }
            if(fabs((epsValue-epsValue_ex)/epsValue)<EPS_J&&fabs((epsVector-epsVector_ex)/epsVector)<EPS_J){//伪收敛
                nojump=0;break;
            }
        } while(1);
        if(nojump){//真收敛
            fprintf(fdata,"%lf  %e\n",G,peigen->value/(2*h*h));
            printf("After %d iterations, G=%lf finished.\n",count,G);
        }
        else{//伪收敛
            fprintf(fdata,"%lf  %e\n",G,peigen->value/(2*h*h));
            
            if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
            freeVector(peigen_ex->pv);
            peigen_ex=peigen;
            peigen=inversePowerMethod(pa,pb,pc,pdis);
            if(normalize(peigen->pv)==0)return 0;
            fprintf(fdata,"%lf  %e\n",G,peigen->value/(2*h*h));
            
            printf("After %d iterations, G=%lf finished.\n",count,G);
        }
    }
    fclose(fdata);
    return 0;
}
