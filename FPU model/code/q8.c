#include "matrix.h"
#define NUM 128
#define BETA 0
#define K 11
#define B 0.5
#define SCALE 1e+3



double dt;
double wk;

double dqH(int j,vector* q){
    if(j<0||j>NUM-1){printf("wrong index: the index of the vector entry should be in range 0 to NUM-1\n");return 0;}
    switch (j)
    {
    case 0:
        return -(q->el[j]-q->el[j+1])-q->el[j]-BETA*pow((q->el[j]-q->el[j+1]),3)+BETA*pow((-q->el[j]),3);
        break;

    case NUM-1:
        return -(q->el[j])+(q->el[j-1]-q->el[j])-BETA*pow((q->el[j]),3)+BETA*pow((q->el[j-1]-q->el[j]),3);
        break;
    
    default:
        return -(q->el[j]-q->el[j+1])+(q->el[j-1]-q->el[j])-BETA*pow((q->el[j]-q->el[j+1]),3)+BETA*pow((q->el[j-1]-q->el[j]),3);
        break;
    }
}

void fprint_vector(FILE* fdata,vector* b){
    int i;
    for ( i = 0; i < b->n; i++)
    {
        fprintf(fdata,"%lf",b->el[i]);
        if(i!=b->n-1)fprintf(fdata,"\t");
    }
    fprintf(fdata,"\n");
}

int main(){
    wk=2*sin(PI*K/(2*(NUM+1)));
    dt=2*PI/wk/SCALE;
    FILE* fdata=fopen("C:/C_files/HW4/q8_beta=0.txt","w");
    int i,j,k;
    vector* q=createNewVector(NUM);    
    vector* p=createNewVector(NUM);

    for ( i = 0; i < NUM; i++)
    {
        q->el[i]=B*cos(PI*K*(i+1-NUM/2)/(NUM+1))/cosh(sqrt(3./2)*B*wk*(i+1-NUM/2));
        p->el[i]=B/cosh(sqrt(3./2)*B*wk*(i+1-NUM/2))*(wk*(1+3./16*wk*wk*B*B)*sin(PI*K*(i+1-NUM/2)/(NUM+1))+sqrt(3./2)*B*cos(PI*K*(i+1-NUM/2)/(NUM+1))*sin(PI*K/(NUM+1))*tanh(sqrt(3./2)*B*wk*(i+1-NUM/2)));
    }
   
    fprintf(fdata,"0\t");
    fprint_vector(fdata,q);  

    vector* pf, *pr;
    pr=createNewVector(NUM);
    pf=createNewVector(NUM);
    for ( j = 0; j < NUM; j++)
    {
        pf->el[j]=p->el[j]+0.5*dt*dqH(j,q);
    }
    
    for(i=1;i<=200*SCALE;i++){
        copy_vector(pf,pr);
        for ( j = 0; j < NUM; j++)
        {
            q->el[j]+=(pr->el[j]*dt);
        }
        for ( j = 0; j < NUM; j++){
            pf->el[j]=pr->el[j]+dt*dqH(j,q);
        }
        average_vector(pr,pf,p);
        if(i%100==0){
            fprintf(fdata,"%lf\t",(double)i/SCALE);
            fprint_vector(fdata,q);  
        }
    }
    return 0;
}
