#include "matrix.h"
#define NUM 32
#define alpha 0.25
#define SCALE 1e+3

double dt;

double dqH(int j,vector* q){
    if(j<0||j>NUM-1){printf("wrong index: the index of the vector entry should be in range 0 to NUM-1\n");return 0;}
    switch (j)
    {
    case 0:
        return -(q->el[j]-q->el[j+1])-q->el[j]-alpha*pow((q->el[j]-q->el[j+1]),2)+alpha*pow((q->el[j]),2);
        break;

    case NUM-1:
        return -(q->el[j])+(q->el[j-1]-q->el[j])-alpha*pow((q->el[j]),2)+alpha*pow((q->el[j-1]-q->el[j]),2);
        break;
    
    default:
        return -(q->el[j]-q->el[j+1])+(q->el[j-1]-q->el[j])-alpha*pow((q->el[j]-q->el[j+1]),2)+alpha*pow((q->el[j-1]-q->el[j]),2);
        break;
    }
}

double eigen_transform(SqMatrix* A,vector* b,int k){
	if(A->n!=b->n)return 0;
    int i;
    double s=0;
    for ( i = 0; i < b->n ; i++)
    {
        s+=A->el[k][i]*b->el[i];
    }
    return s;
}

int main(){
    dt=PI/sin(PI/(2*(NUM+1)))/SCALE;
    FILE* fdata=fopen("C:/C_files/HW4/q5_cmp_sparsify.txt","w");
    SqMatrix* A=createZeroSqMatrix(NUM);
    int i,j,k;
    for ( k = 0; k < NUM; k++)
    {
        for ( j = 0; j < NUM; j++)
        {
            A->el[k][j]=sqrt(2./(NUM+1))*sin(PI*(k+1)*(j+1)/(NUM+1));
        }
    }
    SqMatrix* Acopy=createZeroSqMatrix(NUM);
    copy_SqMatrix(A,Acopy);
	
    vector* q=createNewVector(NUM);
    q->el[0]=4;
    if(Solve_RankN_LinearEq(A,q)){
        print_vector(q);
    }
    vector* p=createNewVector(NUM);

    fprintf(fdata,"0\t");
    
    for ( i=0;i<4;i++)
    {
        double Qk=eigen_transform(Acopy,q,i);
        double Pk=eigen_transform(Acopy,p,i);
        double wk=2*sin((i+1)*PI/2/(NUM+1));
        double Ek=0.5*Pk*Pk+0.5*pow(wk*Qk,2);
        fprintf(fdata,"%lf\t",Ek);
        
    }
    fprintf(fdata,"\n");
    

    vector* pf, *pr;
    pr=createNewVector(NUM);
    pf=createNewVector(NUM);
    for ( j = 0; j < NUM; j++)
    {
        pf->el[j]=0.5*dt*dqH(j,q);
    }
    
    for(i=1;i<=160*SCALE;i++){
        copy_vector(pf,pr);
        for ( j = 0; j < NUM; j++)
        {
            q->el[j]+=(pr->el[j]*dt);
        }
        for ( j = 0; j < NUM; j++){
            pf->el[j]=pr->el[j]+dt*dqH(j,q);
        }
        average_vector(pr,pf,p);
        fprintf(fdata,"%lf\t",(double)i/SCALE);
        for ( k=0;k<4;k++)
        {
            double Qk=eigen_transform(Acopy,q,k);
            double Pk=eigen_transform(Acopy,p,k);
            double wk=2*sin((k+1)*PI/2/(NUM+1));
            double Ek=0.5*Pk*Pk+0.5*pow(wk*Qk,2);
            fprintf(fdata,"%lf\t",Ek);
        }
        fprintf(fdata,"\n");
    }
    return 0;
}
