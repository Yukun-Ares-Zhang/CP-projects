#include "matrix.h"
#define NUM 32
#define alpha 0.25
#define SCALE 1e+2

double dt;

//butcher表
double c[3]={0},b[3]={0},bbar[3],a[3][3]={0};

void create_Butcher(){
    int i;
    for(i=0;i<3;i++){
        c[i]=(3+pow(-1,i)*sqrt(3))/6;
    }
    a[1][0]=(2-sqrt(3))/12;a[2][1]=sqrt(3)/6;
    bbar[0]=(5-3*sqrt(3))/24;bbar[1]=(3+sqrt(3))/12;bbar[2]=(1+sqrt(3))/24;
    b[0]=(3-2*sqrt(3))/12;b[1]=0.5;b[2]=(3+2*sqrt(3))/12;
}

int add_vetor(vector* src1,vector* src2,vector* des){
    if(src1->n!=src2->n||src2->n!=des->n||des->n<0)return 0;
    int i;
    for(i=0;i<des->n;i++){
        des->el[i]=src1->el[i]+src2->el[i];
    }
    return 1;
}

int scalar_mul_vector(vector* src,double k,vector* des){
    if(src->n!=des->n||des->n<0)return 0;
    int i;
    for(i=0;i<des->n;i++){
        des->el[i]=k*src->el[i];
    }
    return 1;
}
int f(vector* g,vector* fg){
    if(g->n!=fg->n||fg->n<0)return 0;
    int j;
    for(j=0;j<fg->n;j++){
        if(j==0){
            fg->el[j] = -(g->el[j]-g->el[j+1])-g->el[j]-alpha*pow((g->el[j]-g->el[j+1]),2)+alpha*pow((g->el[j]),2);
            continue;
        }
        if(j==fg->n-1){
            fg->el[j] = -(g->el[j])+(g->el[j-1]-g->el[j])-alpha*pow((g->el[j]),2)+alpha*pow((g->el[j-1]-g->el[j]),2);
            continue;
        }
        fg->el[j] = -(g->el[j]-g->el[j+1])+(g->el[j-1]-g->el[j])-alpha*pow((g->el[j]-g->el[j+1]),2)+alpha*pow((g->el[j-1]-g->el[j]),2);
    }
    return 1;
}
int onestep_RKN_rank4(vector* q,vector* p,double h){
    if(q->n!=p->n)return 0;
    int i,j;
    vector* g[3],* fg[3];
    vector* temp=createNewVector(q->n);
    for(i=0;i<3;i++){
        g[i]=createNewVector(q->n);
        fg[i]=createNewVector(q->n);
        if(g[i]==NULL||fg[i]==NULL)return 0;
    }
    //生成所有g[i],fg[i]
    for(i=0;i<3;i++){
        for(j=0;j<i;j++){
            scalar_mul_vector(fg[i],a[i][j]*h*h,temp);
            add_vetor(g[i],temp,g[i]);
        }
        scalar_mul_vector(p,c[i]*h,temp);
        add_vetor(g[i],temp,g[i]);
        add_vetor(g[i],q,g[i]);
        f(g[i],fg[i]);
    }
    //计算新q
    scalar_mul_vector(p,h,temp);
    add_vetor(q,temp,q);
    for(i=0;i<3;i++){
        scalar_mul_vector(fg[i],bbar[i]*h*h,temp);
        add_vetor(q,temp,q);
    }
    //计算新p
    for(i=0;i<3;i++){
        scalar_mul_vector(fg[i],b[i]*h,temp);
        add_vetor(p,temp,p);
    }
    //释放中间运算向量
    for(i=0;i<3;i++){
        free(g[i]);
        free(fg[i]);
    }
    free(temp);
    return 1;
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
    create_Butcher();
    dt=PI/sin(PI/(2*(NUM+1)))/SCALE;
    FILE* fdata=fopen("C:/C_files/HW4/q2_RKN4.txt","w");
    SqMatrix* A=createZeroSqMatrix(NUM);
    int i,j,k;
    for ( k = 0; k < NUM; k++)
    {
        for ( j = 0; j < NUM; j++)
        {
            A->el[k][j]=sqrt(2./NUM)*sin(PI*(k+1)*(j+1)/(NUM+1));
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

    print_vector(p);

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
    
    
    for(i=1;i<=5;i++){
        onestep_RKN_rank4(q,p,dt);
        print_vector(q);
        print_vector(p);
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
