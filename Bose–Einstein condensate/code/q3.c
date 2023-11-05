//本程序计算第三问
//输出为“q3.txt”
#include "matrix.h"
#define NUM 100
#define LOW -20
#define UP 20
#define PDIS 0.1
double h;
double trailFunc(double x){//初置波函数
    return 1/pow(10*PI,0.25)*pow(E,-x*x/20);
}
int normalize(vector* pv){//归一化函数，采用中点法数值求积
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
int resetCoeff(vector *phi,vector *pa,vector *pb,vector *pc){//重置系数矩阵
    if(phi==NULL||pa==NULL||pb==NULL||pc==NULL)return 0;
    if(phi->n!=pa->n||pa->n!=pb->n||pb->n!=pc->n){
        printf("Coefficient reset error\n");
        return 0;
    }
    int i;
    for(i=1;i<pa->n;i++){pa->el[i]=-1;}
    for(i=0;i<pb->n;i++){pb->el[i]=2*(1-pow(phi->el[i]*h,2));}
    for(i=0;i<pc->n-1;i++){pc->el[i]=-1;}
    return 1;
}
int main(){
    FILE* fdata=fopen("../data/q3/q3.txt","w");
    h=1/(double)NUM;
    eigen* peigen=(eigen*)malloc(sizeof(eigen));//储存当前特征值特征向量
    peigen->value=-999999;
    peigen->pv=createNewVector((UP-LOW)*NUM-1);
    int i;
    for(i=0;i<peigen->pv->n;i++){//初猜波函数
        peigen->pv->el[i]=trailFunc(LOW+(i+1)*h);
    } 
    if(normalize(peigen->pv)==0)return 0;//归一化
    //创建并构造系数矩阵，pb是主对角元，pa是下次对角元，pc是上次对角元
    vector *pa,*pb,*pc;
    pa=createNewVector(peigen->pv->n);
    pb=createNewVector(peigen->pv->n);
    pc=createNewVector(peigen->pv->n);
    if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
    
    eigen* peigen_ex;//储存上一个特征值特征向量
    peigen_ex=peigen;
    peigen=inversePowerMethod(pa,pb,pc,PDIS);
    if(normalize(peigen->pv)==0)return 0;//对新特征向量归一化
    
    double epsValue,epsVector;//相邻两次迭代特征值地相对偏差和特征向量差的无穷范数
   do{//重复迭代直到相邻两次特征值相对偏差和特征向量距离（无穷范数）小于给定值EPS
        if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
        freeVector(peigen_ex->pv);
        peigen_ex=peigen;
        peigen=inversePowerMethod(pa,pb,pc,PDIS);
        if(normalize(peigen->pv)==0)return 0;
        //计算两次迭代结果的偏差
        epsValue=fabs((peigen->value-peigen_ex->value)/peigen->value);
        epsVector=distance(peigen->pv,peigen_ex->pv);
        printf("%e	%e\n",epsValue,epsVector);//控制台输出便于监测进程
    } while(epsValue>EPS||epsVector>EPS);
    //对求出的波函数积分算总能量
    double energy=0;
    for(i=0;i<peigen->pv->n;i++){
        if(i==0)energy+=peigen->pv->el[i]*h*((peigen->pv->el[i+1]-2*peigen->pv->el[i])/(h*h)+pow(peigen->pv->el[i],3));
        else if(i==peigen->pv->n-1)energy+=peigen->pv->el[i]*h*((peigen->pv->el[i-1]-2*peigen->pv->el[i])/(h*h)+pow(peigen->pv->el[i],3));
        else{
            energy+=peigen->pv->el[i]*h*((peigen->pv->el[i-1]+peigen->pv->el[i+1]-2*peigen->pv->el[i])/(h*h)+pow(peigen->pv->el[i],3));
        }
    }
    fprintf(fdata,"total energy: %e\n",energy);
    fprintf(fdata,"eigenvalue:\n");
    fprintf(fdata,"%e\n",peigen->value/(2*h*h));
    fprintf(fdata,"eigenvector:\n");
    for(i=0;i<peigen->pv->n;i++){
        fprintf(fdata,"%lf  %e\n",LOW+(i+1)*h,peigen->pv->el[i]);
    }
    
    return 0;
}
