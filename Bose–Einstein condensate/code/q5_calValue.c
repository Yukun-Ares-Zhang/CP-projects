//本程序计算g(n-1)取-2~0时基态epsilon，间隔0.1
//输出格式为"q5_values.txt"
#include "matrix.h"
#include <string.h> 
#define NUM 100
#define LOW -20
#define UP 20
#define EPS_J 1e-5//判断特征值特征向量迭代收敛的特征尺度
double G;
double h;
double trailFunc(double x){//初置波函数
    return 1/pow(PI,0.25)*pow(E,-x*x/2);
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
    for(i=0;i<pb->n;i++){pb->el[i]=2+pow((LOW+(i+1)*h)*h,2)-2*G*pow(phi->el[i]*h,2);}
    for(i=0;i<pc->n-1;i++){pc->el[i]=-1;}
    return 1;
}

int main(){
    h=1/(double)NUM;
    eigen* peigen=(eigen*)malloc(sizeof(eigen));//储存当前特征值特征向量
    eigen* peigen_ex;//储存上一个特征值特征向量
    peigen->value=-999999;
    peigen->pv=createNewVector((UP-LOW)*NUM-1);
    int i;
    //创建并构造系数矩阵，pb是主对角元，pa是下次对角元，pc是上次对角元
    vector *pa,*pb,*pc;
    pa=createNewVector(peigen->pv->n);
    pb=createNewVector(peigen->pv->n);
    pc=createNewVector(peigen->pv->n);

    FILE* fdata=fopen("../data/q5/q5_values.txt","w");
    fprintf(fdata,"G    epsilon\n");
    for(G=0;G>=-2.001;G-=0.1){
        for(i=0;i<peigen->pv->n;i++){//初猜波函数
            peigen->pv->el[i]=trailFunc(LOW+(i+1)*h);
        } 
        if(normalize(peigen->pv)==0)return 0;//归一化
        
        int pdis=0;
        if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
        
        peigen_ex=peigen;
        peigen=inversePowerMethod(pa,pb,pc,pdis);
        if(normalize(peigen->pv)==0)return 0;

        int count=1;
        double epsValue,epsVector;
        //相邻两次迭代特征值地相对偏差和特征向量差的无穷范数
        do{//重复迭代直到相邻两次特征值相对偏差和特征向量距离（无穷范数）小于给定值EPS
            if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;   
            freeVector(peigen_ex->pv);
            peigen_ex=peigen;
            peigen=inversePowerMethod(pa,pb,pc,pdis);
            if(normalize(peigen->pv)==0)return 0;
            count++;
            //计算两次迭代结果的偏差
            epsValue=fabs((peigen->value-peigen_ex->value)/peigen->value);
            epsVector=distance(peigen->pv,peigen_ex->pv);
            printf("%e	%e\n",epsValue,epsVector);//控制台输出便于监测进程
            if(epsValue<EPS_J&&epsVector<EPS_J){
            break;
            }
        } while(1);
        
        fprintf(fdata,"%lf  %e\n",G,peigen->value/(2*h*h));
        printf("After %d iterations, G=%lf finished.\n",count,G);//控制台输出便于监测进程
    }
    fclose(fdata);
    return 0;
}
