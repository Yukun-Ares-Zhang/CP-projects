//本程序计算第四问g(n-1)=0,1,4,9,16情况的基态波函数
//输出格式为"q4_G=XXX.txt"
#include "matrix.h"
#include <string.h> 
#define NUM 100
#define LOW -20
#define UP 20

int G;
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
    for(i=0;i<peigen->pv->n;i++){//初猜波函数
        peigen->pv->el[i]=trailFunc(LOW+(i+1)*h);
    } 
    if(normalize(peigen->pv)==0)return 0;
    //创建并构造系数矩阵，pb是主对角元，pa是下次对角元，pc是上次对角元
    vector *pa,*pb,*pc;
    pa=createNewVector(peigen->pv->n);
    pb=createNewVector(peigen->pv->n);
    pc=createNewVector(peigen->pv->n);

    //对每个G值，创建一个txt文件
    FILE* fdata;
    char name[50]="";
    char prefix[50]="../data/q4/q4_G=";
    char fmt[5]=".txt";
    char num[5]="";
    int j;
    for(j=0;j<5;j++){
        G=j*j;//G=0,1,4,9,16
        memset(name,0,strlen(name));
        strncpy(name,prefix,strlen(prefix));
        itoa(G,num,10);
        strncat(name,num,strlen(num));
        strncat(name,fmt,strlen(fmt));
        fdata=fopen(name,"w");
        
        int pdis=2*G;//经验取法，由于特征值是负数，需要取一个位移
        if(resetCoeff(peigen->pv,pa,pb,pc)==0)return 0;    
        peigen_ex=peigen;
        peigen=inversePowerMethod(pa,pb,pc,pdis);
        if(normalize(peigen->pv)==0)return 0;

        int count=1;
        double epsValue,epsVector;//相邻两次迭代特征值地相对偏差和特征向量差的无穷范数
    
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
            if(epsValue<EPS&&epsVector<EPS){
            break;
            }
        } while(1);
        fprintf(fdata,"eigenvalue:\n");
        fprintf(fdata,"%e\n",peigen->value/(2*h*h));
        fprintf(fdata,"eigenvector:\n");
        for(i=0;i<peigen->pv->n;i++){
            fprintf(fdata,"%lf  %e\n",LOW+(i+1)*h,peigen->pv->el[i]);
        }
        printf("After %d iterations, G=%d finished.\n",count,G);//控制台输出便于监测进程
        fclose(fdata);
    }
    return 0;
}
