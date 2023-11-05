/*
本程序计算第二问，给定一条出射光线，求出之后所有光线与物体（反射面和黑体）的接触点
本程序不用手动设置参数，将自动打出参数取为如下六种情况时对应的txt文件，在"data/q2"文件夹中
其中STARTBY描述光线从A/B哪个球出发，（THETA1,THETA2,PHI!,PHI2）表示出射光线方位角
         plot1       plot2       plot3       plot4       plot5       plot6
THETA1   1.          0.5         1.5         0.          1.5         1.57
THETA2   0.3         0.5         1.          1.5         1.5         1.5
PHI1     0.3         3.          3.          3.          3.          3.
PHI2     0.3         0.3         3.          3.          3.          3.
STARTBY  'A'         'B'         'A'         'B'         'A'         'A'
*/
#include "HW7.h"
#define RA 0.2
#define RB 0.2

double THETA1;
double THETA2;
double PHI1;
double PHI2;
char STARTBY;

double const1[6]={1,0.5,1.5,0.,1.5,1.57};
double const2[6]={0.3,0.5,1.,1.5,1.5,1.5};
double const3[6]={0.3,3.,3.,3.,3.,3.};
double const4[6]={0.3,0.3,3.,3.,3.,3.};
char const5[6]={'A','B','A','B','A','A'};

double a1,c,a2,BOUND0,BOUND1,b1,b2;
transition* pt;
vector* rA,* rB;
int initialize(void);
transition* oneTransition(vector* r0, vector* n0);


int main(){
    int i;
    for(i=0;i<6;i++){
        //生成文件名
        char name[50]="../data/q2/q2_plot";
        char num[5]="";
        snprintf(num,3,"%d",i+1);
        char fmt[5]=".txt";
        strncat(name,num,strlen(num));
        strncat(name,fmt,strlen(fmt));
        FILE *fdata=fopen(name,"w");

        THETA1=const1[i];
        THETA2=const2[i];
        PHI1=const3[i];
        PHI2=const4[i];
        STARTBY=const5[i];
        initialize();
        int count=0;
        while(1){
            if(pt->situation==0){break;}
            if(pt->situation==1){break;}
            print_transition(pt);
            fprint_vector(pt->r,fdata);//打印交点
            pt=oneTransition(pt->r,pt->n);
            if(pt->situation==2)count++;
        }
        print_transition(pt);
        fprint_vector(pt->r,fdata);//打印交点
        printf("%d reflection(s)\n",count);
    }
    return 0;
}

int initialize(void){
    rA=createNewVector(3);
    rB=createNewVector(3);
    if(rA==NULL||rB==NULL)return 0;
    rA->el[0]=-2.25;rB->el[0]=2.25;
    a1=2.5;
    c=2.25;
    a2=171./44;
    b1=sqrt(19)/4;
    b2=9*sqrt(15)/11;
    BOUND0=5./18;
    BOUND1=342./121;
    
    pt=(transition*)malloc(sizeof(transition));
    if(pt==NULL){
        printf("Unable to create pt\n");
        return 0;
    }
    pt->situation=2;
    pt->r=createNewVector(3);
    pt->n=createNewVector(3);
    if(pt->r==NULL||pt->n==NULL){
        printf("Unable to create pt->t&pt->n\n");
        return 0;
    }
    
    vector* tau[3];
    int i;
    for(i=0;i<3;i++){
        tau[i]=createNewVector(3);
        if(tau[i]==NULL){
            printf("Unable to create tau[%d]",i);
            return 0;
        }
    }
    tau[0]->el[0]=cos(THETA1)*cos(PHI1);tau[0]->el[1]=cos(THETA1)*sin(PHI1);tau[0]->el[2]=-sin(THETA1);
    tau[1]->el[0]=-sin(PHI1);tau[1]->el[1]=cos(PHI1);
    tau[2]->el[0]=sin(THETA1)*cos(PHI1);tau[2]->el[1]=sin(THETA1)*sin(PHI1);tau[2]->el[2]=cos(THETA1);
    if(STARTBY=='A'){
        scalarMul(RA,tau[2],pt->r);
        vectorAdd(rA,pt->r,pt->r);
    }
    else{
        scalarMul(RB,tau[2],pt->r);
        vectorAdd(rB,pt->r,pt->r);
    }
    scalarMul(sin(THETA2)*cos(PHI2),tau[0],tau[0]);
    scalarMul(sin(THETA2)*sin(PHI2),tau[1],tau[1]);
    scalarMul(cos(THETA2),tau[2],tau[2]);
    for(i=0;i<3;i++){
        vectorAdd(tau[i],pt->n,pt->n);
    }
    return 1;
}
//给定一个初始坐标和方向，计算光线下一次与物质（反射面或黑体球）接触点，并判断被吸收还是被反射，若吸收，记录被谁吸收；若反射，给出反射后的方向向量。
transition* oneTransition(vector* r0, vector* n0){
    transition *pt=(transition*)malloc(sizeof(transition));
    if(r0->n!=3||n0->n!=3){
        pt->situation=0;
        printf("Wrong input vector\n");
        return pt;
    }
    if(normalize(n0)==0){//把n0归一化
        pt->situation=0;
        printf("Unable to normalize n0\n");
        return pt;
    }
    double t[5]={999,999,999,999,999};
    double A,B,C;
    int i;
    double time1,time2,x1,x2;
    int inRange1,inRange2;
    for(i=0;i<5;i++){
        switch(i){
            case 0: //大球交点t[0]
                A=1;
                B=2*scalarDot(r0,n0)-2*c*n0->el[0];
                C=scalarDot(r0,r0)-2*c*r0->el[0];
                //printf("%lf\n",scalarDot(r0,r0));
                //printf("A=%lf,B=%lf,C=%lf\n",A,B,C);
                if(B*B-4*A*C<0)break;
                time1=(-B+sqrt(B*B-4*A*C))/(2*A);
                time2=(-B-sqrt(B*B-4*A*C))/(2*A);
                inRange1=0;inRange2=0;
                x1=r0->el[0]+time1*n0->el[0];
                x2=r0->el[0]+time2*n0->el[0];
                if(time1>EPS&&x1<BOUND1&&x1>=BOUND0)inRange1=1;
                if(time2>EPS&&x2<BOUND1&&x2>=BOUND0)inRange2=1;
                if(inRange1&&inRange2){
                    t[i]=min(time1,time2);
                }
                else if(inRange1){
                    t[i]=time1;
                }
                else if(inRange2){
                    t[i]=time2;
                }
                break;
            case 1: //小椭球球交点t[1]
                A=pow(n0->el[0],2)/(a1*a1)+pow(n0->el[1],2)/(b1*b1)+pow(n0->el[2],2)/(b1*b1);
                B=2*(n0->el[0]*r0->el[0]/(a1*a1)+n0->el[1]*r0->el[1]/(b1*b1)+n0->el[2]*r0->el[2]/(b1*b1));
                C=pow(r0->el[0],2)/(a1*a1)+pow(r0->el[1],2)/(b1*b1)+pow(r0->el[2],2)/(b1*b1)-1;
                if(B*B-4*A*C<0)break;
                time1=(-B+sqrt(B*B-4*A*C))/(2*A);
                time2=(-B-sqrt(B*B-4*A*C))/(2*A);
                inRange1=0;inRange2=0;
                x1=r0->el[0]+time1*n0->el[0];
                x2=r0->el[0]+time2*n0->el[0];
                if(time1>EPS&&x1<BOUND0)inRange1=1;
                if(time2>EPS&&x2<BOUND0)inRange2=1;
                if(inRange1&&inRange2){
                    t[i]=min(time1,time2);
                }
                else if(inRange1){
                    t[i]=time1;
                }
                else if(inRange2){
                    t[i]=time2;
                }
                break;
            case 2: //大椭球球交点t[2]
                A=pow(n0->el[0],2)/(a2*a2)+pow(n0->el[1],2)/(b2*b2)+pow(n0->el[2],2)/(b2*b2);
                B=2*(n0->el[0]*r0->el[0]/(a2*a2)+n0->el[1]*r0->el[1]/(b2*b2)+n0->el[2]*r0->el[2]/(b2*b2));
                C=pow(r0->el[0],2)/(a2*a2)+pow(r0->el[1],2)/(b2*b2)+pow(r0->el[2],2)/(b2*b2)-1;
                if(B*B-4*A*C<0)break;
                time1=(-B+sqrt(B*B-4*A*C))/(2*A);
                time2=(-B-sqrt(B*B-4*A*C))/(2*A);
                inRange1=0;inRange2=0;
                x1=r0->el[0]+time1*n0->el[0];
                x2=r0->el[0]+time2*n0->el[0];
                if(time1>EPS&&x1>=BOUND1)inRange1=1;
                if(time2>EPS&&x2>=BOUND1)inRange2=1;
                if(inRange1&&inRange2){
                    t[i]=min(time1,time2);
                }
                else if(inRange1){
                    t[i]=time1;
                }
                else if(inRange2){
                    t[i]=time2;
                }
                break;
            case 3: //球A交点t[3]
                A=1;
                B=2*scalarDot(r0,n0)+2*c*n0->el[0];
                C=scalarDot(r0,r0)+2*c*r0->el[0]+c*c-RA*RA;
                if(B*B-4*A*C<0)break;
                time1=(-B+sqrt(B*B-4*A*C))/(2*A);
                time2=(-B-sqrt(B*B-4*A*C))/(2*A);
                inRange1=0;inRange2=0;
                if(time1>EPS)inRange1=1;
                if(time2>EPS)inRange2=1;
                if(inRange1&&inRange2){
                    t[i]=min(time1,time2);
                }
                else if(inRange1){
                    t[i]=time1;
                }
                else if(inRange2){
                    t[i]=time2;
                }
                break;
            default: //球B交点t[4]
                A=1;
                B=2*scalarDot(r0,n0)-2*c*n0->el[0];
                C=scalarDot(r0,r0)-2*c*r0->el[0]+c*c-RB*RB;
                if(B*B-4*A*C<0)break;
                time1=(-B+sqrt(B*B-4*A*C))/(2*A);
                time2=(-B-sqrt(B*B-4*A*C))/(2*A);
                inRange1=0;inRange2=0;
                if(time1>EPS)inRange1=1;
                if(time2>EPS)inRange2=1;
                if(inRange1&&inRange2){
                    t[i]=min(time1,time2);
                }
                else if(inRange1){
                    t[i]=time1;
                }
                else if(inRange2){
                    t[i]=time2;
                }
                break;
        }
    }
    int object=0;
    double time=t[0];
    for(i=1;i<5;i++){//取所有交点时刻中最小的一个是真实相交时间，并记录这个物体的编号
        if(t[i]<time){
            time=t[i];
            object=i;
        }
    }
    if(time>100){//时间过长，出现错误
        pt->situation=0;
        printf("Wrong time\n");
        return pt;
    }
    //printf("object=%d,time=%lf\n",object,time); 
    pt->r=createNewVector(3);
    if(pt->r==NULL){
        pt->situation=0;
        printf("Unable to create r\n");
        return pt;
    }
    scalarMul(time,n0,pt->r);//计算接触点坐标
    vectorAdd(r0,pt->r,pt->r);
    if(object==0||object==1||object==2){//反射
        pt->situation=2;
        pt->n=createNewVector(3);//计算出射方向
        if(pt->n==NULL){
            pt->situation=0;
            printf("Unable to create n\n");
            return pt;
        }
        if(object==0){
            pt->n->el[0]=2*(pt->r->el[0]-c);
            pt->n->el[1]=2*pt->r->el[1];
            pt->n->el[2]=2*pt->r->el[2];
        }
        if(object==1){
            pt->n->el[0]=2*pt->r->el[0]/(a1*a1);
            pt->n->el[1]=2*pt->r->el[1]/(b1*b1);
            pt->n->el[2]=2*pt->r->el[2]/(b1*b1);
        }
        if(object==2){
            pt->n->el[0]=2*pt->r->el[0]/(a2*a2);
            pt->n->el[1]=2*pt->r->el[1]/(b2*b2);
            pt->n->el[2]=2*pt->r->el[2]/(b2*b2);
        }
        if(normalize(pt->n)==0){
            pt->situation=0;
            printf("Unable to normalize n\n");
            return pt;
        }
        scalarMul(-2*scalarDot(n0,pt->n),pt->n,pt->n);
        vectorAdd(n0,pt->n,pt->n);
        return pt;
    }
    else{//吸收
        pt->situation=1;
        if(object==3)pt->absorb='A';//记录被谁吸收
        if(object==4)pt->absorb='B';
        return pt;
    }
}
