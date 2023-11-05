
//本程序进行第三问到第五问的数值试验，总试验次数设置为20次，需要手动调置A和B的球半径
#include "HW7.h"

#define RA 0.005
#define RB 0.005
#define INC_NUM 10000//单次试验光线数
#define RFL_MAX 1e+10//最大反射次数
#define EXP_NUM 20//试验数

double theta1,theta2,phi1,phi2;
char startBy;
double a1,c,a2,BOUND0,BOUND1,b1,b2;
transition* pt;
vector* rA,* rB,* tau[3];

int initialize(void);
transition* oneTransition(vector* r0, vector* n0);


int main(){
     //实现自动命名“RA=XXX&RB=XXX.txt”格式
    char name[50]="../data/improved/RA=";
    char adjoint[5]="&RB=";
    char fmt[5]=".txt";
    char ra_str[7]={0};
    char rb_str[7]={0};
    snprintf(ra_str,5,"%lf",RA);
    snprintf(rb_str,5,"%lf",RB);
    strncat(name,ra_str,strlen(ra_str));
    strncat(name,adjoint,strlen(adjoint));
    strncat(name,rb_str,strlen(rb_str));
    strncat(name,fmt,strlen(fmt));
    FILE *fdata=fopen(name,"w");
    
    srand((unsigned)time(NULL));
    int j;
    for(j=0;j<EXP_NUM;j++){
    	int AtoA=0,AtoB=0,BtoA=0,BtoB=0,Atmr=0,Btmr=0;
    int reachA,reachB,tooManyRfl;
	startBy='A';
    recal:reachA=0;
	reachB=0;
	tooManyRfl=0;
    int i;
    for(i=0;i<INC_NUM;i++){
    	initialize();
        int count_rfl=0;
        while(1){
            if(pt->situation==0){printf("Reflection error\n");break;}
            if(pt->situation==1){break;}
            if(count_rfl>100000){tooManyRfl++;break;}
            pt=oneTransition(pt->r,pt->n);
            if(pt->situation==2)count_rfl++;
        }
        if(pt->situation==1){
            if(pt->absorb=='A')reachA++;
            if(pt->absorb=='B')reachB++;
        }
        //printf("%d\n",i);
    }
    //printf("Start from %c:\nTotal:%d\nReach A:%d\nReach B:%d\nToo many reflections:%d\n",startBy,RFL_NUM,reachA,reachB,tooManyRfl);
    if(startBy=='A'){
        AtoA=reachA;AtoB=reachB;Atmr=tooManyRfl;
    }
    else{
        BtoA=reachA;BtoB=reachB;Btmr=tooManyRfl;
    }
    if(startBy=='A'){
        startBy='B';
        goto recal;
    }
    double etaT,sigma;//计算TA/TB及其不确定度
    etaT=sqrt(RB/RA)*pow(((double)BtoA/(INC_NUM-Btmr))/((double)AtoB/(INC_NUM-Atmr)),0.25);
    sigma=etaT*sqrt((INC_NUM-BtoA)/(16.*INC_NUM*BtoA)+(INC_NUM-AtoB)/(16.*INC_NUM*AtoB));

    fprintf(fdata,"%lf\t%lf\n",etaT,sigma);
    printf("%d\n",j);//控制台输出监视进程
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
    //printf("pt=%#x\n",pt); 
    pt->r=createNewVector(3);
    pt->n=createNewVector(3);
    if(pt->r==NULL||pt->n==NULL){
        printf("Unable to create pt->t&pt->n\n");
        return 0;
    }
    //printf("pt->r=%#x,pt->n=%#x\n",pt->r,pt->n);
    int i;
    for(i=0;i<3;i++){
        tau[i]=createNewVector(3);
        if(tau[i]==NULL){
            printf("Unable to create tau[%d]",i);
            return 0;
        }
    }
	//printf("pt=%#x\n",pt);
    //生成随机方向
    theta1=rand()/((double)RAND_MAX);theta1=acos(1-2*theta1);
    phi1=rand()/((double)RAND_MAX+1);phi1=2*PI*phi1;
    theta2=rand()/((double)RAND_MAX);theta2=acos(1-2*theta2)/2;
    phi2=rand()/((double)RAND_MAX+1);phi2=2*PI*phi2;
    //生成初始坐标及初始方向
    tau[0]->el[0]=cos(theta1)*cos(phi1);tau[0]->el[1]=cos(theta1)*sin(phi1);tau[0]->el[2]=-sin(theta1);
    tau[1]->el[0]=-sin(phi1);tau[1]->el[1]=cos(phi1);
    tau[2]->el[0]=sin(theta1)*cos(phi1);tau[2]->el[1]=sin(theta1)*sin(phi1);tau[2]->el[2]=cos(theta1);
    //printf("pt->r=%#x\n",pt->r);
    if(startBy=='A'){
        scalarMul(RA,tau[2],pt->r);
        vectorAdd(rA,pt->r,pt->r);
    }
    else{
        scalarMul(RB,tau[2],pt->r);
        vectorAdd(rB,pt->r,pt->r);
    }
    //printf("pt->n=%#x\n",pt->n);
    scalarMul(sin(theta2)*cos(phi2),tau[0],pt->n);
    scalarMul(sin(theta2)*sin(phi2),tau[1],tau[1]);
    scalarMul(cos(theta2),tau[2],tau[2]);
    for(i=1;i<3;i++){
        vectorAdd(tau[i],pt->n,pt->n);
    }
    pt->situation=2;
    return 1;
}

transition* oneTransition(vector* r0, vector* n0){
    transition *pt=(transition*)malloc(sizeof(transition));
    if(r0->n!=3||n0->n!=3){
        pt->situation=0;
        printf("Wrong input vector\n");
        return pt;
    }
    if(normalize(n0)==0){
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
    for(i=1;i<5;i++){
        if(t[i]<time){
            time=t[i];
            object=i;
        }
    }
    if(time>100){
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
    scalarMul(time,n0,pt->r);
    vectorAdd(r0,pt->r,pt->r);
    if(object==0||object==1||object==2){
        pt->situation=2;
        pt->n=createNewVector(3);
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
    else{
        pt->situation=1;
        if(object==3)pt->absorb='A';
        if(object==4)pt->absorb='B';
        return pt;
    }
}
