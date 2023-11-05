#ifndef __VECTOR_H__
#define __VECTOR_H__

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "constants.h"

typedef struct vector{
    int n;
    double* el;
}vector;

vector* createNewVector(int n){
    vector* pv;
    pv=(vector*)malloc(sizeof(vector));
    if(pv==NULL)return NULL;
    pv->n=n;
    pv->el=(double*)malloc(sizeof(double)*n);
    if(pv->el==NULL)return NULL;
    int i;
    for ( i = 0; i < n; i++)
    {
        pv->el[i]=0;
    }
    return pv;
}

vector* newOnesVector(int n){
    vector* pv=(vector*)malloc(sizeof(vector));
    if(pv==NULL)return NULL;
    pv->n=n;
    pv->el=(double*)calloc(n,sizeof(double));
    int i;
    for(i=0;i<n;i++){
        pv->el[i]=1;
    }
    return pv;
}

int copy_vector(vector* src,vector* des){
    if(src->n!=des->n)return 0;
    int i;
    for(i=0;i<src->n;i++){
        des->el[i]=src->el[i];
    }
    return 1;
}

int average_vector(vector* src1,vector* src2,vector* des){
    if(src1->n!=src2->n||des->n!=src2->n)return 0;
    int i;
    for(i=0;i<src1->n;i++){
        des->el[i]=(src1->el[i]+src2->el[i])/2;
    }
    return 1;
}

void print_vector(vector* b){
    int i;
    printf("(");
    for ( i = 0; i < b->n; i++)
    {
        printf("%lf",b->el[i]);
        if(i!=b->n-1)printf(",");
    }
    printf(")");
    printf("\n");
}
void fprint_vector(vector* b,FILE* fdata){
    int i;
    for(i=0;i<b->n;i++){
        fprintf(fdata,"%lf",b->el[i]);
        if(i!=b->n-1)fprintf(fdata,"\t");
        else fprintf(fdata,"\n");
    }
}
void freeVector(vector* pv){
    free(pv->el);
    free(pv);
}

double maxEl(vector* pv){
    int i;
    double max=pv->el[0];
    if(pv->n==1)return max;
    for(i=1;i<pv->n;i++){
        if(pv->el[i]>max){
            max=pv->el[i];
        }
    }
    return max;
}

double infNorm(vector* pv){
    int i;
    double max=fabs(pv->el[0]);
    if(pv->n==1)return max;
    for(i=1;i<pv->n;i++){
        if(fabs(pv->el[i])>max)max=fabs(pv->el[i]);
    }
    return max;
}
double distance(vector* pv1,vector* pv2){
    if(pv1==NULL||pv2==NULL)return 999999;
    if(pv1->n!=pv2->n)return 999999;
    vector* pd=createNewVector(pv1->n);
    int i;
    for(i=0;i<pv1->n;i++){
        pd->el[i]=pv1->el[i]-pv2->el[i];
    }
    double d=infNorm(pd);
    freeVector(pd);
    return d;
}

int scalarMul(double a,vector* src,vector* des){
    if(src==NULL||src->el==NULL||des==NULL||des->el==NULL||src->n!=des->n)return 0;
    int i;
    for(i=0;i<src->n;i++){
        des->el[i]=a*src->el[i];
    }
    return 1;
}

int vectorAdd(vector* add1,vector* add2, vector* des){
    if(add1->n<=0||add1->n!=add2->n||des->n!=add2->n)return 0;
    int i;
    for(i=0;i<des->n;i++){
        des->el[i]=add1->el[i]+add2->el[i];
    }
    return 1;
}
int vectorMinus(vector* minuend,vector* subtrahend,vector* des){
    if(minuend->n<=0||minuend->n!=subtrahend->n||des->n!=minuend->n)return 0;
    int i;
    for(i=0;i<des->n;i++){
        des->el[i]=minuend->el[i]-subtrahend->el[i];
    }
    return 1;
}
double scalarDot(vector* pa,vector* pb){
    if(pa->n<=0||pa->n!=pb->n){
        printf("scalarDot wrong\n");
        return 999999;
    }
    double s=0;
    int i;
    for(i=0;i<pa->n;i++){
        s+=pa->el[i]*pb->el[i];
    }
    return s;
}
double vector2Norm(vector* pv){
    if(pv->n<=0)return -1;
    double s=scalarDot(pv,pv);
    if(s<0)return -1;
    s=sqrt(s);
    return s;
}
int normalize(vector* pv){
    if(pv->n<1)return 0;
    int i;
    double s=vector2Norm(pv);
    if(s<EPS){
        printf("Elements too small\n");
        return 0;
    }
    for(i=0;i<pv->n;i++){
        pv->el[i]/=s;
    }
    return 1;
}
#endif