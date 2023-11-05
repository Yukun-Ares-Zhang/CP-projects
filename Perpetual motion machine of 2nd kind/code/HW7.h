#ifndef __HW7_H__
#define __HW7_H__

#include "vector.h"
#include<time.h>
#include<string.h>

typedef struct transition{
    int situation;
    char absorb;
    vector* r;
    vector* n;
}transition;

double min(double a,double b){
    return a<b?a:b;
}

void print_transition(transition* pt){
    printf("%d\n",pt->situation);
    if(pt->situation==1){
        printf("absorb=%c\n",pt->absorb);
        printf("r=");
        print_vector(pt->r);
    }
    if(pt->situation==2){
        printf("r=");
        print_vector(pt->r);
        printf("n=");
        print_vector(pt->n);
    }
}

#endif