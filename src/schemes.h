#ifndef SCHEMES
#define SCHEMES

#include "field.h"

void advanceTimeHyperbolicGudunov(
    phase phase1, phase phase2, 
    double dt, double dx,
    int n, int nhc, 
    double**& Qn, double**& Qnp1,
    double**& Snf, double**& En, double**& Enpmh
);

#endif
