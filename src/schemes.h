#ifndef SCHEMES
#define SCHEMES

#include "field.h"

void advanceTimeHyperbolicGudunov(
    phase phase1, phase phase2, 
    double dt, double dx,
    int nt, 
    double**& Qn, double**& Qnp1,
    double**& En, double***& Enpmh
);

#endif
