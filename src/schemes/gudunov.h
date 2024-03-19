#ifndef GUDUNOV
#define GUDUNOV

#include "../field.h"
#include "../runtimeparameters.h"

void advanceTimeHyperbolicGudunov(
    phase phase1, phase phase2, 
    RunTimeParameters sim_par,
    int n, int nhc, 
    double**& Qn, double**& Qnp1,
    double**& Snf, double**& En, double**& Enpmh
);

#endif
