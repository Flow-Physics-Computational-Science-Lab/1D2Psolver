#ifndef RELAX
#define RELAX

#include "../field.h"

void relaxationVelocity(
    int n, int nhc, 
    double**& Qn, double**& Qh, double**& Qhu
);
void relaxationPressure(
    phase phase1, phase phase2,
    int n, int nhc, 
    double**& Qn, double**& Qh, double**& Qhu, double**& Qhup
);

#endif
