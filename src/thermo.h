#ifndef THERMO
#define THERMO

#include "field.h"

double stiffenedGasSoundSpeed(
    phase phase, double a, double a_rho, double a_rho_u, double a_rho_et
);
double computePressureInterface(phase phase1, phase phase2, double*& Qn);
double computeVelocityInterface(phase phase1, phase phase2, double*& Qn);

#endif
