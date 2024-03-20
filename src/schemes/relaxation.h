#ifndef RELAX
#define RELAX

void relaxationVelocity(int n, int nhc, double**& Qn, double**& Qh, double**& Qhu);
void relaxationPressure(int nt, double**& Qn, double**& Qnp1);

#endif
