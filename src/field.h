#ifndef FIELD
#define FIELD

struct phase {
    double rho;
    double p;
    double u;
    double alpha;
    double gamma; // Parameter 1 for Stiffened Gas EoS;
    double pi;    // Parameter 2 for Stiffened Gas EoS;
    double epsilon;
};

void initializeWaterAirShockTube(
    phase phase1, phase phase2, double xI, 
    int n, int nhc, double*& xcs,
    double**& Q0
);

void phaseToConservative(
    phase phase1, phase phase2, int n, int nhc,
    double**& Qn
);
//void conservativeToPhase(phase& phase1, phase& phase2, double**& Qn);

void computeBCs(int n, int nhc, double**& Qn);
void computeFluxCellCenters(phase phase1, phase phase2, int nt, double**& Qn, double**& En);

#endif
