#include <iostream>
#include "field.h"
#include "schemes.h"

/* Driver method to advance in time the hyperbolic part of the 2-phase formula-
 * tion. Overall structure:
 *      1) compute fluxes at cell centers, E_{i}^{n};
 *      2) extrapolate to cell faces, E_{i \pm 1/2}^{n};
 *      3) advance in time:
 *
 *      Q_{i}^{n+1} = Q_{i}^{n} - \frac{dt}{dx} (E_{i+1/2}^{n} - E_{i-1/2}^{n})
 *                  + dt H_{i}^{n} \Delta
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure with EoS parameters;
 *  struct phase phase2 : phase1 structure with EoS parameters;
 *  double       dt     : time increment;
 *  double       dx     : grid spacing;
 *  int          nt     : total number of cells, counting halos (n+2*nhc-1);
 *  double**&    Qn     : reference to array of conservative variables at time
 *                        step n;
 *  double**&    Qnp1   : reference to array of conservative variables at time
 *                        step n+1;
 *  double**&    En     : reference to array of fluxes based on Qn at time step n;
 *  double**&    Enf    : reference to array of fluxes extrapolated to cell faces 
 *                        based on En at time step n;
 *
 * TODO:
 * -----
 *  - check how to use it as local variable without having to allocate memory for
 *  it everytime;
 */
void advanceTimeHyperbolicGudunov(
    phase phase1, phase phase2, 
    double dt, double dx, 
    int nt, 
    double**& Qn, double**& Qnp1,
    double**& En, double***& Enf
)
{
    //std::cout << "advanceTimeHyperbolicGudunov" << std::endl;

    // Compute fluxes at cell centers;
    computeFluxCellCenters(phase1, phase2, nt, Qn, En);
}
