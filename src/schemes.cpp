#include <iostream>
#include "field.h"
#include "schemes.h"
#include "thermo.h"

/* Method to extrapolate fluxes to cell faces considering wave speeds. 
 *
 * TODO:
 *  - fix below;
 * Define Snf as maximum wave speeds in the form:
 *      Snf[cell i][west or east faces][phase 1 or 2];
 *
 * And Enf in the form:
 *      Enf[cell i][flux of property j][west or east faces];
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure with EoS parameters;
 *  struct phase phase2 : phase1 structure with EoS parameters;
 *  int          nt     : total number of cells, counting halos (n+2*nhc-1);
 *  double**&    Qn     : reference to array of conservative variables at time
 *                        step n;
 *  double***&   Snf    : reference to array of wave speeds at cell faces;
 *  double**&    En     : reference to array of fluxes based on Qn at time step n;
 *  double***&   Enf    : reference to array of fluxes extrapolated to cell faces 
 *                        based on En at time step n;
 */
void extrapolateGudunovFluxCellFaces(
    phase phase1, phase phase2, 
    int nt, 
    double**& Qn, 
    double**& Snf, double**& En, double**& Enf
)
{
    std::cout << "extrapolateGudunovFluxCellFaces" << std::endl;
    
    double c_1, c_2;
    // Compute wave speeds at cell faces:
    Snf[0][0] = 0;
    c_1 = stiffenedGasSoundSpeed(
        phase1, Qn[0][0], Qn[0][1], Qn[0][2], Qn[0][2]
    );
    c_2 = stiffenedGasSoundSpeed(
        phase2, (1.0-Qn[0][0]), Qn[0][4], Qn[0][5], Qn[0][6]
    );
    //Snf[0][1] = 
}

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
 *  double**&    Snf    : reference to array of wave speeds at cell faces at time 
 *                        step n;
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
    double**& Snf, double**& En, double**& Enf
)
{
    //std::cout << "advanceTimeHyperbolicGudunov" << std::endl;

    // Compute fluxes at cell centers:
    computeFluxCellCenters(phase1, phase2, nt, Qn, En);

    // Extrapolate fluxes at cell faces:
    extrapolateGudunovFluxCellFaces(phase1, phase2, nt, Qn, Snf, En, Enf);
}
