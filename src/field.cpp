#include <iostream>
#include <cmath>
#include "field.h"
#include "allocate.h"
#include "log.h"

/* Method to initialize the water-air shock tube simulation. The inteface is de-
 * fined by xI, before it mostly phase 1 and after it mostly phase 2, except for
 * a small epsilon.
 *
 * Details on the floating point operations to replace boolean logic:
 *
 * a  = x < x0 = -\frac{(x-x0) - |x-x0|}{2 |x-x0|} \\ boolean array of x < x0
 * Qn = a*phase1 + (1-a)*phase2
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : structure phase1 defined in "field.h";
 *  struct phase phase2 : structure phase2 defined in "field.h";
 *  double       xI     : location of the interface;
 *  int          n      : number of nodes;
 *  int          nhc    : number of halo cells;
 *  double*&     xcs    : reference to array of cell centers padded with halo 
 *                        cells;
 *  double**&    Q0     : reference to array of conservative variables;
 *
 * TODO:
 * -----
 *  - test if one should do the "boolean" operation already inside the for loop
 *  to assign the values for \alpha_1;
 */
void initializeWaterAirShockTube(
    phase phase1, phase phase2, double xI, 
    int n, int nhc, double*& xcs, 
    double **&Q0)
{
    std::cout << "initializeAirShockTube" << std::endl;

    // Compute boolean array for location of phase 1, which is before xI:
    double* bp1;
    allocate1d(n+2*nhc-1, bp1);
    for (int i=0; i<(n+2*nhc-1); i++) {
        bp1[i] = -((xcs[i]-xI)-fabs(xcs[i]-xI))/(2*fabs(xcs[i]-xI));
    }
    
    // Compute alpha1 array in Q0[:][0], interior cells then fill halos:
    for (int i=nhc; i<(n+nhc-1); i++) {
        Q0[i][0] = bp1[i]*phase1.alpha + (1-bp1[i])*phase1.epsilon;
    }
    for (int i=0; i<nhc; i++) {
        Q0[nhc-(i+1)  ][0] = Q0[nhc+i][0];
        Q0[nhc+(i+n-1)][0] = Q0[nhc+(-i+n-2)][0];
    }
    phaseToConservative(phase1, phase2, n, nhc, Q0);
    //print2DArray(n+2*nhc-1, 1, Q0);
    deallocate1d(bp1);
}

/* Method to convert from struct phase to conservative vector.
 *
 * Definition of vector:
 *
 * Q = (\alpha_1, 
 *      \alpha_1\rho_1, \alpha_1\rho_1\u_1, \alpha_1_\rho_1 e_t^1,
 *      \alpha_2\rho_2, \alpha_2\rho_2\u_2, \alpha_2_\rho_2 e_t^2,
 * )
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : structure phase1 defined in "field.h";
 *  struct phase phase2 : structure phase2 defined in "field.h";
 *  int          n      : number of nodes;
 *  int          nhc    : number of halo cells;
 *  double**&    Qn     : reference to array of conservative variables;
 */
void phaseToConservative(phase phase1, phase phase2, int n, int nhc, double**& Qn)
{
    std::cout << "phaseToConservative" << std::endl;
    // Compute interior cell values:
    double e1, e2;
    for (int i=nhc; i<(n+nhc-1); i++) {
        Qn[i][1] = Qn[i][0]*phase1.rho;
        Qn[i][2] = Qn[i][0]*phase1.rho*phase1.u;
        e1       = phase1.p/(phase1.gamma - 1.0) \
                 + phase1.gamma*phase1.pi/(phase1.gamma-1);
        Qn[i][3] = Qn[i][0]*phase1.rho*(e1 + 0.5*phase1.u*phase1.u);
        Qn[i][4] = (1.0 - Qn[i][0])*phase2.rho;
        Qn[i][5] = (1.0 - Qn[i][0])*phase2.rho*phase2.u;
        e2       = phase2.p/(phase2.gamma - 1.0) \
                 + phase2.gamma*phase2.pi/(phase2.gamma-1);
        Qn[i][6] = (1.0 - Qn[i][0])*phase2.rho*(e2 + 0.5*phase2.u*phase2.u);
    }

    // Compute BCs:
    // No-slip for velocity (left wall for sure, hard coded no-gradient for right
    // wall as well);
    // No-gradient for scalars;
    computeBCs(n, nhc, Qn);
}

/* Method to fill up the halo cells. The definition of Q follows from above.
 * Here, the boundary conditions are defined as:
 *      1) no-gradient for scalars (volume fraction, pressure, temperature);
 *      2) no-slip/impermeability for velocity;
 *
 * Parameters:
 * -----------
 *  int       n   : number of nodes;
 *  int       nhc : number of halo cells;
 *  double**& Qn  : reference to array of conservative variables;
 *
 * TODO:
 * -----
 *  - improve for custom boundary condition types;
 */
void computeBCs(int n,  int nhc, double**& Qn)
{
    std::cout << "computeBCs" << std::endl;

    // Compute BCs:
    // No-slip for velocity;
    // No-gradient for scalars;
    for (int i=0; i<nhc; i++) {
        Qn[nhc-(i+1)  ][0] =  Qn[nhc+i][0];
        Qn[nhc+(i+n-1)][0] =  Qn[nhc+(-i+n-2)][0];
        Qn[nhc-(i+1)  ][1] =  Qn[nhc+i][1];
        Qn[nhc+(i+n-1)][1] =  Qn[nhc+(-i+n-2)][1];
        Qn[nhc-(i+1)  ][2] = -Qn[nhc+i][2];
        Qn[nhc+(i+n-1)][2] = -Qn[nhc+(-i+n-2)][2];
        Qn[nhc-(i+1)  ][3] =  Qn[nhc+i][3];
        Qn[nhc+(i+n-1)][3] =  Qn[nhc+(-i+n-2)][3];
        Qn[nhc-(i+1)  ][4] =  Qn[nhc+i][4];
        Qn[nhc+(i+n-1)][4] =  Qn[nhc+(-i+n-2)][4];
        // No-slip for right wall:
        //Qn[nhc-(i+1)  ][5] = -Qn[nhc+i][5];
        //Qn[nhc+(i+n-1)][5] = -Qn[nhc+(-i+n-2)][5];
        // No-gradient:
        Qn[nhc-(i+1)  ][5] =  Qn[nhc+i][5];
        Qn[nhc+(i+n-1)][5] =  Qn[nhc+(-i+n-2)][5];
        Qn[nhc-(i+1)  ][6] =  Qn[nhc+i][6];
        Qn[nhc+(i+n-1)][6] =  Qn[nhc+(-i+n-2)][6];
    }
}

/* Method to compute fluxes at cell centers. All algebraic operations for phasic
 * pressure consider an stiffened gas EoS.
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : structure phase1 defined in "field.h";
 *  struct phase phase2 : structure phase2 defined in "field.h";
 *  int          nt     : total number of cells padded with halos;
 *  double**&    Qn     : reference to array of conservative variables;
 *  double**&    En     : reference to array of fluxes based on Qn;
 *
 * TODO:
 * -----
 *  - implement arbitraty EoS to compute \alpha_k p_k;
 */
void computeFluxCellCenters(phase phase1, phase phase2, int nt, double **&Qn, double **&En)
{
    //std::cout << "computeFluxCellCenters" << std::endl;
    
    // Loop through cells. Be aware that q_0 is in fact \alpha_1.
    double u_k, a_kp_k;
    for (int i=0; i<nt; i++) {
        // Phase 1:
        En[i][0] = Qn[i][2];
        u_k = Qn[i][2]/Qn[i][1];
        a_kp_k = (phase1.gamma-1.0)*(Qn[i][3]-0.5*Qn[i][2]*u_k) \
               - Qn[i][0]*phase1.gamma*phase1.pi;
        En[i][1] = Qn[i][2]*u_k + a_kp_k;
        En[i][2] = Qn[i][3]*u_k + a_kp_k*u_k;
        // Phase 2:
        En[i][3] = Qn[i][5];
        u_k = Qn[i][5]/Qn[i][4];
        a_kp_k = (phase2.gamma-1.0)*(Qn[i][6]-0.5*Qn[i][5]*u_k) \
               - (1.0-Qn[i][0])*phase2.gamma*phase2.pi;
        En[i][4] = Qn[i][5]*u_k + a_kp_k;
        En[i][5] = Qn[i][6]*u_k + a_kp_k*u_k;
    }
}

/* Method to update Qn based on Qn+1.
 *
 * Parameters:
 * -----------
 *  int       nt   : total number of cells padded with halos;
 *  double**& Qn   : reference to array of conservative variables;
 *  double**& Qnp1 : reference to array of fluxes based on Qn;
 */
void updateQn(int nt, double**& Qn, double**& Qnp1)
{
    //std::cout << "updateQn" << std::endl;

    for (int i=0; i<nt; i++) {
        for (int j=0; j<7; j++) {
            Qn[i][j] = Qnp1[i][j];
        }
    }
}
