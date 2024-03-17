#include <iostream>
#include <cmath>
#include "field.h"
#include "schemes.h"
#include "thermo.h"

/* Method to maximum wave speed at cell face. 
 *
 * Definition: S_{i+1/2} = \max \{
 *      |\lambda_i^+|, |\lambda_i^-|, |\lambda_{i+1}^+|, |\lambda_{i+1}^-|
 *  \}
 *
 * Note: I'm comparing with \lambda_{i(+1)}^0 and \lambda_{i(+1)}^{I}. I'm not
 * so sure if this is the correct approach.
 *
 * Parameters:
 * -----------
 * double*& u_I        : velocity of the interface at cell i and i+1, "lambda_I_0";
 * double*& u_l        : velocity of phase l at cell i and i+1, "lambda_l_0";
 * double*& lambda_l_p : "positive" wave speed of phase l at cell i and i+1;
 * double*& lambda_l_m : "negative" wave speed of phase l at cell i and i+1;
 *
 * Returns:
 * --------
 *  double max_wave_speed : maximum in absolute value of all wave speeds;
 *
 * TODO:
 * -----
 *  - improve this maximum operation;
 */
double maxWaveSpeedGudunovCellFace(
    double*& u_I, double*& u_l, double*& lambda_l_p, double*& lambda_l_m
)
{
    //std::cout << "maxWaveSpeedGudunovCellFace" << std::endl;

    double max_wave_speed = 0.0;
    max_wave_speed = fmax(max_wave_speed, fabs(u_I[0]));
    max_wave_speed = fmax(max_wave_speed, fabs(u_I[1]));
    max_wave_speed = fmax(max_wave_speed, fabs(u_l[0]));
    max_wave_speed = fmax(max_wave_speed, fabs(u_l[1]));
    max_wave_speed = fmax(max_wave_speed, fabs(lambda_l_p[0]));
    max_wave_speed = fmax(max_wave_speed, fabs(lambda_l_p[1]));
    max_wave_speed = fmax(max_wave_speed, fabs(lambda_l_m[0]));
    max_wave_speed = fmax(max_wave_speed, fabs(lambda_l_m[1]));
    return max_wave_speed;
}

/* Method to compute Gudunov wave speeds at cell faces. 
 *
 * Parameters:
 * -----------
 *  struct phase phase1           : phase1 structure with EoS parameters;
 *  struct phase phase2           : phase1 structure with EoS parameters;
 *  double*&     Qi               : reference to array of conservative variables 
 *                                  at cell i;
 *  double*&     Qip1             : reference to array of conservative variables 
 *                                  at cell i+1;
 *  double*&     local_wave_speed : reference to array of local wave speed for
 *                                  phases 1 and 2 at the cell face;
 */
void computeWaveSpeedGudunovCellFace(
    phase phase1, phase phase2,
    double*& Qi, double*& Qip1,
    double*& local_wave_speed
)
{
    //std::cout << "computeWaveSpeedGudunov" << std::endl;

    // Compute speed of sound;
    double *c_l = new double[2];

    c_l[0] = stiffenedGasSoundSpeed(
        phase1, Qi[0], Qi[1], Qi[2], Qi[3]
    );
    c_l[1] = stiffenedGasSoundSpeed(
        phase1, Qip1[0], Qip1[1], Qip1[2], Qip1[3]
    );
    
    // Compute wave speed for phase 1:
    double *lambda_l_p = new double[2], \
           *lambda_l_m = new double[2], \
           *u_l        = new double[2], \
           *u_I        = new double[2];
    
    u_l[0] = Qi[2]/Qi[1];
    u_l[1] = Qip1[2]/Qip1[1];

    u_I[0] = (Qi[2]+Qi[5])/(Qi[1]+Qi[4]);
    u_I[1] = (Qip1[2]+Qip1[5])/(Qip1[1]+Qip1[4]);

    lambda_l_p[0] = u_l[0] + c_l[0];
    lambda_l_p[1] = u_l[1] + c_l[1];
    lambda_l_m[0] = u_l[0] - c_l[0];
    lambda_l_m[1] = u_l[1] - c_l[1];

    local_wave_speed[0] = maxWaveSpeedGudunovCellFace(
        u_I, u_l, lambda_l_p, lambda_l_m
    );

    // Compute wave speed for phase 2:
    c_l[0] = stiffenedGasSoundSpeed(
        phase2, (1.0-Qi[0]), Qi[4], Qi[5], Qi[6]
    );
    c_l[1] = stiffenedGasSoundSpeed(
        phase2, (1.0-Qip1[0]), Qip1[4], Qip1[5], Qip1[6]
    );

    u_l[0] = Qi[5]/Qi[4];
    u_l[1] = Qip1[5]/Qip1[4];

    lambda_l_p[0] = u_l[0] + c_l[0];
    lambda_l_p[1] = u_l[1] + c_l[1];
    lambda_l_m[0] = u_l[0] - c_l[0];
    lambda_l_m[1] = u_l[1] - c_l[1];

    local_wave_speed[1] = maxWaveSpeedGudunovCellFace(
        u_I, u_l, lambda_l_p, lambda_l_m
    );

    delete[] c_l;
    delete[] u_l;
    delete[] u_I;
    delete[] lambda_l_p;
    delete[] lambda_l_m;
}

/* Method to extrapolate fluxes to cell faces considering wave speeds. 
 *
 * Define Snf as maximum wave speeds in the form:
 *      Snf[face i][phase 1 or 2 / compared with interface value];
 *
 * And Enf in the form:
 *      Enf[face i][flux of property j];
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure with EoS parameters;
 *  struct phase phase2 : phase1 structure with EoS parameters;
 *  int          n      : number of nodes;
 *  int          nhc    : number of halo cells;
 *  double**&    Qn     : reference to array of conservative variables at time
 *                        step n;
 *  double**&    Snf    : reference to array of wave speeds at cell faces;
 *  double**&    En     : reference to array of fluxes based on Qn at time step n;
 *  double**&    Enf    : reference to array of fluxes extrapolated to cell faces 
 *                        based on En at time step n;
 */
void extrapolateGudunovFluxCellFaces(
    phase phase1, phase phase2, 
    int n, int nhc, 
    double**& Qn, 
    double**& Snf, double**& En, double**& Enf
)
{
    //std::cout << "extrapolateGudunovFluxCellFaces" << std::endl;
    
    // Compute wave speeds at cell faces:
    // HERE - check wave speeds required only in faces inside the domain and its
    // boundary
    double* local_wave_speed = new double[2];
    for (int i=nhc; i<(n+nhc); i++) {
        // Will compare \lambda_l with \lambda_I, check if required.
        computeWaveSpeedGudunovCellFace(
            phase1, phase2, 
            Qn[i], Qn[i+1],
            local_wave_speed
        );
        Snf[i][0] = local_wave_speed[0];
        Snf[i][1] = local_wave_speed[1];
    }

    // Extrapolate fluxes at cell faces based on the Gudunov approach:
    // HERE - check fluxes required only in faces inside the domain and its
    // boundary
    for (int i=n; i<(n+nhc); i++) {
        for (int j=0; j<3; j++) {
            Enf[i][j] = 0.5*(En[i][j] + En[i+1][j] \
                            -Snf[i][0]*(Qn[i+1][j+1] - Qn[i][j+1]));
        }
        for (int j=3; j<6; j++) {
            Enf[i][j] = 0.5*(En[i][j] + En[i+1][j] \
                            -Snf[i][1]*(Qn[i+1][j+1] - Qn[i][j+1]));
        }
    }

    delete[] local_wave_speed;
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
    int n, int nhc, 
    double**& Qn, double**& Qnp1,
    double**& Snf, double**& En, double**& Enf
)
{
    //std::cout << "advanceTimeHyperbolicGudunov" << std::endl;

    // Compute fluxes at cell centers:
    computeFluxCellCenters(phase1, phase2, n+2*nhc-1, Qn, En);

    // Extrapolate fluxes at cell faces:
    extrapolateGudunovFluxCellFaces(phase1, phase2, n, nhc, Qn, Snf, En, Enf);

    // Advance in time:
    for (int i=0; i<n; i++) {
    }
}
