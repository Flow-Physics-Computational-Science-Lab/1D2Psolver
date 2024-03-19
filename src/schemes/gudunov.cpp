#include <iostream>
#include <cmath>
#include <stdexcept>
#include "../field.h"
#include "gudunov.h"
#include "../thermo.h"
#include "../runtimeparameters.h"

double maxWaveSpeedGudunovCellFace(
    double*& u_I, double*& u_l, double*& lambda_l_p, double*& lambda_l_m
);

void computeWaveSpeedGudunovCellFace(
    phase phase1, phase phase2,
    double*& Qi, double*& Qip1,
    double*& local_wave_speed
);

void extrapolateGudunovFluxCellFaces(
    phase phase1, phase phase2, 
    RunTimeParameters sim_par,
    int n, int nhc, 
    double**& Qn, 
    double**& Snf, double**& En, double**& Enf
);

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
 *  struct phase phase1  : phase1 structure with EoS parameters;
 *  struct phase phase2  : phase1 structure with EoS parameters;
 *  struct RunTimeParameters
 *               sim_par : sim_par structure with dt and dx parameters;
 *  int          nt      : total number of cells, counting halos (n+2*nhc-1); double**&    Qn      : reference to array of conservative variables at time
 *                         step n;
 *  double**&    Qnp1    : reference to array of conservative variables at time
 *                         step n+1;
 *  double**&    Snf     : reference to array of wave speeds at cell faces at time 
 *                         step n;
 *  double**&    En      : reference to array of fluxes based on Qn at time step n;
 *  double**&    Enf     : reference to array of fluxes extrapolated to cell faces 
 *                         based on En at time step n;
 *
 * TODO:
 * -----
 *  - check how to use it as local variable without having to allocate memory for
 *  it everytime;
 *  - improve exception for CFL check;
 */
void advanceTimeHyperbolicGudunov(
    phase phase1, phase phase2, 
    RunTimeParameters sim_par,
    int n, int nhc, 
    double**& Qn, double**& Qnp1,
    double**& Snf, double**& En, double**& Enf
)
{
    //std::cout << "advanceTimeHyperbolicGudunov" << std::endl;

    // Compute fluxes at cell centers:
    computeFluxCellCenters(phase1, phase2, n+2*nhc-1, Qn, En);

    // Extrapolate fluxes at cell faces:
    try {
        extrapolateGudunovFluxCellFaces(phase1, phase2, sim_par, n, nhc, Qn, Snf, En, Enf);
    } catch (const std::runtime_error) {
        throw std::runtime_error("Insufficient dt for prescribed CFL.");
    }

    // Advance in time:
    double p_I, u_I;
    double dt = sim_par.dt, dx = sim_par.dx;
    for (int i=nhc; i<(n+nhc-1); i++) {
        p_I = computePressureInterface(phase1, phase2, Qn[i]);
        u_I = computeVelocityInterface(phase1, phase2, Qn[i]);
        // Volume fraction of phase 1:
        Qnp1[i][0] = Qn[i][0] \
            - (dt/(2.0*dx))*((Qn[i][2]/Qn[i][1])*(Qn[i+1][0]-Qn[i-1][0])
                             -Snf[i-nhc+1][0]*(Qn[i+1][0]-Qn[i][0])
                             +Snf[i-nhc][0]*(Qn[i][0]-Qn[i-1][0]));
        // Phase 1:
        Qnp1[i][1] = Qn[i][1] - (dt/dx)*(Enf[i-nhc+1][0]-Enf[i-nhc][0]);
        Qnp1[i][2] = Qn[i][2] - (dt/dx)*(Enf[i-nhc+1][1]-Enf[i-nhc][1]) \
                   + (dt/(2.0*dx))*p_I*(Qn[i+1][0]-Qn[i-1][0]);
        Qnp1[i][3] = Qn[i][3] - (dt/dx)*(Enf[i-nhc+1][2]-Enf[i-nhc][2]) \
                   + (dt/(2.0*dx))*p_I*u_I*(Qn[i+1][0]-Qn[i-1][0]);
        // Phase 2:
        Qnp1[i][4] = Qn[i][4] - (dt/dx)*(Enf[i-nhc+1][3]-Enf[i-nhc][3]);
        Qnp1[i][5] = Qn[i][5] - (dt/dx)*(Enf[i-nhc+1][4]-Enf[i-nhc][4]) \
                   + (dt/(2.0*dx))*p_I*(Qn[i-1][0]-Qn[i+1][0]);
        Qnp1[i][6] = Qn[i][6] - (dt/dx)*(Enf[i-nhc+1][5]-Enf[i-nhc][5]) \
                   + (dt/(2.0*dx))*p_I*u_I*(Qn[i-1][0]-Qn[i+1][0]);
    }
    computeBCs(n, nhc, Qnp1);
}

/* Method to extrapolate fluxes to cell faces considering wave speeds. 
 *
 * Define Snf as maximum wave speeds in the form:
 *      Snf[face i][phase 1 or 2 / compared with interface value];
 *
 * And Enf in the form:
 *      Enf[face i][flux of property j];
 *
 * Also checks if dt is sufficient for the prescribed CFL.
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure with EoS parameters;
 *  struct phase phase2 : phase1 structure with EoS parameters;
 *  struct RunTimeParameters
 *               sim_par : sim_par structure with dt and dx parameters;
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
    RunTimeParameters sim_par,
    int n, int nhc, 
    double**& Qn, 
    double**& Snf, double**& En, double**& Enf
)
{
    //std::cout << "extrapolateGudunovFluxCellFaces" << std::endl;
    
    // Compute wave speeds at cell faces:
    double* local_wave_speed = new double[2];
    double overall_max = 0.0;
    double local_max;
    for (int i=(nhc-1); i<(n+nhc-1); i++) {
        // Will compare \lambda_l with \lambda_I, check if required.
        computeWaveSpeedGudunovCellFace(
            phase1, phase2, 
            Qn[i], Qn[i+1],
            local_wave_speed
        );
        Snf[i-nhc+1][0] = local_wave_speed[0];
        Snf[i-nhc+1][1] = local_wave_speed[1];
        // Check if CFL is satisfied for required dt:
        local_max = fmax(local_wave_speed[0], local_wave_speed[1]);
        overall_max = fmax(overall_max, local_max);
    }

    double local_CFL_max = overall_max*sim_par.dt/sim_par.dx;
    
    std::cout << "max(CFL) = " << local_CFL_max << std::endl;

    if (local_CFL_max > sim_par.CFL) {
        delete[] local_wave_speed;
        throw std::runtime_error("Insufficient dt for prescribed CFL.");
    }

    // Extrapolate fluxes at cell faces based on the Gudunov approach:
    for (int i=(nhc-1); i<(n+nhc-1); i++) {
        for (int j=0; j<3; j++) {
            Enf[i-nhc+1][j] = 0.5*(En[i][j] + En[i+1][j] \
                                  -Snf[i-nhc+1][0]*(Qn[i+1][j+1] - Qn[i][j+1]));
        }
        for (int j=3; j<6; j++) {
            Enf[i-nhc+1][j] = 0.5*(En[i][j] + En[i+1][j] \
                                  -Snf[i-nhc+1][1]*(Qn[i+1][j+1] - Qn[i][j+1]));
        }
    }

    delete[] local_wave_speed;
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
