#include <iostream>
#include <cmath>
#include "relaxation.h"
#include "../field.h"

/* Method to relax each phase velocity to its equilibrium value, ie, the mixture 
 * velocity and readjust the phasic energy due to the relaxation of the velocity.
 *
 * Parameters:
 * -----------
 *  int          n       : total number of nodes; 
 *  int          nhc     : number of halo cells; 
 *  double**&    Qn      : reference to array of conservative variables at time
 *                         step n;
 *  double**&    Qh      : reference to array of conservative variables after 
 *                         hyperbolic integration;
 *  double**&    Qhu     : reference to array of conservative variables after 
 *                         hyperbolic integration and velocity relaxation;
 * 
 * TODO:
 * -----
 *  - understand the update in internal energy for the V_i^n term, is it the one
 *  before or after the hyperbolic integration?
 *  - check sign in the internal energy update: plus or minus for e_1_hu, etc;
 */
void relaxationVelocity(int n, int nhc, double **&Qn, double **&Qh, double **&Qhu)
{
    //std::cout << "relaxationVelocity" << std::endl;

    double u_I_n, u_1_h, u_2_h, u_I_hu;
    double e_1_h, e_2_h, e_1_hu, e_2_hu;
    for (int i=nhc; i<(n+nhc-1); i++) {
        // Compute equilibrium velocity after relaxation:
        u_I_n  = (Qn[i][2]+Qn[i][5])/(Qn[i][1]+Qn[i][4]);
        u_1_h  =  Qh[i][2]/Qh[i][1];
        u_2_h  =  Qh[i][5]/Qh[i][4];
        u_I_hu = (Qh[i][2]+Qh[i][5])/(Qh[i][1]+Qh[i][4]);

        // Compute readjusted internal energy due to velocity relaxation:
        e_1_h  =  Qh[i][3]/Qh[i][1] - 0.5*u_1_h*u_1_h;
        e_2_h  =  Qh[i][6]/Qh[i][4] - 0.5*u_2_h*u_2_h;
        e_1_hu =  e_1_h + 0.5*(u_I_hu-u_1_h)*(u_I_n-u_1_h);
        e_2_hu =  e_2_h - 0.5*(u_I_hu-u_2_h)*(u_I_n-u_2_h);
        //e_1_hu =  e_1_h - 0.5*(u_I_hu-u_1_h)*(u_I_n-u_1_h);
        //e_2_hu =  e_2_h + 0.5*(u_I_hu-u_2_h)*(u_I_n-u_2_h);

        // Update conservative variables array:
        Qhu[i][0] = Qh[i][0];
        Qhu[i][1] = Qh[i][1];
        Qhu[i][2] = Qh[i][1]*u_I_hu;
        Qhu[i][3] = Qh[i][1]*(e_1_hu + 0.5*u_I_hu*u_I_hu);
        Qhu[i][4] = Qh[i][4];
        Qhu[i][5] = Qh[i][4]*u_I_hu;
        Qhu[i][6] = Qh[i][4]*(e_2_hu + 0.5*u_I_hu*u_I_hu);
    }
    
    computeBCs(n, nhc, Qhu);
}

double functionAlphaStar(
    phase phase1, phase phase2, 
    double*& Qhu, 
    double a_s, double*& e_k_s
);

/* Method to relax each phase pressure to its equilibrium value, .
 *
 * Bisection method to solve for:
 *
 *      f(\alpha_1*) = p_1*(\rho_1*, e_1*) - p_2*(\rho_2*, e_2*) = 0
 * 
 * where:
 *      p_k* = (\gamma_k -1)\rho_k* e_1* - \gamma_k \pi_k
 *      \rho_k* = (\alpha_k^0 \rho_k^0)/\alpha_k* (?)
 *      e_k* = e_k^0 \pm \frac{\bar{p_I}}{\alpha_k^0\rho_k^0}(\alpha_k*-\alpha_k^0)
 *      \bar{p_I} = 0.5*(p_I* + p_I^0)
 *  `
 * Parameters:
 * -----------
 *  int          n       : total number of nodes; 
 *  int          nhc     : number of halo cells; 
 *  double**&    Qn      : reference to array of conservative variables at time
 *                         step n;
 *  double**&    Qh      : reference to array of conservative variables after 
 *                         hyperbolic integration;
 *  double**&    Qhu     : reference to array of conservative variables after 
 *                         hyperbolic integration and velocity relaxation;
 *  double**&    Qhup    : reference to array of conservative variables after 
 *                         hyperbolic integration, velocity and pressure relaxa-
 *                         tion;
 */
void relaxationPressure(
    phase phase1, phase phase2,
    int n, int nhc, 
    double **&Qn, double **&Qh, double **&Qhu, double **&Qhup
)
{
    std::cout << "relaxationPressure" << std::endl;

    // Variables for bisection method:
    double err, max_err;
    int it, max_it;
    max_err = 1.0e-8;
    max_it  = 50;
    double a, b, c;
    double fa, fb, fc;

    // Variables after relaxation:
    double u_k_s;
    double *e_k_s = new double[2];

    for (int i=nhc; i<(n+nhc-1); i++) {
        it  = 1;
        a   = 0.0+1.0e-8;
        b   = 1.0-1.0e-8;
        err = fabs(a-b);
        while (err > max_err && it < max_it) {
            c = 0.5*(a+b); // \alpha_1*
            fc = functionAlphaStar(phase1, phase2, Qhu[i], c, e_k_s); 
            fa = functionAlphaStar(phase1, phase2, Qhu[i], a, e_k_s); 
            fb = functionAlphaStar(phase1, phase2, Qhu[i], b, e_k_s); 
            if (fa * fb >= 0.0) {
                std::cout << "No root or multiple roots at " << i-nhc+1;
                c = Qhu[i][0];
                break;
            } else if (fc * fa < 0.0) {
                b = c;
                err = fabs(a-b);
                it += 1;
            } else if (fc * fb < 0.0) {
                a = c;
                err = fabs(a-b);
                it += 1;
            } else {
                std::cout << "Unknown error." << std::endl;
            }
        } // c is the \alpha_1*
        
        fc = functionAlphaStar(phase1, phase2, Qhu[i], c, e_k_s); 
        // Update conservative variables array:
        Qhup[i][0] = c;
        Qhup[i][1] = Qhu[i][1];
        Qhup[i][2] = Qhu[i][2];
        u_k_s = Qhu[i][2]/Qhu[i][1];
        Qhup[i][3] = Qhu[i][1]*(e_k_s[0] + 0.5*u_k_s*u_k_s);
        Qhup[i][4] = Qhu[i][4];
        Qhup[i][5] = Qhu[i][5];
        u_k_s = Qhu[i][5]/Qhu[i][4];
        Qhup[i][6] = Qhu[i][4]*(e_k_s[1] + 0.5*u_k_s*u_k_s);
    }

    computeBCs(n, nhc, Qhup);
    
    delete[] e_k_s;
}

/* Method to evaluate function f(\alpha_1*).
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure with EoS parameters;
 *  struct phase phase2 : phase1 structure with EoS parameters;
 *  double*&     Qhu    : reference to array of conservative variables after 
 *                        hyperbolic integration and velocity relaxation at cell i;
 *  double       a_s    : \alpha_1* guessed;
 *  double*&     e_k_s  : reference of array pointer for internal energies after
 *                        pressure relaxation;
 */
double functionAlphaStar(
    phase phase1, phase phase2, 
    double*& Qhu, 
    double a_s, double*& e_k_s
)
{
    //std::cout << "functionAlphaStar" << std::endl;
    
    // Initial state, ie, after hyperbolic integration and velocity relaxation:
    double rho_1_0, u_1_0, e_1_0, a_1_0, p_1_0;
    double rho_2_0, u_2_0, e_2_0, a_2_0, p_2_0;
    double p_I_0;

    rho_1_0 = Qhu[1]/Qhu[0];
    u_1_0   = Qhu[2]/Qhu[1];
    e_1_0   = Qhu[3]/Qhu[1] - 0.5*u_1_0*u_1_0;
    p_1_0   = (phase1.gamma - 1.0)*rho_1_0*e_1_0 + phase1.gamma*phase1.pi;
    rho_2_0 = Qhu[4]/(1.0-Qhu[0]);
    u_2_0   = Qhu[5]/Qhu[4];
    e_2_0   = Qhu[6]/Qhu[4] - 0.5*u_2_0*u_2_0;
    p_2_0   = (phase2.gamma - 1.0)*rho_2_0*e_2_0 + phase2.gamma*phase2.pi;
    p_I_0   = Qhu[0]*p_1_0 + (1.0-Qhu[0])*p_2_0;

    double var;
    double rho_1_s, e_1_s, p_1_s;
    double rho_2_s, e_2_s, p_2_s;
    double p_I_s; // approximating as p_1_s, since at equilibrium p_I = p_1 = p_2;
    rho_1_s = Qhu[1]/a_s;
    rho_2_s = Qhu[4]/(1.0-a_s);
    // Estimate of p_1_s:
    var     = (phase1.gamma - 1.0)*rho_1_s*(a_s - Qhu[0])/(2.0*Qhu[1]);
    p_1_s   = ((phase1.gamma - 1.0)*rho_1_s*e_1_0 + phase1.gamma*phase1.pi \
              + var*p_I_0)/(1.0-var);
    p_I_s   = 0.5*(p_I_0 + p_1_s);
    // Compute e_k* and correct p_k*
    e_1_s   = e_1_0 + p_I_s*(a_s - Qhu[0])/Qhu[1];
    e_2_s   = e_2_0 - p_I_s*(a_s - Qhu[0])/Qhu[4];
    // Update phasic pressures:
    p_1_s   = (phase1.gamma - 1.0)*rho_1_s*e_1_s + phase1.gamma*phase1.pi;
    p_2_s   = (phase2.gamma - 1.0)*rho_2_s*e_2_s + phase2.gamma*phase2.pi;

    e_k_s[0] = e_1_s;
    e_k_s[1] = e_2_s;

    return p_1_s - p_2_s;
}
