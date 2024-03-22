#include <iostream>
#include <cmath>
#include "thermo.h"
#include "field.h"

/* Method to compute speed of sound for the specified phase considering a stiffe-
 * ned gas equation of state.
 *
 * Parameters:
 * -----------
 *  struct phase phase    : phase structure defined in "field.h";
 *  double       a        : volume fraction of phase;
 *  double       a_rho    : phasic mass (volume fractions * phasic density);
 *  double       a_rho_u  : phasic momentum (volume fractions * phasic density *
 *                          phasic velocity);
 *  double       a_rho_et : phasic total energy (volume fractions * phasic den-
 *                          sity * phasic total energy);
 */
double stiffenedGasSoundSpeed(
    phase phase, 
    double a, double a_rho, double a_rho_u, double a_rho_et
)
{
    //std::cout << "stiffenedGasSoundSpeed" << std::endl;

    double c2;
    c2 = (phase.gamma-1)*(phase.gamma*(a_rho_et-0.5*a_rho_u*a_rho_u/a_rho) \
       - a*phase.gamma*phase.pi)/a_rho;

    // Prevent negative speed of sound:
    double epsilon = 1.0e-8;
    c2 = -(c2-fabs(c2))/(2.0*fabs(c2))*epsilon + (c2+fabs(c2))/2.0;

    return sqrt(c2);
}

void computePressures(phase phase1, phase phase2, double *&Qn, double *&p)
{
    //std::cout << "computePressures" << std::endl;
    
    double rho_1, u_1, e_1, a_1, p_1;
    double rho_2, u_2, e_2, a_2, p_2;
    double g_1, pi_1, g_2, pi_2;

    g_1   = phase1.gamma;
    pi_1  = phase1.pi;
    g_2   = phase2.gamma;
    pi_2  = phase2.pi;

    a_1   = Qn[0];
    rho_1 = Qn[1]/Qn[0];
    u_1   = Qn[2]/Qn[1];
    e_1   = Qn[3]/Qn[1] - 0.5*u_1*u_1;
    //p_1   = (g_1 - 1.0)*rho_1*e_1 - g_1*pi_1;
    p_1   = (g_1 - 1.0)*(Qn[3] - 0.5*Qn[2]*u_1)/a_1 - g_1*pi_1;
    a_2   = 1.0-Qn[0];
    rho_2 = Qn[4]/(1.0-Qn[0]);
    u_2   = Qn[5]/Qn[4];
    e_2   = Qn[6]/Qn[4] - 0.5*u_2*u_2;
    //p_2   = (g_2 - 1.0)*rho_2*e_2 - g_2*pi_2;
    p_2   = (g_2 - 1.0)*(Qn[6] - 0.5*Qn[5]*u_2)/a_2 - g_2*pi_2;

    // Assign to values to pointer:
    p[0] = p_1;
    p[1] = p_2;
}

/* Method to interface pressure based on stiffened gas equation of state and the
 * following closure:
 *
 *      p_I = \sum_{l} \alpha_l p_l
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure defined in "field.h";
 *  struct phase phase2 : phase2 structure defined in "field.h";
 *  double*&     Qn     : reference to array of conservative variables at cell i;
 */
double computePressureInterface(phase phase1, phase phase2, double *&Qn)
{
    //std::cout << "computePressureInterface" << std::endl;

    double a1_p1, a2_p2;
    a1_p1 = (phase1.gamma-1)*(Qn[3]-0.5*Qn[2]*Qn[2]/Qn[1]) \
            - Qn[0]*phase1.gamma*phase1.pi;
    a2_p2 = (phase2.gamma-1)*(Qn[6]-0.5*Qn[5]*Qn[5]/Qn[4]) \
            - (1.0-Qn[0])*phase2.gamma*phase2.pi;

    // Avoid negative pressures:
    double epsilon = 1.0e-8;
    a1_p1 = -(a1_p1-fabs(a1_p1))/(2.0*fabs(a1_p1))*epsilon + (a1_p1+fabs(a1_p1))/2.0;
    a2_p2 = -(a2_p2-fabs(a2_p2))/(2.0*fabs(a2_p2))*epsilon + (a2_p2+fabs(a2_p2))/2.0;

    return a1_p1 + a2_p2;
}

/* Method to interface velocity based on the following closure:
 *
 *      u_I = \sum_{l} \alpha_l \rho_l u_l / \sum_{l} \alpha_l \rho_l
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure defined in "field.h";
 *  struct phase phase2 : phase2 structure defined in "field.h";
 *  double*&     Qn     : reference to array of conservative variables at cell i;
 */
double computeVelocityInterface(phase phase1, phase phase2, double *&Qn)
{
    //std::cout << "computeVelocityInterface" << std::endl;

    return (Qn[2]+Qn[5])/(Qn[1]+Qn[4]);
}
