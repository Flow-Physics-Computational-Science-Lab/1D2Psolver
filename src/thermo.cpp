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
    c2 = (phase.gamma-1)*((phase.gamma-2)*(a_rho_et-0.5*a_rho_u*a_rho_u/a_rho) \
       - a*phase.gamma*phase.pi)/a_rho;

    // Prevent negative speed of sound:
    double epsilon = 1.0e-8;
    c2 = -(c2-fabs(c2))/(2.0*fabs(c2))*epsilon + (c2+fabs(c2))/2.0;

    return sqrt(c2);
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
    a1_p1 = (phase1.gamma-1)*(Qn[3]-0.5*Qn[2]*Qn[2]/Qn[1] \
            - Qn[0]*phase1.gamma*phase1.pi);
    a2_p2 = (phase2.gamma-1)*(Qn[6]-0.5*Qn[5]*Qn[5]/Qn[4] \
            - (1.0-Qn[0])*phase2.gamma*phase2.pi);

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
