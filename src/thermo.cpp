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
    std::cout << "stiffenedGasSoundSpeed" << std::endl;

    double c2;
    c2 = (phase.gamma-1)*((phase.gamma-2)*(a_rho_et-0.5*a_rho_u*a_rho_u/a_rho) \
       - a*phase.gamma*phase.pi)/a_rho;

    // Prevent negative speed of sound:
    double epsilon = 1.0e-8;
    c2 = -(c2-fabs(c2))/(2.0*fabs(c2))*epsilon + (c2+fabs(c2))/2.0;

    return sqrt(c2);
}
