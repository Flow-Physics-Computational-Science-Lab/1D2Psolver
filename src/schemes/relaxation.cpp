#include <iostream>
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
    std::cout << "relaxationVelocity" << std::endl;

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
