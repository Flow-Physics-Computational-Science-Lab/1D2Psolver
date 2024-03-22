#include <iostream>
#include <cmath>
#include "relaxation.h"
#include "../field.h"
#include "../thermo.h"

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
        // with V_i^n as u_I_n (stage prior to integration), check sings, I think
        // it flips when writing the equation for internal energy;
        //e_1_hu =  e_1_h + 0.5*(u_I_hu-u_1_h)*(u_I_n-u_1_h);
        //e_2_hu =  e_2_h - 0.5*(u_I_hu-u_2_h)*(u_I_n-u_2_h);
        // signs flipped but with u_I_n;
        e_1_hu =  e_1_h - 0.5*(u_I_hu-u_1_h)*(u_I_n-u_1_h);
        e_2_hu =  e_2_h + 0.5*(u_I_hu-u_2_h)*(u_I_n-u_2_h);
        // with V_i^n as u_I_hu:
        //e_1_hu =  e_1_h + 0.5*(u_I_hu-u_1_h)*(u_I_hu-u_1_h);
        //e_2_hu =  e_2_h - 0.5*(u_I_hu-u_2_h)*(u_I_hu-u_2_h);
        //e_1_hu =  e_1_h - 0.5*(u_I_hu-u_1_h)*(u_I_hu-u_1_h);
        //e_2_hu =  e_2_h + 0.5*(u_I_hu-u_2_h)*(u_I_hu-u_2_h);

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
void relaxationPressure1(
    phase phase1, phase phase2,
    int n, int nhc, 
    double **&Qn, double **&Qh, double **&Qhu, double **&Qhup
)
{
    //std::cout << "relaxationPressure" << std::endl;

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
                std::cout << i-nhc+1 << ", ";
                //std::cout << "No root or multiple roots at " << i-nhc+1;
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
    std::cout << std::endl;

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
    e_1_s   = e_1_0 - p_I_s*(a_s - Qhu[0])/Qhu[1];
    e_2_s   = e_2_0 + p_I_s*(a_s - Qhu[0])/Qhu[4];
    // Update phasic pressures:
    p_1_s   = (phase1.gamma - 1.0)*rho_1_s*e_1_s + phase1.gamma*phase1.pi;
    p_2_s   = (phase2.gamma - 1.0)*rho_2_s*e_2_s + phase2.gamma*phase2.pi;

    e_k_s[0] = e_1_s;
    e_k_s[1] = e_2_s;

    return p_1_s - p_2_s;
}

double functionPresStar(
    phase phase1, phase phase2, 
    double*& Qhu, 
    double p_s
);

double functionPrimePresStar(
    phase phase1, phase phase2, 
    double*& Qhu, 
    double p_s
);

void computePropertiesStar(
    phase phase1, phase phase2,
    double*& Qhu,
    double p_s,
    double*& a_s, double*& e_k_s
);

void bisectionPressure(
    phase phase1, phase phase2,
    double*& Qhu, double*& Qhup,
    double*& p
);

void newtonPressure(
    phase phase1, phase phase2,
    double*& Qhu, double*& Qhup,
    double*& p
);

/* Method to relax each phase pressure to its equilibrium value, .
 *
 * Bisection method to solve for:
 *
 *      g(p*) = \alpha_1* + \alpha_2* - 1.0 = 0.0
 *
 * See Saurel et al (2009).
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
 *  double**&    Qhup    : reference to array of conservative variables after 
 *                         hyperbolic integration, velocity and pressure relaxa-
 *                         tion;
 */
void relaxationPressure2(
    phase phase1, phase phase2,
    int n, int nhc, 
    double **&Qn, double **&Qh, double **&Qhu, double **&Qhup
)
{
    //std::cout << "relaxationPressure" << std::endl;

    // Variables for bisection method:
    double err, max_err;
    double a, b;

    // Variables after relaxation:
    double* p = new double[2];

    for (int i=nhc; i<(n+nhc-1); i++) {
        computePressures(phase1, phase2, Qhu[i], p);
        //a   = fmin(p[0], p[1]); // - 100.0;
        //b   = fmax(p[0], p[1]); // + 100.0;
        //err = fabs(a-b);
        //c   = 0.5*(a+b); // p*
        /*
        if (err > max_err) {
            bisectionPressure(phase1, phase2, Qhu[i], Qhup[i], p);
        } else {
            Qhup[i][0] = Qhu[i][0];
            Qhup[i][1] = Qhu[i][1];
            Qhup[i][2] = Qhu[i][2];
            Qhup[i][3] = Qhu[i][3];
            Qhup[i][4] = Qhu[i][4];
            Qhup[i][5] = Qhu[i][5];
            Qhup[i][6] = Qhu[i][6];
        }
        */
        newtonPressure(phase1, phase2, Qhu[i], Qhup[i], p);
    }
    //std::cout << std::endl;

    computeBCs(n, nhc, Qhup);
    
    delete[] p;
}

void bisectionPressure(
    phase phase1, phase phase2,
    double*& Qhu, double*& Qhup,
    double*& p
)
{
    // Variables for bisection method:
    double err, max_err;
    int it, max_it;
    max_err = 100;
    max_it  = 50;
    double a, b, c;
    double fa, fb, fc;

    double* a_k_s = new double[2];
    double  u_k_s;
    double* e_k_s = new double[2];

    a   = fmin(p[0], p[1]); // - 100.0;
    b   = fmax(p[0], p[1]); // + 100.0;
    err = fabs(a-b);
    it  = 1;
    while (err > max_err && it < max_it) {
        c = 0.5*(a+b); // p*
        fa = functionPresStar(phase1, phase2, Qhu, a); 
        fc = functionPresStar(phase1, phase2, Qhu, c); 
        fb = functionPresStar(phase1, phase2, Qhu, b); 
        if (fa * fb >= 0.0) {
            //std::cout << "No root or multiple roots at " << i-nhc+1;
            //c = fmax(p[0], p[1]);
            computePropertiesStar(phase1, phase2, Qhu, c, a_k_s, e_k_s);
            std::cout << a << ", " << b << std::endl;
            std::cout << a_k_s[0] << ", " << a_k_s[1] << std::endl;
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
    } // c is the p*
    
    computePropertiesStar(phase1, phase2, Qhu, c, a_k_s, e_k_s);
    // Update conservative variables array:
    Qhup[0] = a_k_s[0];
    Qhup[1] = Qhu[1];
    Qhup[2] = Qhu[2];
    u_k_s   = Qhu[2]/Qhu[1];
    Qhup[3] = Qhu[1]*(e_k_s[0] + 0.5*u_k_s*u_k_s);
    Qhup[4] = Qhu[4];
    Qhup[5] = Qhu[5];
    u_k_s   = Qhu[5]/Qhu[4];
    Qhup[6] = Qhu[4]*(e_k_s[1] + 0.5*u_k_s*u_k_s);

    delete[] a_k_s;
    delete[] e_k_s;
}

void newtonPressure(
    phase phase1, phase phase2,
    double*& Qhu, double*& Qhup,
    double*& p
)
{
    // Variables for newton raphson method:
    double err, max_err;
    int it, max_it;
    max_err = 1.0e-3;
    max_it  = 50;
    double a, b, c;
    double fa, fb, fc;
    double p_s_n, p_s_np1;
    double f_p_s_n, fp_p_s_n;

    double* a_k_s = new double[2];
    double  u_k_s;
    double* e_k_s = new double[2];

    a   = fmin(p[0], p[1]); // - 100.0;
    b   = fmax(p[0], p[1]); // + 100.0;
    err = fabs(a-b);
    it  = 1;
    //p_s_n = 0.5*(a+b); // p*
    if (err > max_err) {
        p_s_n = computePressureInterface(phase1, phase2, Qhu); // p*
    } else {
        p_s_n = b;
    }
    while (err > max_err && it < max_it) {
        f_p_s_n  = functionPresStar(phase1, phase2, Qhu, p_s_n); 
        fp_p_s_n = functionPrimePresStar(phase1, phase2, Qhu, p_s_n); 
        p_s_np1 = p_s_n - f_p_s_n/fp_p_s_n;
        err = fabs(p_s_n-p_s_np1);
        p_s_n = p_s_np1;
        it += 1;
    }
    if (p_s_n < a || p_s_n > b) {
        std::cout << "Diverged" << std::endl;
    }
    
    computePropertiesStar(phase1, phase2, Qhu, p_s_n, a_k_s, e_k_s);
    // Update conservative variables array:
    Qhup[0] = a_k_s[0];
    Qhup[1] = Qhu[1];
    Qhup[2] = Qhu[2];
    u_k_s   = Qhu[2]/Qhu[1];
    Qhup[3] = Qhu[1]*(e_k_s[0] + 0.5*u_k_s*u_k_s);
    Qhup[4] = Qhu[4];
    Qhup[5] = Qhu[5];
    u_k_s   = Qhu[5]/Qhu[4];
    Qhup[6] = Qhu[4]*(e_k_s[1] + 0.5*u_k_s*u_k_s);

    delete[] a_k_s;
    delete[] e_k_s;
}

/* Method to evaluate function g(p*).
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure with EoS parameters;
 *  struct phase phase2 : phase1 structure with EoS parameters;
 *  double*&     Qhu    : reference to array of conservative variables after 
 *                        hyperbolic integration and velocity relaxation at cell i;
 *  double       p_s    : p* guessed;
 */
double functionPresStar(
    phase phase1, phase phase2, 
    double*& Qhu, 
    double p_s
)
{
    //std::cout << "functionPresStar" << std::endl;
    
    // Initial state, ie, after hyperbolic integration and velocity relaxation:
    double rho_1_0, u_1_0, e_1_0, a_1_0, p_1_0;
    double rho_2_0, u_2_0, e_2_0, a_2_0, p_2_0;
    double g_1, pi_1, g_2, pi_2;
    double p_I_0;

    g_1     = phase1.gamma;
    pi_1    = phase1.pi;
    g_2     = phase2.gamma;
    pi_2    = phase2.pi;

    a_1_0   = Qhu[0];
    //rho_1_0 = Qhu[1]/Qhu[0];
    u_1_0   = Qhu[2]/Qhu[1];
    //e_1_0   = Qhu[3]/Qhu[1] - 0.5*u_1_0*u_1_0;
    //p_1_0   = (g_1 - 1.0)*rho_1_0*e_1_0 - g_1*pi_1;
    p_1_0   = (g_1 - 1.0)*(Qhu[3] - Qhu[2]*u_1_0)/a_1_0 - g_1*pi_1;
    a_2_0   = 1.0-Qhu[0];
    //rho_2_0 = Qhu[4]/(1.0-Qhu[0]);
    u_2_0   = Qhu[5]/Qhu[4];
    //e_2_0   = Qhu[6]/Qhu[4] - 0.5*u_2_0*u_2_0;
    //p_2_0   = (g_2 - 1.0)*rho_2_0*e_2_0 - g_2*pi_2;
    p_2_0   = (g_2 - 1.0)*(Qhu[6] - Qhu[5]*u_2_0)/a_2_0 - g_2*pi_2;
    p_I_0   = a_1_0*p_1_0 + a_2_0*p_2_0;

    // Relaxed state
    double a_1_s, a_2_s, p_I_b;
    p_I_b = 0.5*(p_I_0 + p_s);
    a_1_s = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_I_0) \
          / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_I_0);
    a_2_s = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_I_0) \
          / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_I_0);
    //a_1_s = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_I_b) \
    //      / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_I_b);
    //a_2_s = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_I_b) \
    //      / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_I_b);
    //a_1_s = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_s) \
    //      / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_s);
    //a_2_s = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_s) \
    //      / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_s);

    return a_1_s + a_2_s - 1.0;
}

/* Method to evaluate function g'(p*).
 *
 * Parameters:
 * -----------
 *  struct phase phase1 : phase1 structure with EoS parameters;
 *  struct phase phase2 : phase1 structure with EoS parameters;
 *  double*&     Qhu    : reference to array of conservative variables after 
 *                        hyperbolic integration and velocity relaxation at cell i;
 *  double       p_s    : p* guessed;
 */
double functionPrimePresStar(
    phase phase1, phase phase2, 
    double*& Qhu, 
    double p_s
)
{
    //std::cout << "functionPresStar" << std::endl;
    
    // Initial state, ie, after hyperbolic integration and velocity relaxation:
    double rho_1_0, u_1_0, e_1_0, a_1_0, p_1_0;
    double rho_2_0, u_2_0, e_2_0, a_2_0, p_2_0;
    double g_1, pi_1, g_2, pi_2;
    double p_I_0;

    g_1     = phase1.gamma;
    pi_1    = phase1.pi;
    g_2     = phase2.gamma;
    pi_2    = phase2.pi;

    a_1_0   = Qhu[0];
    //rho_1_0 = Qhu[1]/Qhu[0];
    u_1_0   = Qhu[2]/Qhu[1];
    //e_1_0   = Qhu[3]/Qhu[1] - 0.5*u_1_0*u_1_0;
    //p_1_0   = (g_1 - 1.0)*rho_1_0*e_1_0 - g_1*pi_1;
    p_1_0   = (g_1 - 1.0)*(Qhu[3] - Qhu[2]*u_1_0)/a_1_0 - g_1*pi_1;
    a_2_0   = 1.0-Qhu[0];
    //rho_2_0 = Qhu[4]/(1.0-Qhu[0]);
    u_2_0   = Qhu[5]/Qhu[4];
    //e_2_0   = Qhu[6]/Qhu[4] - 0.5*u_2_0*u_2_0;
    //p_2_0   = (g_2 - 1.0)*rho_2_0*e_2_0 - g_2*pi_2;
    p_2_0   = (g_2 - 1.0)*(Qhu[6] - Qhu[5]*u_2_0)/a_2_0 - g_2*pi_2;
    p_I_0   = a_1_0*p_1_0 + a_2_0*p_2_0;

    // Relaxed state
    double a_1_s, a_2_s, p_I_b;
    p_I_b = 0.5*(p_I_0 + p_s);
    a_1_s = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_I_0) \
          / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_I_0);
    a_2_s = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_I_0) \
          / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_I_0);
    //a_1_s = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_I_b) \
    //      / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_I_b);
    //a_2_s = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_I_b) \
    //      / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_I_b);
    //a_1_s = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_s) \
    //      / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_s);
    //a_2_s = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_s) \
    //      / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_s);

    double fp_p_s_n;
    
    fp_p_s_n = -a_1_s/(p_s + g_1*pi_1 + (g_1 - 1.0)*p_I_0) \
               -a_2_s/(p_s + g_2*pi_2 + (g_2 - 1.0)*p_I_0);

    return fp_p_s_n;
}

void computePropertiesStar(
    phase phase1, phase phase2,
    double*& Qhu,
    double p_s,
    double*& a_k_s, double*& e_k_s
)
{
    //std::cout << "computePropertiesStar" << std::endl;
    
    // Initial state, ie, after hyperbolic integration and velocity relaxation:
    double rho_1_0, u_1_0, e_1_0, a_1_0, p_1_0;
    double rho_2_0, u_2_0, e_2_0, a_2_0, p_2_0;
    double g_1, pi_1, g_2, pi_2;
    double p_I_0;

    g_1     = phase1.gamma;
    pi_1    = phase1.pi;
    g_2     = phase2.gamma;
    pi_2    = phase2.pi;

    a_1_0   = Qhu[0];
    //rho_1_0 = Qhu[1]/Qhu[0];
    u_1_0   = Qhu[2]/Qhu[1];
    //e_1_0   = Qhu[3]/Qhu[1] - 0.5*u_1_0*u_1_0;
    //p_1_0   = (g_1 - 1.0)*rho_1_0*e_1_0 - g_1*pi_1;
    p_1_0   = (g_1 - 1.0)*(Qhu[3] - Qhu[2]*u_1_0)/a_1_0 - g_1*pi_1;
    a_2_0   = 1.0-Qhu[0];
    //rho_2_0 = Qhu[4]/(1.0-Qhu[0]);
    u_2_0   = Qhu[5]/Qhu[4];
    //e_2_0   = Qhu[6]/Qhu[4] - 0.5*u_2_0*u_2_0;
    //p_2_0   = (g_2 - 1.0)*rho_2_0*e_2_0 - g_2*pi_2;
    p_2_0   = (g_2 - 1.0)*(Qhu[6] - Qhu[5]*u_2_0)/a_2_0 - g_2*pi_2;
    p_I_0   = a_1_0*p_1_0 + a_2_0*p_2_0;

    // Relaxed state:
    double a_1_s, e_1_s, p_1_s, nu_1_s;
    double a_2_s, e_2_s, p_2_s, nu_2_s;
    double p_I_b, p, rho_e;
    
    p_I_b   = 0.5*(p_I_0 + p_s);
    a_1_s = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_I_0) \
          / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_I_0);
    //a_2_s = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_I_0) \
    //      / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_I_0);
    //a_1_s   = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_I_b) \
    //        / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_I_b);
    //a_2_s   = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_I_b) \
    //        / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_I_b);
    //a_1_s = (p_1_0*a_1_0 + g_1*pi_1*a_1_0 + (g_1 - 1.0)*a_1_0*p_s) \
    //      / (p_s + g_1*pi_1 + (g_1 - 1.0)*p_s);
    //a_2_s = (p_2_0*a_2_0 + g_2*pi_2*a_2_0 + (g_2 - 1.0)*a_2_0*p_s) \
    //      / (p_s + g_2*pi_2 + (g_2 - 1.0)*p_s);
    a_2_s = 1.0 - a_1_s;

    // Reinitialization step as in Saurel et al (2009):
    // Compute mixture pressure:
    //p = ((Qhu[3]+Qhu[6])-(a_1_s*g_1*pi_1/(g_1-1.0) + (1.0-a_1_s)*g_2*pi_2/(g_2-1.0))) \
    //  / (a_1_s/(g_1-1.0) + (1.0-a_1_s)/(g_2-1.0));
    // Need to remove kinetic energy here!!!
    rho_e = (Qhu[3]+Qhu[6]) - 0.5*(Qhu[2]+Qhu[5])*(Qhu[2]+Qhu[5])/(Qhu[1]+Qhu[4]);
    p = (rho_e - (a_1_s*g_1*pi_1/(g_1-1.0) + a_2_s*g_2*pi_2/(g_2-1.0))) \
      / (a_1_s/(g_1-1.0) + a_2_s/(g_2-1.0));

    // Need to use EOS:
    nu_1_s = a_1_s/Qhu[1];
    e_1_s  = p*nu_1_s/(g_1 - 1.0) + g_1*pi_1*nu_1_s/(g_1 - 1.0);
    //nu_2_s = (1.0-a_1_s)/Qhu[4];
    nu_2_s = a_2_s/Qhu[4];
    e_2_s  = p*nu_2_s/(g_2 - 1.0) + g_2*pi_2*nu_2_s/(g_2 - 1.0);

    a_k_s[0] = a_1_s;
    a_k_s[1] = a_2_s;
    e_k_s[0] = e_1_s;
    e_k_s[1] = e_2_s;
}
