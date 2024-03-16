#include <iostream>
#include <string>
#include "allocate.h"
#include "array.h"
#include "grid.h"
#include "field.h"
#include "schemes.h"
#include "log.h"

struct RunTimeParameters {
    int nit, nrest;
    double dt;
};

int main() 
{
    // 1) Parameters:
    // a) Grid:
    double *xs, *xhs,*xcs;
    int n = 101;
    double x0=0.0, xn=1.0;
    double dx=(xn-x0)/(n-1);
    const char* file_grid = "./out/grid";
    int nhc = 1;
    // b) Initial condition:
    double epsilon = 1.0e-8;
    phase phase1_0, phase2_0;
    phase1_0 = (phase){ 
        .rho = 1.0e3, .p = 1.0e9, .u = 0.0, .alpha=(1-epsilon), 
        .gamma = 4.4, .pi = 6.0e8, .epsilon = epsilon
    };
    phase2_0 = (phase){ 
        .rho = 5.0e1, .p = 1.0e5, .u = 0.0, .alpha=(1-epsilon), 
        .gamma = 1.4, .pi = 0.0, .epsilon = epsilon
    };
    //printPhase(phase1_0);
    //printPhase(phase2_0);
    double xI = 0.7;
    // c) Run time parameters:
    double **Qn, **Qnp1, **Snf, **En, **Enf;
    RunTimeParameters sim_par;
    sim_par = (RunTimeParameters){.nit = 101, .nrest = 100, .dt = 1.0e-5};
    //std::string file_Qn;
    char buffer[20];

    // 2) Initialize grid:
    linSpace(x0, xn, n, xs);
    writeBinary1DArray(file_grid, n, xs);
    addHaloCells(nhc, n, xs, xhs);
    //const char* file_grid_whc = "./out/grid_whc";
    //writeBinary1DArray(file_grid_whc, n+2*nhc, xhs);
    computeCellCenters(n+2*nhc, xhs, xcs);

    // 3) Initialize flow field:
    allocate2d(n+2*nhc-1, 7, Qn);
    allocate2d(n+2*nhc-1, 7, Qnp1);
    initializeWaterAirShockTube(phase1_0, phase2_0, xI, n, nhc, xcs, Qn);
    //const char* file_Q0 = "./out/Q00000";
    sprintf(buffer, "./out/Q%05d", 0);
    std::string file_Qn(buffer);
    writeBinary2DArray(file_Qn, n+2*nhc-1, 7, Qn);
    //writeCSV2DArray(file_Q0, n+2*nhc-1, 7, Qn);
    
    // 4) Run simulation:
    allocate2d(n+2*nhc-2, 2, Snf);
    allocate2d(n+2*nhc-1, 6, En);
    allocate2d(n+2*nhc-2, 6, Enf);
    for (int i=0; i<sim_par.nit; i++) {
        // Hyperbolic Gudunov:
        advanceTimeHyperbolicGudunov(
            phase1_0, phase2_0, 
            sim_par.dt, dx,
            n+2*nhc-1, 
            Qn, Qnp1,
            Snf, En, Enf
        );
        // Pressure and velocity relaxation:
        //relaxationPresVel();
        // Update Qn as the computed Qn+1:
        //updateQn();
        if (i%sim_par.nrest == 0) {
            //file_Qn = "./out/Q" + std::format("{i:05d}"); 
            sprintf(buffer, "./out/Q%05d", i);
            std::string file_Qn(buffer);
            writeBinary2DArray(file_Qn, n+2*nhc-1, 7, Qn); 
        }
    }

    // 5) Deallocate arrays:
    deallocate1d(xs);
    deallocate1d(xcs);
    deallocate1d(xhs);
    deallocate2d(Qn);
    deallocate2d(Qnp1);
    deallocate2d(Snf);
    deallocate2d(En);
    deallocate2d(Enf);

    std::cout << "end" << std::endl;
    return 0;
}
