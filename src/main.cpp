#include <iostream>
#include <stdexcept>
#include <string>
#include "allocate.h"
#include "array.h"
#include "grid.h"
#include "field.h"
#include "runtimeparameters.h"
#include "schemes/gudunov.h"
#include "schemes/relaxation.h"
#include "log.h"

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
    double **Qn, **Qh, **Qhu, **Qhup, **Qnp1, **Snf, **En, **Enf;
    RunTimeParameters sim_par;
    sim_par = (RunTimeParameters){
        .nit = 11, .nrest = 1, .dt = 1.0e-10, .dx = dx, .CFL = 0.8};
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
    allocate2d(n+2*nhc-1, 7, Qh);
    allocate2d(n+2*nhc-1, 7, Qhu);
    allocate2d(n+2*nhc-1, 7, Qhup);
    allocate2d(n+2*nhc-1, 7, Qnp1);
    initializeWaterAirShockTube(phase1_0, phase2_0, xI, n, nhc, xcs, Qn);
    //const char* file_Q0 = "./out/Q00000";
    sprintf(buffer, "./out/Q%05d", 0);
    std::string file_Qn(buffer);
    writeBinary2DArray(file_Qn, n+2*nhc-1, 7, Qn);
    //writeCSV2DArray(file_Q0, n+2*nhc-1, 7, Qn);
    
    // 4) Run simulation:
    allocate2d(n, 2, Snf);
    //allocate2d(n+2*nhc-2, 2, Snf);
    allocate2d(n+2*nhc-1, 6, En);
    allocate2d(n, 6, Enf);
    //allocate2d(n+2*nhc-2, 6, Enf);
    std::cout << "Iterations:" << std::endl;
    for (int i=1; i<sim_par.nit; i++) {
        std::cout << i << ", ";

        // Hyperbolic Gudunov:
        try {
            advanceTimeHyperbolicGudunov(
                phase1_0, phase2_0, 
                sim_par,
                n, nhc, 
                Qn, Qh,
                Snf, En, Enf
            );
        } catch (const std::runtime_error) {
            deallocate1d(xs);
            deallocate1d(xcs);
            deallocate1d(xhs);
            deallocate2d(Qn);
            deallocate2d(Qh);
            deallocate2d(Qhu);
            deallocate2d(Qhup);
            deallocate2d(Qnp1);
            deallocate2d(Snf);
            deallocate2d(En);
            deallocate2d(Enf);
            throw std::runtime_error("Insufficient dt for prescribed CFL.");
            //return 1;
        }

        // TODO:
        // -----
        //  - improve this multiple copy of Qn+1 to Qn;
        // Update Qn as the computed Qn+1:
        //updateQn(n+2*nhc-1, Qn, Qnp1);

        // Velocity and pressure relaxation:
        relaxationVelocity(n, nhc, Qn, Qh, Qhu);
        //updateQn(n+2*nhc-1, Qn, Qnp1);
        //relaxationPressure();

        updateQn(n+2*nhc-1, Qn, Qhu);
        //updateQn(n+2*nhc-1, Qn, Qnp1);

        // Write file:
        if (i%sim_par.nrest == 0) {
            //file_Qn = "./out/Q" + std::format("{i:05d}"); 
            sprintf(buffer, "./out/Qh%05d", i);
            writeBinary2DArray(std::string(buffer), n+2*nhc-1, 7, Qh); 
            sprintf(buffer, "./out/Qhu%05d", i);
            //std::string file_Qn(buffer);
            writeBinary2DArray(std::string(buffer), n+2*nhc-1, 7, Qhu); 
        }
    }
    std::cout << std::endl;

    // 5) Deallocate arrays:
    deallocate1d(xs);
    deallocate1d(xcs);
    deallocate1d(xhs);
    deallocate2d(Qn);
    deallocate2d(Qh);
    deallocate2d(Qhu);
    deallocate2d(Qhup);
    deallocate2d(Qnp1);
    deallocate2d(Snf);
    deallocate2d(En);
    deallocate2d(Enf);

    std::cout << "end" << std::endl;
    return 0;
}
