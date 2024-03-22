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
#include "thermo.h"

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
    double xI = 0.7;
    // c) Run time parameters:
    double **Qn, **Qh, **Qhu, **Qhup, **Qnp1, **Snf, **En, **Enf;
    RunTimeParameters sim_par;
    sim_par = (RunTimeParameters){
        //.nit = 11, .nrest = 1, .dt = 1.0e-10, .dx = dx, .CFL = 0.8};
        //.nit = 101, .nrest = 10, .dt = 1.0e-10, .dx = dx, .CFL = 0.8};
        .nit = 1001, .nrest = 100, .dt = 1.0e-10, .dx = dx, .CFL = 0.8};
        //.nit = 2501, .nrest = 500, .dt = 1.0e-10, .dx = dx, .CFL = 0.8};
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
    sprintf(buffer, "./out/Q%05d", 0);
    writeBinary2DArray(std::string(buffer), n+2*nhc-1, 7, Qn);

    // 4) Run simulation:
    allocate2d(n, 2, Snf);
    allocate2d(n+2*nhc-1, 6, En);
    allocate2d(n, 6, Enf);
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

        // Velocity and pressure relaxation:
        relaxationVelocity(n, nhc, Qn, Qh, Qhu);
        relaxationPressure2(phase1_0, phase2_0, n, nhc, Qn, Qh, Qhu, Qhup);

        //updateQn(n+2*nhc-1, Qn, Qh);
        //updateQn(n+2*nhc-1, Qn, Qhu);
        updateQn(n+2*nhc-1, Qn, Qhup);

        // Write file:
        if (i%sim_par.nrest == 0) {
        //if (i%sim_par.nrest == 0 && i>659 && i<671) {
            sprintf(buffer, "./out/Q%05d_h", i);
            writeBinary2DArray(std::string(buffer), n+2*nhc-1, 7, Qh); 
            sprintf(buffer, "./out/Q%05d_hu", i);
            writeBinary2DArray(std::string(buffer), n+2*nhc-1, 7, Qhu); 
            sprintf(buffer, "./out/Q%05d_hup", i);
            writeBinary2DArray(std::string(buffer), n+2*nhc-1, 7, Qhup); 
            sprintf(buffer, "./out/Q%05d", i);
            writeBinary2DArray(std::string(buffer), n+2*nhc-1, 7, Qn); 
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
