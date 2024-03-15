#include <iostream>
#include "allocate.h"
#include "array.h"
#include "grid.h"
#include "field.h"
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
    double **Qn, **Qnp1;

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
    const char* file_Q0 = "./out/Q00000";
    writeBinary2DArray(file_Q0, n+2*nhc-1, 7, Qn);
    //writeCSV2DArray(file_Q0, n+2*nhc-1, 7, Qn);
    
    // 4) Run simulation:

    // 5) Deallocate arrays:
    deallocate1d(xs);
    deallocate1d(xcs);
    deallocate1d(xhs);
    deallocate2d(Qn);
    deallocate2d(Qnp1);

    std::cout << "end" << std::endl;
    return 0;
}
