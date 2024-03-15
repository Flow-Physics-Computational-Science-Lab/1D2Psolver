#include <iostream>
#include "grid.h"
#include "allocate.h"

/* Method to pad the 1d grid with nhc halo cells.
 *
 * Parameters:
 * -----------
 *  int      nhc : number of halo cells to add;
 *  int      n   : number of nodes;
 *  double*& xs  : x locations of grid nodes;
 *  double*& xhs : output array of grid nodes padded with halo cells;
 */
void addHaloCells(int nhc, int n, double *&xs, double *&xhs)
{
    std::cout << "addHaloCells" << std::endl;
    // Allocate memory for xhs array:
    allocate1d(n+2*nhc, xhs); 

    // Copy interior nodes:
    for (int i=0; i<n; i++) {
        xhs[nhc+i] = xs[i];
    }

    // Add halo nodes:
    for (int i=0; i<nhc; i++) {
        xhs[nhc-(i+1)] = xs[0]   - (i+1)*(xs[1]-xs[0]);
        xhs[nhc+(i+n)] = xs[n-1] + (i+1)*(xs[n-1]-xs[n-2]);
    }
}

/* Method to compute cell centers based on 1d array.
 *
 * Parameters:
 * -----------
 *  int      n   : number of nodes;
 *  double*& xs  : x locations of grid nodes;
 *  double*& xcs : output array of cell centers with n-1 cells;
 */
void computeCellCenters(int n, double*& xs, double*& xcs)
{
    std::cout << "computeCellCenters" << std::endl;
    // Allocate memory for xcs array:
    allocate1d(n-1, xcs);
    
    // Average of location of grid nodes;
    for (int i=0; i<(n-1); i++) {
        xcs[i] = 0.5*(xs[i+1] + xs[i]);
    }
}
