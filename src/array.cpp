#include <iostream>
#include "allocate.h"
#include "array.h"

/* Method to create an 1d array of zeros.
 *
 * Parameters:
 * -----------
 *  int      n   : number of rows;
 *  double*& arr : reference of 1d array of pointers;
 */
void arrayZeros(int n, double *&arr)
{
    std::cout << "arrayZeros" << std::endl;

    // Uses allocate method defined elsewhere;
    allocate1d(n, arr);

    // Sets all items to zero;
    for (int i=0; i<n; i++) {
        arr[i] = 0.0;
    }
}

/* Method to create an 1d array linearly spaced from x0 to xn with n nodes.
 *
 * Parameters:
 * -----------
 *  double   x0  : first item;
 *  double   xn  : last item;
 *  int      n   : number of items;
 *  double*& arr : reference of 1d array of pointers;
 */
void linSpace(double x0, double xn, int n, double *&arr)
{
    std::cout << "linSpace" << std::endl;

    // Uses allocate method defined elsewhere;
    allocate1d(n, arr);
    
    // Compute dx
    double dx;
    dx = (xn-x0)/(n-1);

    // Assings values for interior nodes;
    for (int i=1; i<(n-1); i++) {
        arr[i] = x0 + i*dx; 
    }
    // Assign values for exterior nodes;
    arr[0]   = x0;
    arr[n-1] = xn;
}
