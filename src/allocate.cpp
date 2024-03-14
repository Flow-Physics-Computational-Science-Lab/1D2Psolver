#include <iostream>
#include "allocate.h"

/* Method to allocate memory for a 1D array.
 *
 * Parameters:
 * -----------
 *  int      n   : number of rows;
 *  double*& arr : reference of 1d array of pointers;
 */
void allocate1d(int n, double*& arr)
{
    std::cout << "allocate1d" << std::endl;
    // Create an array of pointer with length n;
    double* newarr;
    newarr = new double[n];
    // Assign memory allocated to the reference of the array of pointers;
    arr = newarr;
}

/* Method to allocate memory for a 2D array. Using contiguous memory allocation.
 *
 * Parameters:
 * -----------
 *  int      nx   : number of rows;
 *  int      ny   : number of columns;
 *  double**& arr : reference of 2d array of pointers;
 */
void allocate2d(int nx, int ny, double**& arr)
{
    std::cout << "allocate2d" << std::endl;
    // Create a two dimensional array of pointers with shape nx*ny;
    double** newarr = new double*[nx];
    newarr[0] = new double[nx*ny];
    for (int i=1; i<nx; i++) {
        newarr[i] = &newarr[0][i*ny];
    }
    // Assign memory allocated to the reference of the array of pointers;
    arr = newarr;
}

/* Method to deallocate memory for a 1D array.
 *
 * Parameters:
 * -----------
 *  double*& arr : reference of 1d array of pointers;
 */
void deallocate1d(double*& arr)
{
    std::cout << "deallocate1d" << std::endl;
    delete[] arr;
}

/* Method to deallocate memory for a 2D array.
 *
 * Parameters:
 * -----------
 *  double*& arr : reference of 1d array of pointers;
 */
void deallocate2d(double**& arr)
{
    std::cout << "deallocate2d" << std::endl;
    delete[] arr[0];
    delete[] arr;
}
