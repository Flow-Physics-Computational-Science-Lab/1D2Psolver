#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "log.h"
#include "field.h"

/* Method to print an 1d array.
 *
 * Parameters:
 * -----------
 *  int      n   : number of items;
 *  double*& arr : reference of 1d array of pointers;
 */
void print1DArray(int n, double *&arr) 
{
    std::cout << "print1DArray" << std::endl;
    for (int i=0; i<n; i++) {
        std::cout << arr[i] << ", ";
    }
    std::cout << std::endl;
}

/* Method to print an 2d array.
 *
 * Parameters:
 * -----------
 *  int      nx  : number of rows;
 *  int      iy  : column index to print;
 *  double*& arr : reference of 1d array of pointers;
 */
void print2DArray(int nx, int iy, double**& arr) 
{
    std::cout << "print2DArray" << std::endl;
    for (int i=0; i<nx; i++) {
        std::cout << arr[i][iy] << ", ";
    }
    std::cout << std::endl;
}

/* Method to write binary file of an 1d array.
 *
 * Parameters:
 * -----------
 *  const char*& fname : reference of string literal of file name;
 *  int          n     : number of items;
 *  double*&     arr   : reference of 1d array of pointers;
 */
void writeBinary1DArray(const char*& fname, int n, double*& arr)
{
    std::cout << "writeBinary1DArray" << std::endl;
    
    // Open a binary file for writing
    std::ofstream out(fname, std::ios::binary);

    if (!out.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    // Write the array to the binary file
    out.write(reinterpret_cast<const char*>(arr), n*sizeof(double));

    // Close the file
    out.close();
}

/* Method to write binary file of a 2d array.
 *
 * Parameters:
 * -----------
 *  const char*& fname : reference of string literal of file name;
 *  int          nx    : number of rows;
 *  int          ny    : number of columns;
 *  double**&    arr   : reference of 2d array of pointers;
 *
 * TODO:
 * -----
 *  - the way this is going to output is in form of a ravel row after row, I 
 *  would still need to reshape the array in python for example, however, that
 *  required prior knowledge of the dimensions of the array. Hence, ideally I
 *  should output in some other way that the metadata or the shape of the array
 *  is still written to the output.
 */
void writeBinary2DArray(std::string fname, int nx, int ny, double**& arr)
{
    std::cout << "writeBinaryArray" << std::endl;
    
    // Open a binary file for writing
    std::ofstream out(fname, std::ios::binary);

    if (!out.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    // Write the array to the binary file
    //out.write(reinterpret_cast<const char*>(arr), (nx*ny)*sizeof(double));

    // Write each row of the 2D array to the binary file
    for (int i=0; i<nx; i++) {
        out.write(reinterpret_cast<const char*>(arr[i]), ny*sizeof(double));
    }

    // Close the file
    out.close();
}

/* Method to write csv file of a 2d array.
 *
 * Parameters:
 * -----------
 *  const char*& fname : reference of string literal of file name;
 *  int          nx    : number of rows;
 *  int          ny    : number of columns;
 *  double**&    arr   : reference of 2d array of pointers;
 */
void writeCSV2DArray(const char*& fname, int nx, int ny, double**& arr)
{
    std::cout << "writeCSVArray" << std::endl;
    
    // Open a binary file for writing
    std::ofstream out(fname);

    if (!out.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }
    
    out << std::fixed << std::setprecision(10);

    // Loop through array and output to file:
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            out << arr[i][j] << ","; 
        }
        out << std::endl;
    }

    // Close the file
    out.close();
}

/* Method to print all variables inside the structure phase.
 *
 * Parameters:
 * -----------
 *  struct phase phase : structure phase defined in "field.h";
 */
void printPhase(phase phase)
{
    std::cout << "printPhase" << std::endl;
    //std::cout << std::fixed ;
    std::cout << "{rho = "  << phase.rho   << ", ";
    std::cout << "p = "     << phase.p     << ", ";
    std::cout << "u = "     << phase.u     << ", ";
    std::cout << std::setprecision(10) << "alpha = " << phase.alpha << ", ";
    std::cout << "gamma = " << phase.gamma << ", ";
    std::cout << "pi = "    << phase.pi    << "};" << std::endl;
}
