#ifndef LOG
#define LOG

#include "field.h"

void print1DArray      (int n, double*& arr);
void print2DArray      (int nx, int iy, double**& arr);
void writeBinary1DArray(const char*& fname, int n, double*& arr);
void writeBinary2DArray(const char*& fname, int nx, int ny, double**& arr);
void writeCSV2DArray   (const char*& fname, int nx, int ny, double**& arr);
void printPhase        (phase phase);

#endif
