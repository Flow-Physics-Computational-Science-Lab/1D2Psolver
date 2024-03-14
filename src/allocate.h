#ifndef ALLOCATE
#define ALLOCATE

void allocate1d  (int n, double*& arr);
void allocate2d  (int nx, int ny, double**& arr);
void deallocate1d(double*& arr);
void deallocate2d(double**& arr);

#endif
