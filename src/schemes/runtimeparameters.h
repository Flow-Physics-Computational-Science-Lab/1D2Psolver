#ifndef RUNPAR
#define RUNPAR

struct RunTimeParameters {
    int nit, nrest;
    double dt, dx;
    double CFL;
};

#endif
