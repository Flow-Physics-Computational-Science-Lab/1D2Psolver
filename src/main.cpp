#include <iostream>
#include "allocate.h"
#include "array.h"

int main() 
{
    double* xs;
    int n = 101;

    //arrayZeros(n, xs);
    linSpace(0.0, 1.0, n, xs);
    deallocate1d(xs);

    std::cout << "end" << std::endl;
    return 0;
}
