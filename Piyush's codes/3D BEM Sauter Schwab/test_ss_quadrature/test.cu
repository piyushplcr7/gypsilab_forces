#include <cuda_runtime.h>
#include <stdio.h>
#include <cuda.h>
#include <math.h>

#include <cassert>
#include <iostream>
#include <memory>
#include <eigen3/Eigen/Dense>

__global__ void add1(double *a, int Nx, int Ny, double b)
{
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {

            //*(*(a + i) + j) += b;
            a[i + Nx * j] += b;
        }
    }
}