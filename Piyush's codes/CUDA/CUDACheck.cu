#include <stdio.h>
#include <cuda_runtime.h>

int main() {
    int device;
    cudaGetDevice(&device);

    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, device);

    printf("Maximum block dimensions: %d x %d x %d\n",
           props.maxThreadsDim[0], props.maxThreadsDim[1], props.maxThreadsDim[2]);
    printf("Maximum grid dimensions: %d x %d x %d\n",
           props.maxGridSize[0], props.maxGridSize[1], props.maxGridSize[2]);

    return 0;
}