#include <cuda_runtime.h>
#include <stdio.h>
#include <cuda.h>

#include <cassert>
#include <iostream>
#include <memory>
#include <eigen3/Eigen/Dense>

__device__ void EigenCheck(int *a, int N) {
    for (unsigned i = 0 ; i < N ; ++i) {
        printf("a[%d] = %d \n",i,a[i]);
    }

    printf("\n Checking Eigen identity matrix \n" );
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(N,N);
    Id(N-1,N-1) = N;

    for (unsigned i = 0 ; i < N ; ++i) {
        for (unsigned j = 0 ; j < N ; ++j) {
            printf("%.2f ",Id(i,j));
        }
        printf("\n");
    }

    printf("\n Checking another Eigen matrix \n" );
    Eigen::MatrixXd A(N,2);

    for (unsigned i = 0 ; i < N ; ++i) {
        A(i,0) = a[i];
    }
    A(0,0) = 15;
    for (unsigned i = 0 ; i < N ; ++i) {
        printf("%d ",A(i,0));
    }
    printf("\n\n");

    printf("\n Checking Eigen vector \n" );
    Eigen::VectorXi vec(N);

    for (unsigned i = 0 ; i < N ; ++i) {
        vec(i) = a[i];
    }

    for (unsigned i = 0 ; i < N ; ++i) {
        printf("%d ",vec(i));
    }
    printf("\n\n");
}

__device__ void LolCheck(int *a, int N) {
    
    printf("\n Checking another Eigen matrix \n" );
    Eigen::Map<Eigen::MatrixXi> A(a,N,1);

    for (unsigned i = 0 ; i < N ; ++i) {
        for (unsigned j = 0 ; j < 1 ; ++j) {
            printf("%d ",A(i,j));
        }
        printf("\n");
    }
}

__device__ void manualcheck(int N) {
    Eigen::MatrixXd A(N,1);

    //A << 1,2,3,4,5,6,7,8,9;
    A(0,0) = 3030;
    for (unsigned i = 0 ; i < N ; ++i) {
        for (unsigned j = 0 ; j < 1 ; ++j) {
            printf("%.2f ",A(i,j));
        }
        printf("\n");
    }
}


__global__ void testfn(int *a, int *b, int *c,int N) {
    printf("Inside block, thread: %d,%d with value N = %d \n",blockIdx.x,threadIdx.x,N);
    
    /* for (int i = 0 ; i < N ; ++i) {
        printf("loop idx %d \n",i);
    } */
    
    for (int i = 0 ; i < N ; ++i) {
        //printf("loop idx %d \n",i);
        c[i] = a[i] + b[i];
        if (threadIdx.x == 0) {
            printf("c[%i] = %i\n\n",i,c[i]);
        }
    }
    if (threadIdx.x == 0) {

    //EigenCheck(a,N);
    LolCheck(a,N);
    //manualcheck(N);
    }

}

int main() {
    // Adding vectors 
    int* a = new int[5];
    int* b = (int*)std::malloc(5 * sizeof(int));
    int *c = (int*)std::malloc(5*sizeof(int));
    for (int i = 0 ; i < 5 ; ++i) {
        a[i] = 1+i;
        b[i] = 2+i;
        std::cout << a[i] << " ";
    } 
    std::cout << std::endl;
    // Pointers for vectors on device
    int* da, *db, *dc;
    // Allocating memory on device
    cudaMalloc((void**)&da,5*sizeof(int));
    cudaMalloc((void**)&db,5*sizeof(int));
    cudaMalloc((void**)&dc,5*sizeof(int));

    // Transferring the data
    cudaMemcpy(da,a,5*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,5*sizeof(int),cudaMemcpyHostToDevice);

    // Adding vectors on the gpu
    int N = 5;
    //testfn<<<1,5>>> (a,b,c,N); // Problematic one
    testfn<<<1,3>>> (da,db,dc,N);
    //helloFromGPU<<<5,1>>>();

    cudaMemcpy(c,dc,5*sizeof(int),cudaMemcpyDeviceToHost);

    for (int i = 0 ; i < 5 ; ++i) {
        std::cout << c[i] << " ";
    }
    std::cout << std::endl;

    // Free CPU memory
    delete[] a;
    std::free(b);
    std::free(c);

    // Free GPU memory
    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
 
    return 0;
}