#include "mex.hpp"
#include "mexAdapter.hpp"

class MexFunction : public matlab::mex::Function {
public:
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    // get the vector and matrix from the inputs
    matlab::data::TypedArray<double> inVector = std::move(inputs[0]);
    matlab::data::TypedArray<double> inMatrix = std::move(inputs[1]);

    // get the size of the vector and matrix
    size_t vectorSize = inVector.getDimensions()[0];
    size_t numRows = inMatrix.getDimensions()[0];
    size_t numCols = inMatrix.getDimensions()[1];

    // get raw pointers from the TypedArray, assuming the inputs are 1D vector
    // and 2D matrix respectively
    double *vectorData = inVector.release().get();
    double *matrixData = inMatrix.release().get();

    for (size_t i = 0; i < vectorSize; ++i) {
      printf("vector elem %ld = %f\n", i, vectorData[i]);
    }

    for (size_t i = 0; i < numRows; ++i) {
      for (size_t j = 0; j < numCols; ++j) {
        printf("Matrix(%ld, %ld) = %f\n", i, j, matrixData[i + j * numRows]);
      }
    }

    // create a MATLAB scalar from the result
    matlab::data::ArrayFactory factory;
    outputs[0] = factory.createScalar(vectorData[0] + matrixData[0]);
  }
};