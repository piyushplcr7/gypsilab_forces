#include "mex.hpp"
#include "mexAdapter.hpp"
#include <eigen3/Eigen/Dense>
#include <vector>

class MexFunction : public matlab::mex::Function {
public:
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    printf("Entered the function \n");
    // Parsing the inputs into matlab data type
    matlab::data::TypedArray<int> _TrialDim = std::move(inputs[0]);
    printf("wot \n");
    int TrialDim = _TrialDim[0];
    printf("TrialDim = %d", TrialDim);
    /* matlab::data::TypedArray<int> _TestDim = std::move(inputs[1]);
    int TestDim = _TestDim[0];
    matlab::data::TypedArray<int> _NTriangles = std::move(inputs[2]);
    int NTriangles = _NTriangles[0];
    matlab::data::TypedArray<int> _NVertices = std::move(inputs[3]);
    int NVertices = _NVertices[0];
    matlab::data::TypedArray<int> _NInteractions = std::move(inputs[4]);
    int NInteractions = _NInteractions[0];
    matlab::data::TypedArray<int> I = std::move(inputs[5]);
    matlab::data::TypedArray<int> J = std::move(inputs[6]);
    matlab::data::TypedArray<int> relation = std::move(inputs[7]);
    matlab::data::TypedArray<double> _W0 = std::move(inputs[8]);
    double *W0 = &*_W0.begin();
    matlab::data::TypedArray<double> _X0 = std::move(inputs[9]);
    double *X0 = &*_X0.begin();
    matlab::data::TypedArray<int> _Nq0 = std::move(inputs[10]);
    int Nq0 = _Nq0[0];
    matlab::data::TypedArray<double> _W1 = std::move(inputs[11]);
    double *W1 = &*_W1.begin();
    matlab::data::TypedArray<double> _X1 = std::move(inputs[12]);
    double *X1 = &*_X1.begin();
    matlab::data::TypedArray<int> _Nq1 = std::move(inputs[13]);
    int Nq1 = _Nq1[0];
    matlab::data::TypedArray<double> _W2 = std::move(inputs[14]);
    double *W2 = &*_W2.begin();
    matlab::data::TypedArray<double> _X2 = std::move(inputs[15]);
    double *X2 = &*_X2.begin();
    matlab::data::TypedArray<int> _Nq2 = std::move(inputs[16]);
    int Nq2 = _Nq2[0];
    matlab::data::TypedArray<double> _W3 = std::move(inputs[17]);
    double *W3 = &*_W3.begin();
    matlab::data::TypedArray<double> _X3 = std::move(inputs[18]);
    double *X3 = &*_X3.begin();
    matlab::data::TypedArray<int> _Nq3 = std::move(inputs[19]);
    int Nq3 = _Nq3[0];
    matlab::data::TypedArray<double> trial_vec = std::move(inputs[20]);
    matlab::data::TypedArray<double> test_vec = std::move(inputs[21]);
    matlab::data::TypedArray<int> Elements = std::move(inputs[22]);
    matlab::data::TypedArray<double> Vertices = std::move(inputs[23]);
    matlab::data::TypedArray<double> Normals = std::move(inputs[24]);
    matlab::data::TypedArray<double> Areas = std::move(inputs[25]);
    matlab::data::TypedArray<int> Elt2DofTest = std::move(inputs[26]);
    matlab::data::TypedArray<int> Elt2DofTrial = std::move(inputs[27]);
    matlab::data::TypedArray<int> _TrialSpace = std::move(inputs[28]);
    int TrialSpace = _TrialSpace[0];
    matlab::data::TypedArray<int> _TestSpace = std::move(inputs[29]);
    int TestSpace = _TestSpace[0];
    matlab::data::TypedArray<int> _TrialOperator = std::move(inputs[30]);
    int TrialOperator = _TrialOperator[0];
    matlab::data::TypedArray<int> _TestOperator = std::move(inputs[31]);
    int TestOperator = _TestOperator[0];
    matlab::data::TypedArray<int> _NRSFTrial = std::move(inputs[32]);
    int NRSFTrial = _NRSFTrial[0];
    matlab::data::TypedArray<int> _NRSFTest = std::move(inputs[33]);
    int NRSFTest = _NRSFTest[0]; */
  }
};
