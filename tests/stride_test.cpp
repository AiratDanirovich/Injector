#include <array>
#include <iostream>

#include <Eigen/Core>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace Eigen;

// Tests Cartesian grid
TEST_CASE("EigenStride")
{
    int array[30ull];
    for (int i = 0; i < 30; ++i)
        array[i] = i;

    cout << "Column-major:\n"
         << Map<Matrix<int, 3, 10>>(array) << endl;
    // cout << "Row-major:\n"
    //      << Map<Matrix<int, 1, 4, RowMajor>>(array) << endl;
     cout << "Row-major:\n"
          << Map<Matrix<int, 10, 1>, Unaligned, Stride<1, 3>>(array) << endl;

          using Stride_t = InnerStride<Dynamic>;
          using Map_t = Map<VectorXi, Unaligned, Stride_t>;

          Stride_t stride{3};

          Map_t map{array+1, 10, stride};
     cout << "Row-major 2:\n";
          cout << map << endl;

    // cout << "Row-major using stride:\n"
    //      << Map<Matrix<int, 2, 4>, Unaligned, Stride<1, 4>>(array) << endl;
}
