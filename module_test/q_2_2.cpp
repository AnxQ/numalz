//
// Created by 安啸琪 on 2019/12/25.
//

#include <iostream>
#include "../Matrix.h"
#include "../LinearAlgebra.h"

using namespace std;
using namespace LinearAlgebra;

int main() {
    Mat M({
                  {-0.0001, 5.096, 5.101, 1.853},
                  {0.0,     3.737, 3.740, 3.392},
                  {0.0,     0.0,   0.006, 5.254},
                  {0.0,     0.0,   0.0,   4.567}
          });
    cout << "Eigens:" << endl; for (auto eig: eig(M)) cout << eig << "  "; cout << endl;
    cout << "Max Eigen by Model: " << eig_max(M, Mat({1,1,1,1})) << endl;
    cout << "Min Eigen by Model: " << eig_min(M, Mat({1,1,1,1}),true);

    return 0;
};