//
// Created by 安啸琪 on 2019/12/24.
//

#include <iostream>
#include "../Matrix.h"
#include "../LinearAlgebra.h"

using namespace std;
using namespace LinearAlgebra;

int main() {
    Mat M({
        {9.1, 3.0, 3.6, 4.0},
        {4.2, 5.3, 4.7, 1.6},
        {3.2, 1.7, 9.4, 1.1},
        {6.1, 4.9, 3.5, 6.2}
    });
    auto qr_res = qr(M);
    cout << "QR in Householder:" << endl << qr_res.first << qr_res.second;
    cout << "Eigens: "; for (auto eig :eig(M)) cout << eig << " ";

    return 0;
};