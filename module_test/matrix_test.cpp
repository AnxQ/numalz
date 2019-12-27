#include <iostream>
#include "../Matrix.h"
#include "../LinearAlgebra.h"
#include <chrono>
using namespace std;
using namespace chrono;

auto last = system_clock::now();
auto t = []() -> float {
    auto now = system_clock::now();
    auto dur = duration_cast<microseconds>(now - last).count();
    last = now;
    return dur / 1000.0f;
};

int main() {
    Mat M({
        {9.1, 3.0, 3.6, 4.0},
        {4.2, 5.3, 4.7, 1.6},
        {3.2, 1.7, 9.4, 1.0},
        {6.1, 4.9, 3.5, 6.2}
    });
    cout << "Hello, World!" << std::endl;

    // Basic functions
    // Matrix indices range from [0, 0] -> [m - 1, n - 1]
    cout << "Eye matrix: " << endl << Mat::eye(5);
    cout << "Upper triangular: " << endl << M.triu(1);
    cout << "Lower triangular: " << endl << M.tril(1);
    cout << "M[:2,:1]: " << endl << M.slice(sq(0, 2, 0, 1));
    cout << "Diagonal: " << endl << M.diag();
    cout << "Transpose: " << endl << M.transpose();
    cout << "Cofactor(1, 2): " << endl << M.cofactor(1,2);

    // Basic Calculating functions
    cout << "M . M': " << endl << M * M.transpose();
    cout << "M + M': " << endl << M + M.transpose();
    cout << "M * 10: " << endl << M * 10;
    cout << "L1 Norm: " << M.norm(NormType::L1) << endl;
    cout << "F Norm: " << M.norm(NormType::Forbenius) << endl;

    using namespace LinearAlgebra;

    // Decomposing functions
    auto h1 = householder(M.slice(sq(0, 3, 0, 0)), 0);
    cout << "Householder by M[:,0]:" << endl << h1 * M;

    auto g1 = givens(0, 1, M.slice(sq(0, 3, 0, 0)));
    cout << "Givens by M[:,0]:" << endl << g1 * M;

    t(); auto qr_res = qr(M);
    cout << "QR in Householder: " << t() << "ms"<< endl << qr_res.first << qr_res.second;

    t(); qr_res = qr(M, QRMethod::Givens);
    cout << "QR in Givens: " << t() << "ms" << endl << qr_res.first << qr_res.second;

    t(); auto _inv = inv(M);
    cout << "M^-1: " << t() << "ms" << endl << _inv;


    t(); auto lu_res = lu(M);
    cout << "LU in Doolittle: " << t() << "ms" << endl << lu_res.first << lu_res.second;

    // Calculating functions
    t(); auto _det = det(M);
    cout << "Determinant: " << t() << "ms" << endl << _det << endl;

    t(); auto eigs = eig(M);
    cout << "Eigens: " << t() << "ms" << endl; for (auto eig :eigs) cout << eig << " "; cout << endl;

    return 0;
}