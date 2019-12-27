//
// Created by 安啸琪 on 2019/12/22.
//

#ifndef NUMALZ_SOLVER_H
#define NUMALZ_SOLVER_H

#include <utility>
#include "Matrix.h"
#include "LinearAlgebra.h"

using namespace LinearAlgebra;

class S_exception : exception {
    string msg;
public:
    explicit S_exception(string msg) : msg(std::move(msg)) {}

    const char *what() { return msg.c_str(); }
};

class BaseSolver {
public:
    Mat A;
    Mat b;
    size_t dim;

    BaseSolver(Mat A, Mat b) : A(std::move(A)), b(std::move(b)), dim(b.shape().first) {
        auto shape = A.shape();
        if (shape.first != dim || shape.first != shape.second)
            throw S_exception("Input dimension mismatch");
    }

    virtual Mat solve() = 0;
};

class GaussBaseSolver : public BaseSolver {
public:
    GaussBaseSolver(Mat A, Mat b) : BaseSolver(std::move(A), std::move(b)) {}

    virtual pair<Mat, Mat> Triangularize() = 0;

    virtual Mat BackSubstitute(const pair<Mat, Mat> &t_Ab) const {
        size_t dim = b.shape().first;
        Mat res(b.shape());
        for (int i = dim - 1; i >= 0; --i)
            res[i][0] = (t_Ab.second[i][0] - (t_Ab.first.slice(sq(i, i, 0, dim - 1)) * res)[0][0]) / t_Ab.first[i][i];
        return res;
    }

    Mat solve() override {
        return BackSubstitute(Triangularize());
    }
};

class GaussSolver : public GaussBaseSolver {
public:
    GaussSolver(Mat A, Mat b) : GaussBaseSolver(std::move(A), std::move(b)) {}

    pair<Mat, Mat> Triangularize() override {
        auto Ab = hstack(A, b);
        Mat l_inv;
        for (size_t i = 0; i < dim; ++i) {
            // Form a Gauss transforming matrix
            l_inv = Mat::eye(dim);
            for (size_t j = i + 1; j < dim; ++j) l_inv[j][i] = -Ab[j][i] / Ab[i][i];
            Ab = l_inv * Ab;
        }
        return make_pair(
                Ab.slice(sq(0, dim - 1, 0, dim - 1)),
                Ab.slice(sq(0, dim - 1, dim, dim)));
    }
};

class GaussColPivotSolver : public GaussBaseSolver {
public:
    GaussColPivotSolver(Mat A, Mat b) : GaussBaseSolver(std::move(A), std::move(b)) {}

    pair<Mat, Mat> Triangularize() override {
        auto Ab = hstack(A, b);
        Mat l_inv;
        for (size_t i = 0; i < dim; ++i) {
            // After line swaping, form a Gauss transforming matrix
            l_inv = Mat::eye(dim);
            size_t max_col_i = i;
            for (size_t j = i + 1; j < dim; ++j)
                max_col_i = fabs(Ab[max_col_i][i]) > fabs(Ab[j][i]) ? max_col_i : j;

            if (max_col_i != i)
                Ab.swap(max_col_i, i);
            for (size_t j = i + 1; j < dim; ++j) l_inv[j][i] = -Ab[j][i] / Ab[i][i];
            Ab = l_inv * Ab;
        }
        return make_pair(
                Ab.slice(sq(0, dim - 1, 0, dim - 1)),
                Ab.slice(sq(0, dim - 1, dim, dim)));
    }
};

class IterBaseSolver : public BaseSolver {
public:
    double precision;
    int max_iter;
    Mat x;

    IterBaseSolver(Mat A, Mat b, Mat x_init = Mat(), double precision = 1e-6, int max_iter = 100)
            : BaseSolver(std::move(A), std::move(b)), x(std::move(x_init)), precision(precision), max_iter(max_iter) {
        if (x.empty())
            x = Mat(A.shape().first, 1, 1);
    }

    virtual pair<Mat, Mat> IterFactor() const = 0;

    virtual Mat Iterate(const pair<Mat, Mat> &iter_factors) {
        using namespace LinearAlgebra;
        if (eig_max(iter_factors.first, x).real() > 1)
            throw S_exception("Iteration method not convergence.");
        double loss = 1e10;
        for (int i = 0; i < max_iter && loss > precision; ++i) {
            x = iter_factors.first * x + iter_factors.second;
            loss = (A * x - b).norm(L2);
        }
        return x;
    }

    Mat solve() override {
        return Iterate(IterFactor());
    }
};

class JacobiIterSolver : public IterBaseSolver {
public:
    JacobiIterSolver(Mat A, Mat b, Mat x_init = Mat(), double precision = 1e-6, int max_iter = 100)
            : IterBaseSolver(std::move(A), std::move(b), std::move(x_init), precision, max_iter) {}

    pair<Mat, Mat> IterFactor() const override {
        using namespace LinearAlgebra;
        auto inv_D = inv(A.diag());
        return make_pair(Mat::eye(dim) - inv_D * A, inv_D * b);
    }
};

class GaussSeidelIterSolver : public IterBaseSolver {
public:
    GaussSeidelIterSolver(Mat A, Mat b, Mat x_init = Mat(), double precision = 1e-6, int max_iter = 100)
            : IterBaseSolver(std::move(A), std::move(b), std::move(x_init), precision, max_iter) {}

    pair<Mat, Mat> IterFactor() const override {
        using namespace LinearAlgebra;
        auto inv_D_L = inv(A.diag() + A.tril(1));
        return make_pair(Mat::eye(dim) - inv_D_L * A, inv_D_L * b);
    }
};

class SORIterSolver : public IterBaseSolver {
public:
    double omega;

    SORIterSolver(Mat _A, Mat b, double omega = 1, Mat x_init = Mat(), double precision = 1e-6, int max_iter = 100)
            : IterBaseSolver(std::move(_A), std::move(b), std::move(x_init), precision, max_iter), omega(omega) {
        if (this->omega == 1 ) {
            if (is_symmetric(A) && is_positive_def(A) && is_3_diag(A)) {
                auto Bj = Mat::eye(dim) - inv(A.diag()) * A;
                auto rho = eig_max(Bj, x).real();
                this->omega = 2 / (1 + sqrt(1 - rho * rho));
            } else {
                cout << "Warning: optimized SOR factor not find. "
                        "Current method is identical to the Gauss Seidel" << endl;
            }
        }
    }

    pair<Mat, Mat> IterFactor() const override {
        using namespace LinearAlgebra;
        auto D = A.diag();
        auto nU = A.triu(1);
        auto inv_D_L_omega = l_inv(D + omega * A.tril(1));
        return make_pair(inv_D_L_omega * ((1-omega) * D - omega * nU), omega * (inv_D_L_omega * b));
    }
};


#endif //NUMALZ_SOLVER_H
