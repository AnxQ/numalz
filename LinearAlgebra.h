//
// Created by 安啸琪 on 2019/12/22.
//

#ifndef NUMALZ_LINEARALGEBRA_H
#define NUMALZ_LINEARALGEBRA_H

#include <numeric>
#include <complex>
#include <utility>
#include <vector>
#include <string>
#include "Matrix.h"

namespace LinearAlgebra {

    enum QRMethod {
        Givens,
        Householder,
        GramSchmidt
    };

    enum InvMethod {
        QR,
        LU,
        Adjugate
    };

    enum LUMethod {
        Doolittle,
    };

    template<typename _T>
    _T det(const Matrix<_T> &);

    class LA_exception : exception {
        string msg;
    public:
        explicit LA_exception(string msg) : msg(std::move(msg)) {}

        const char *what() { return msg.c_str(); }
    };

    template<typename _T>
    bool is_symmetric(const Matrix<_T> &m) { return m.transpose() == m; }

    template<typename _T>
    bool is_positive_def(const Matrix<_T> &m) {
        for (int i = 0; i < min(m.shape().first, m.shape().second); ++i)
            if (det(m.slice(sq(0, i, 0, i))) < 0)
                return false;
        return true;
    }

    template<typename _T>
    bool is_3_diag(const Matrix<_T> &m, const double acc = 0) {
        for (size_t i = 0; i < m.shape().first; ++i) {
            for (size_t _i = 0; _i < i - 1; ++_i)
                if (fabs(m[i][_i]) > acc) return false;
            for (size_t _i = i + 2; _i < m.shape().second; ++_i)
                if (fabs(m[i][_i]) > acc) return false;
        }
        return true;
    }

    template<typename _T>
    static Matrix<_T> householder(Matrix<_T> x, int bias = 0) {
        if (bias < 0 || bias > x.shape().first)
            throw LA_exception("Bias must be in range.");
        if (x.shape().second > 1)
            throw LA_exception("Must specific a vector for generating Household transform matrix.");
        auto dim = x.shape().first;
        auto trans_x = x.slice(make_pair(bias, dim - 1), make_pair(0, 0));
        auto e = Matrix<_T>(dim - bias, 1);
        e[0][0] = _T(1);
        _T k = -sgn(trans_x[0][0]) * trans_x.norm(NormType::L2);
        auto u = trans_x - k * e;
        auto pre = 2 * u * u.transpose() / pow(u.norm(NormType::L2), 2.0);
        auto H = Matrix<_T>::eye(x.shape().first);
        for (size_t i = bias; i < dim; ++i) {
            for (size_t j = bias; j < dim; ++j)
                H[i][j] = H[i][j] - pre[i - bias][j - bias];
        }
        return H;
    }

    template<typename _T>
    static Matrix<_T> givens(int i, int k, const Matrix<_T> &x) {
        if (x.shape().second > 1)
            throw LA_exception("Must specific a vector for generating Household transform matrix.");
        if (k < i || i < 0 || k > x.shape().first)
            throw LA_exception("i and k must be in range.");
        auto trans = Matrix<_T>::eye(x.shape().first);
        double c = 0, s = 0, t = 0;
        if (x[k][0] == 0)
            c = 1, s = 0;
        else if (fabs(x[k][0]) >= fabs(x[i][0])) {
            t = x[i][0] / x[k][0];
            s = sgn(x[k][0]) * pow(1 + t * t, -0.5);
            c = s * t;
        } else {
            t = x[k][0] / x[i][0];
            c = sgn(x[i][0]) * pow(1 + t * t, -0.5);
            s = c * t;
        }
        trans[i][i] = c;
        trans[i][k] = s;
        trans[k][i] = -s;
        trans[k][k] = c;
        return trans;
    }

    template<typename _T>
    pair<Matrix<_T>, Matrix<_T>> qr(const Matrix<_T> &M, int method = QRMethod::Householder) {
        M.assert_square("qr");
        size_t m = M.shape().first;
        auto R = M;
        auto Q = Matrix<_T>::eye(m);
        if (method == QRMethod::Householder) {
            for (size_t i = 0; i < m - 1; ++i) {
                auto H = householder(R.slice(sq(0, m - 1, i, i)), i);
                R = H * R;
                Q = H * Q;
            }
            return make_pair(Q.transpose(), R);
        } else if (method == QRMethod::Givens) {
            for (size_t i = 0; i < m - 1; ++i) {
                auto G = Matrix<_T>::eye(m);
                auto Vi = R.slice(sq(0, m - 1, i, i));
                for (size_t k = i + 1; k < m; ++k) {
                    G = givens(i, k, Vi) * G;
                    Vi = G * Vi;
                }
                R = G * R;
                Q = G * Q;
            }
            return make_pair(Q.transpose(), R);
        } else if (method == QRMethod::GramSchmidt) {
            // TODO: pending implement.
        }
        throw LA_exception("QR method unknown");
    }

    template<typename _T>
    pair<Matrix<_T>, Matrix<_T>> lu(const Matrix<_T> &M) {
        auto dim = M.shape().first;
        Mat L = Mat::eye(dim), U = Mat(dim, dim);
        U[0] = M[0];
        for (size_t i = 1; i < dim; ++i) L[i][0] = M[i][0] / U[0][0];
        for (size_t i = 1; i < dim; ++i) {
            for (size_t j = i; j < dim; ++j)
                U[i][j] = M[i][j] - (L.slice({i}, Row) * U.slice({j}, Col))[0][0];
            for (size_t j = i + 1; j < dim; ++j)
                L[j][i] = (M[j][i] - (L.slice({j}, Row) * U.slice({i}, Col))[0][0]) / U[i][i];
        }
        return make_pair(L, U);
    }

    template<typename _T>
    _T _det_recursive(const Matrix<_T> &M) {
        size_t m = M.shape().first;
        if (m == 1)
            return M[0][0];
        Matrix<_T> tmp_l = M.slice(sq(1, m - 1, 1, m - 1)),
                tmp_r = M.slice(sq(1, m - 1, 0, m - 2));
        _T sum = M[0][0] * _det_recursive(tmp_l) +
                 pow(-1, m - 1) * M[0][m - 1] * _det_recursive(tmp_r);
        for (size_t i = 1; i < m - 1; ++i) {
            Matrix<_T> tmp = hstack(M.slice(sq(1, m - 1, 0, i - 1)),
                                    M.slice(sq(1, m - 1, i + 1, m - 1)));
            sum += pow(-1, i) * M[0][i] * _det_recursive(tmp);
        }
        return sum;
    }

    template<typename _T>
    _T det(const Matrix<_T> &M) {
        M.assert_square("det");
        return _det_recursive(M);
    }

    // Inverse an upper triangular matrix, rolling calculation
    template<typename _T>
    Matrix<_T> u_inv(const Matrix<_T> &M) {
        Matrix<_T> res(M.shape());
        M.assert_square("u_inv");
        for (int k = 0; k < M.shape().first; ++k) {
            for (int i = 0, j = k; j < M.shape().first; ++i, ++j) {
                _T sum = 0;
                for (int m = i + 1; m <= j; ++m) sum += res[m][j] * M[i][m];
                res[i][j] = (i == j ? 1 : -sum) / M[i][i];
            }
        }
        return res;
    }

    template<typename _T>
    Matrix<_T> l_inv(const Matrix<_T> &M) {
        return u_inv(M.transpose()).transpose();
    }

    // Inverse any non-singular square matrix
    template<typename _T>
    Matrix<_T> inv(const Matrix<_T> &M, int method = InvMethod::QR) {
        M.assert_square("inv");
        auto dim = M.shape().first;
        _T m_det = det(M);
        if (m_det == _T(0))
            throw LA_exception("Calling inv on a singular matrix.");
        Matrix<_T> res;
        if (method == InvMethod::QR) {
            // A = QR, A^-1 = R^-1*Q^-1 = R^-1*Q'
            auto _qr = qr(M);
            res = u_inv(_qr.second) * (_qr.first.transpose());
        } else if (method == InvMethod::LU) {
            // A = LU, A^-1 = U^-1 * (L'^-1)'
            auto _lu = lu(M);
            auto lt = _lu.second.transpose();
            res = u_inv(lt) * u_inv(_lu.first);
        } else if (method == InvMethod::Adjugate) {
            // Be careful: very slow.
            Matrix<_T> M_adj(dim, dim);
            Matrix<_T> M_cof;
            for (size_t i = 0; i < dim; ++i)
                for (size_t j = 0; j < dim; ++j) {
                    M_cof = M.cofactor(i, j);
                    M_adj[j][i] = det(M_cof);
                }
            res = std::move(M_adj / m_det);
        } else {
            throw LA_exception("Inv method unknown");
        }
        return res;
    }

    template<typename _T>
    class Eigen : public complex<_T> {
    public:
        Eigen() = default;

        Eigen(_T re, _T im) : complex<_T>(re, im) {}

        friend std::ostream &operator<<(std::ostream &out, const Eigen<_T> &eigen) {
            if (fabs(eigen.imag()) < 1e-10)
                out << eigen.real();
            else
                out << eigen.real() << (eigen.imag() < 0 ? "" : "+") << eigen.imag() << "i";
            return out;
        }
    };

    typedef Eigen<double> Eig;

    template<typename _T>
    auto const abs_max = [](const _T &a, const _T &b) { return fabs(a) < fabs(b); };

    // Calculate min eigen val by power iteration
    template<typename _T>
    Eigen<_T>
    eig_max(const Matrix<_T> &M, const Matrix<_T> &init,
            const int &max_iter = 100, const double &acc = 1e-10) {
        // Power iteration
        auto v = init;
        vector<_T> v_col;
        _T k = 1e100, k_pre;
        for (int i = 0; i < max_iter && fabs(k - k_pre) > acc; ++i) {
            v = M * v;
            v_col = v.col();
            k_pre = k;
            k = *max_element(v_col.begin(), v_col.end(), abs_max<_T>);
            v = v / k;
        }
        return Eigen<_T>(k, 0);
    }

    // Calculate min eigen val by inverting power iteration
    template<typename _T>
    Eigen<_T>
    eig_min(const Matrix<_T> &M, const Matrix<_T> &init,
            const bool balance = false,
            const int &max_iter = 100,
            const double &acc = 1e-10) {
        auto v = init;
        auto dim = M.shape().first;
        Matrix<_T> D = Mat::eye(dim), DM = M;
        vector<_T> v_col, v_row, S(dim);
        _T k = 1e100, k_pre;
        if (balance) {
            for (size_t j = 0; j < dim; ++j)
                v_row = M.row(j), S[j] = 1 / fabs(*max_element(v_row.begin(), v_row.end(), abs_max<_T>));
            D = Mat::diag(S);
            DM = D * M;
        }
        for (int i = 0; i < max_iter && fabs(k - k_pre) > acc; ++i) {
            v = u_inv(DM) * D * v;
            v_col = v.col();
            k_pre = k;
            k = *max_element(v_col.begin(), v_col.end(), abs_max<_T>);
            v = v / k;
        }
        return Eigen<_T>(1 / k, 0);
    }

    // Calculate all eigen val by QR Hessenberg iteration
    template<typename _T>
    vector<Eigen<_T>> eig(const Matrix<_T> &M, int max_iter = 100) {
        M.assert_square("eig");
        vector<Eigen<_T>> eigens;
        auto n = M.shape().first;
        // QR iterating
        auto _qr = qr(M);
        auto hess = _qr.second * _qr.first;
        for (int i = 0; i < max_iter; ++i) {
            _qr = qr(hess);
            hess = _qr.second * _qr.first;
        }
        cout << "Hessenberg Matrix: " << endl << hess;
        // Do not ignore the none-convergence value, make it complex.
        for (size_t i = 0; i < M.shape().first - 1; ++i) {
            if (fabs(hess[i + 1][i]) < 1e-10)
                eigens.emplace_back(Eigen<_T>(hess[i][i], 0));
            else {
                _T b = -(hess[i][i] + hess[i + 1][i + 1]), c = (hess[i][i] * hess[i + 1][i + 1] -
                                                                hess[i][i + 1] * hess[i + 1][i]),
                        delta = b * b - 4 * c;
                if (delta > 0) {
                    eigens.emplace_back(Eigen<_T>(-b / 2 + sqrt(delta) / 2, 0));
                    eigens.emplace_back(Eigen<_T>(-b / 2 - sqrt(delta) / 2, 0));
                } else {
                    eigens.emplace_back(Eigen<_T>(-b / 2, sqrt(-delta) / 2));
                    eigens.emplace_back(Eigen<_T>(-b / 2, -sqrt(-delta) / 2));
                }
                i++;
            }
        }
        if (eigens.size() < n)
            eigens.emplace_back(Eigen<_T>(hess[n - 1][n - 1], 0));
        return eigens;
    }
}


#endif //NUMALZ_LINEARALGEBRA_H
