//
// Created by 安啸琪 on 2019/12/12.
//

#ifndef NUMALZ_MATRIX_H
#define NUMALZ_MATRIX_H

#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>

using namespace std;

#define sq(m, m_, n, n_) make_pair(m, m_), make_pair(n, n_)
enum NormType {
    L1, L2, Inf, Forbenius
};
enum AxisType {
    Row, Col
};

template<typename _T>
inline bool sgn(const _T &x) {
    if (fabs((double) x) < 1e-6) return _T(1);
    return x < _T(0) ? _T(-1) : _T(1);
}

class M_exception : exception {
    string msg;
public:
    explicit M_exception(string msg) : msg(std::move(msg)) {}

    const char *what() { return msg.c_str(); }
};

template<typename _T = double>
class Matrix : public vector<vector<_T>> {
private:
    size_t m = 0, n = 0;
public:
    Matrix() = default;

    Matrix(const vector<vector<_T>> &in) : vector<vector<_T>>(in) {
        m = vector<vector<_T>>::size();
        n = m > 0 ? (*this)[0].size() : 0;
    }

    // When input a vector, treat as vertical
    explicit Matrix(const vector<_T> &in) {
        vector<vector<_T>>::resize(m = in.size());
        for (size_t i = 0; i < in.size(); ++i) {
            (*this)[i].resize(n = 1, in[i]);
        }
    }

    Matrix(int m, int n, _T init) : m(m), n(n) { vector<vector<_T>>::resize(m, vector<_T>(n, init)); }

    Matrix(int m, int n) : Matrix(m, n, _T(0)) {}

    Matrix(pair<size_t, size_t> s, _T init) : Matrix(s.first, s.second, _T(0)) {}

    explicit Matrix(pair<size_t, size_t> s) : Matrix(s.first, s.second) {}

    pair<size_t, size_t> shape() const { return make_pair(m, n); }

    Matrix transpose() const {
        Matrix res(n, m);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                res[j][i] = (*this)[i][j];
        return res;
    }

    Matrix for_each(function<_T(_T)> func) {
        Matrix res(m, n);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                res[i][j] = func((*this)[i][j]);
        return res;
    }

    Matrix abs() {
        if (typeid(_T) == typeid(int)) {
            return for_each([](_T i) {
                return _T(std::abs((int) i));
            });
        } else {
            return for_each([](_T i) {
                return _T(std::fabs((float) i));
            });
        }
    }

    void de_zero(double threshold = 1e-6) {
        // Accuracy problem: if abs less than 1e-6, make it zero
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                if (fabs((double) (*this)[i][j]) < threshold)
                    (*this)[i][j] = 0;
    }

    friend Matrix operator+(const Matrix &L, const Matrix &R) {
        auto shape = L.shape();
        if (shape != R.shape())
            throw M_exception("Matrix size mismatch");
        Matrix<_T> res = Matrix(shape);
        for (size_t i = 0; i < shape.first; ++i)
            for (size_t j = 0; j < shape.second; ++j)
                res[i][j] = L[i][j] + R[i][j];
        return res;
    }

    friend Matrix operator-(const Matrix &L, const Matrix &R) {
        auto shape = L.shape();
        if (shape != R.shape())
            throw M_exception("Matrix size mismatch");
        Matrix<_T> res = Matrix(shape);
        for (size_t i = 0; i < shape.first; ++i)
            for (size_t j = 0; j < shape.second; ++j)
                res[i][j] = L[i][j] - R[i][j];
        return res;
    }

    friend Matrix operator*(const Matrix &L, const Matrix &R) {
        auto l_shape = L.shape(), r_shape = R.shape();
        _T temp;
        if (l_shape.second != r_shape.first)
            throw M_exception("Inner dimension mismatch");
        Matrix<_T> res = Matrix(l_shape.first, r_shape.second);
        for (size_t i = 0; i < l_shape.first; ++i)
            for (size_t j = 0; j < r_shape.second; ++j) {
                temp = 0;
                for (size_t k = 0; k < l_shape.second; ++k) temp += L[i][k] * R[k][j];
                res[i][j] = temp;
            }
//        res.de_zero();
        return res;
    }

    friend Matrix operator*(const Matrix &L, const _T &r) {
        auto shape = L.shape();
        Matrix<_T> res = Matrix(shape);
        for (size_t i = 0; i < shape.first; ++i)
            for (size_t j = 0; j < shape.second; ++j)
                res[i][j] = L[i][j] * r;
//        res.de_zero();
        return res;
    }

    friend Matrix operator*(const _T &l, const Matrix &R) {
        return R * l;
    }

    friend Matrix operator/(const Matrix &L, const _T &r) {
        auto shape = L.shape();
        Matrix<_T> res = Matrix(shape);
        for (size_t i = 0; i < shape.first; ++i)
            for (size_t j = 0; j < shape.second; ++j)
                res[i][j] = L[i][j] / r;
//        res.de_zero();
        return res;
    }

    friend bool operator==(const Matrix &L, const Matrix &R) {
        if (L.shape() != R.shape())
            return false;
        for (size_t i = 0; i < L.shape().first; ++i)
            if (L[i] != R[i])
                return false;
        return true;
    }

    friend std::ostream &operator<<(std::ostream &out, const Matrix<_T> &m) {
        auto f = out.flags();
        int wide = 5;
        // If a number have noticeable decimal part, adjust the output wide to fit.
        for (auto &row: m) {
            for (auto &item: row)
                if (std::abs(item - (int) item) > 1e-5) {
                    out.precision(5), wide = 11;
                    break;
                }
            if (wide == 11) break;
        }
        out << "[";
        for (auto iter = m.begin(); iter != m.end(); ++iter) {
            out << (iter == m.begin() ? "[" : " [");
            for (auto item: *iter) {
                out.width(wide);
                out << item;
            }
            out << (iter == (m.end() - 1) ? "]" : "]\n");
        }
        out << "]" << endl;
        out.flags(f);
        return out;
    }



    // Slice into specific row vector / column vector
    Matrix slice(const vector<size_t> &indices, int axis = AxisType::Col) const {
        if (axis == AxisType::Row) {
            size_t _m = indices.size();
            Matrix<_T> res = Matrix(_m, n);
            for (size_t i = 0; i < indices.size(); ++i)
                res[i] = (*this)[indices[i]];
            return res;
        } else if (axis == AxisType::Col) {
            size_t _n = indices.size();
            Matrix<_T> res = Matrix(m, _n);
            for (size_t i = 0; i < m; ++i)
                for (size_t j = 0; j < indices.size(); ++j)
                    res[i][j] = (*this)[i][indices[j]];
            return res;
        }
        return Matrix<_T>();
    }

    // Slice a part of the matrix
    Matrix slice(const pair<size_t, size_t> m_range, pair<size_t, size_t> n_range) const {
        auto check_range = [](pair<size_t, size_t> range) {
            return range.first <= range.second && range.first >= 0;
        };
        if (!check_range(m_range) || !check_range(n_range) || m_range.second >= m || n_range.second >= n)
            throw M_exception("Slice indices must be in range.");
        Matrix res = Matrix(m_range.second - m_range.first + 1, n_range.second - n_range.first + 1);
        for (auto i = m_range.first; i <= m_range.second; ++i)
            for (auto j = n_range.first; j <= n_range.second; ++j)
                res[i - m_range.first][j - n_range.first] = (*this)[i][j];
        return res;
    }

    void assert_square(const string &invoker) const {
        if (m == n)
            return;
        else
            throw M_exception("Attempt to call " + invoker + " on non-square matrix.");
    }

    Matrix diag() const {
        Matrix<_T> res = Matrix<_T>(this->shape());
        for (size_t i = 0; i < m; ++i)
            res[i][i] = (*this)[i][i];
        return res;
    }

    Matrix tril(int bias = 0) const {
        Matrix<_T> res = Matrix<_T>(this->shape());
        if (std::abs(bias) >= std::max(m, n))
            return res;
        for (int i = min(bias, 0); i < m; ++i)
            for (int j = 0; j < min(i + 1 - bias, (int) n); ++j)
                res[i][j] = (*this)[i][j];
        return res;
    }

    Matrix triu(int bias) const {
        Matrix<_T> res = Matrix<_T>(this->shape());
        if (std::abs(bias) >= std::max(m, n))
            return res;
        for (int i = min(m - 1, m - bias - 1); i >= 0; --i)
            for (int j = n - 1; j >= std::max(i + bias, 0); --j)
                res[i][j] = (*this)[i][j];
        return res;
    }

    vector<_T> row(const size_t &i = 0) const {
        if (i >= m)
            throw M_exception("Row index out of range.");
        return (*this)[i];
    }

    vector<_T> col(const size_t &j = 0) const {
        if (j >= n)
            throw M_exception("Col index out of range.");
        return this->transpose()[j];
    }

    void swap(int a, int b, int axis = AxisType::Row) {
        if (axis == AxisType::Row) {
            if (a < 0 || b < 0 || a > m || b > m)
                throw M_exception("Row index out of range.");
            for (size_t j = 0; j < n; ++j)
                std::swap((*this)[a][j], (*this)[b][j]);
        } else if (axis == AxisType::Col) {
            if (a < 0 || b < 0 || a > n || b > n)
                throw M_exception("Col index out of range.");
            for (size_t i = 0; i < m; ++i)
                std::swap((*this)[i][a], (*this)[i][b]);
        } else {
            throw M_exception("Axis unkown");
        }
    }

    static Matrix<_T> ones(size_t m, size_t n) {
        return Matrix<_T>(m, n, 1);
    }

    static Matrix<_T> zeros(size_t m, size_t n) {
        return Matrix<_T>(m, n, 0);
    }

    static Matrix<_T> diag(const vector<_T> &d_vec) {
        size_t dim = d_vec.size();
        Matrix<_T> res = Matrix<_T>(dim, dim);
        for (size_t i = 0; i < dim; ++i)
            res[i][i] = d_vec[i];
        return res;
    }

    static Matrix<_T> eye(size_t m) {
        return diag(vector<_T>(m, (_T) 1.0f));
    }

    _T trace() {
        assert_square("trace");
        _T tr;
        for (size_t i = 0; i < m; ++i)
            tr += (*this)[i][i];
        return tr;
    }

    Matrix cofactor(size_t i, size_t j) const {
        Matrix res(m - 1, n - 1);
        for (size_t k = 0; k < m; ++k) {
            if (k == i) continue;
            for (size_t l = 0; l < n; ++l) {
                if (l == j) continue;
                res[k < i ? k : k - 1][l < j ? l : l - 1] = (*this)[k][l];
            }
        }
        return pow(-1, i + j) * res;
    }

    double norm(int type = NormType::Forbenius) {
        if (!(m == 1 || n == 1))
            assert_square("norm");
        else if (type == Forbenius)
            throw M_exception("Calling Forbenius norm on vector.");
        size_t i = 0;
        auto m_abs = this->abs();
        _T sum = 0;
        vector<_T> vec_norm(n);
        switch (type) {
            case NormType::L1:
                for (vector<_T> &col: m_abs.transpose()) {
                    vec_norm[i++] = accumulate(col.begin(), col.end(), _T(0));
                }
                return *max_element(vec_norm.begin(), vec_norm.end());
            case NormType::Inf:
                for (vector<_T> &col: m_abs) {
                    vec_norm[i++] = accumulate(col.begin(), col.end(), _T(0));
                }
                return *max_element(vec_norm.begin(), vec_norm.end());
            case NormType::Forbenius:
                for (vector<_T> &col: m_abs)
                    for (_T &item: col)
                        sum += item * item;
                return sqrt(sum);
            case NormType::L2:
                if (m == 1 || n == 1) {
                    for (vector<_T> &col: (m == 1 ? m_abs.transpose() : m_abs))
                        sum += col[0] * col[0];
                    return sqrt(sum);
                } else {
                    // TODO: Not implemented
                    return 0;
                }
            default:
                break;
        }
        return _T(0);
    }
};

template<typename _T>
Matrix<_T> vstack(const Matrix<_T> &U, const Matrix<_T> &D) {
    if (U.shape().second != D.shape().second)
        throw M_exception("Matrix horizon size mismatch.");
    size_t u_m = U.shape().first, d_m = D.shape().first;
    size_t m = u_m + d_m, n = D.shape().second;
    Matrix<_T> res(m, n);
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            res[i][j] = i < u_m ? U[i][j] : D[i - u_m][j];
    return res;
}

template<typename _T>
Matrix<_T> hstack(const Matrix<_T> &L, const Matrix<_T> &R) {
    if (L.shape().first != R.shape().first)
        throw M_exception("Matrix horizon size mismatch.");
    size_t l_n = L.shape().second, r_n = R.shape().second;
    size_t m = L.shape().first, n = l_n + r_n;
    Matrix<_T> res(m, n);
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            res[i][j] = j < l_n ? L[i][j] : R[i][j - l_n];
    return res;
}

typedef Matrix<double> Mat;

#endif //NUMALZ_MATRIX_H
