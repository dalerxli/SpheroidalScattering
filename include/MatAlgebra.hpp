#ifndef __MATALGERA__
#define __MATALGERA__

#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_double std::complex<double>
#define lapack_complex_float std::complex<float>

#include <cstddef>
#include <string>
#include <utility>
#include <memory>
#include <fstream>

#include "lapacke.h"
#include <Eigen/Dense>

//#include <boost/python/numpy.hpp>
//namespace p = boost::python;
//namespace np = boost::python::numpy;

template<typename T>
class Matrix {
    public:
    Matrix() = delete;
    Matrix(int m, int n) {
        n_row = m;
        n_col = n;
        data = new T[m*n];
        for(int i = 0; i < m*n; ++i) {
            data[i] = static_cast<T>(0.0);
        }
    }

    Matrix(const Matrix& a) {
        //std::cout << "inside copy constructor" << std::endl;
        n_row = a.GetNRow();
        n_col = a.GetNCol();

        T* a_data = a.GetData();

        data = new T[n_row*n_col];
        for(int i = 0; i < n_row*n_col; ++i) {
            data[i] = a_data[i];
        }
    }

    Matrix& operator=(const Matrix& a) {
        //std::cout << "inside the = operator" << std::endl;
        n_row = a.GetNRow();
        n_col = a.GetNCol();

        T* a_data = a.GetData();

        data = new T[n_row*n_col];
        for(int i = 0; i < n_row*n_col; ++i) {
            data[i] = a_data[i];
        }

        return this;
    }

    ~Matrix() {
        delete[] data;
        data = nullptr;
    }

    T* GetData() const {
        return data;
    }

    int GetNRow() const {
        return n_row;
    }

    int GetNCol() const {
        return n_col;
    }

    T& operator[](const std::pair<int, int> ind) const{
        return data[ind.first*n_col + ind.second];
    }

    T& operator()(const int i, const int j) const {
        return data[i*n_col + j];
    }

    void Reshape(int m, int n) {
        assert(m*n == n_row*n_col);
        n_row = m;
        n_col = n;
    }

    friend Matrix operator-(const Matrix& A, const Matrix& B) {
        double A_M = A.GetNRow();
        double A_N = A.GetNCol();
        double B_M = B.GetNRow();
        double B_N = B.GetNCol();

        assert(A_M == B_M && A_N == B_N);

        Matrix C(A_M, A_N);

        for(int i = 0; i < A_M; ++i) {
            for(int j = 0; j < A_N; ++j) {
                C(i, j) = A[{i, j}] - B[{i, j}];
            }
        }

        return C;
    }

    friend Matrix operator*(const Matrix& A, const Matrix& B) {
        double A_M = A.GetNRow();
        double A_N = A.GetNCol();
        double B_M = B.GetNRow();
        double B_N = B.GetNCol();

        assert( A_N == B_M);

        Matrix C(A_M, B_N);

        for(int i = 0; i < A_M; ++i) {
            for(int j = 0; j < B_N; ++j) {
                for(int k = 0; k < A_N; ++k) {
                    C(i, j) += A[{i, k}] * B[{k, j}];
                }
            }
        }

        return C;
    }

    double GetNorm() {
        double norm = 0.0;
        for(int i = 0; i < n_row*n_col; ++i) {
            norm += std::norm(data[i]);
        }
        return std::sqrt(norm);
    }

    void Print() {
        for(int i = 0; i < n_row; ++i) {
            for(int j = 0; j < n_col; ++j) {
                std::cout << data[i*n_col + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    void WriteToFile(std::string fileName) {
        std::ofstream fileOut(fileName.c_str(), std::ios::out | std::ios::binary);
        assert(fileOut.is_open());
        std::size_t dataSize = n_row * n_col * sizeof(T);
        fileOut.write((char*)(data), dataSize);
        fileOut.close();
    }

    private:
    int n_row;
    int n_col;
    T* data = nullptr;

};


template<typename T>
Matrix<T> ReadMatrixFromFile(std::string fileName) {
    std::ifstream fileIn(fileName.c_str(), std::ios::in | std::ios::binary);
    assert(fileIn.is_open());
    fileIn.seekg(0, fileIn.end);
    std::size_t fileSize = fileIn.tellg();
    assert(fileSize % sizeof(T) == 0);
    int n_col = fileSize / sizeof(T);
    Matrix<T> A(1, n_col);
    T* a_data = A.GetData();
    fileIn.seekg(0, fileIn.beg);
    fileIn.read((char*)(a_data), fileSize);
    fileIn.close();

    return A;
}


Matrix<std::complex<double>> SolveLinear(Matrix<std::complex<double>>& A, Matrix<std::complex<double>>& B,
                                        bool useExpertVersion = false) {

    const int n_row = A.GetNRow();
    const int n_col = A.GetNCol();
    const int n_rhs = B.GetNCol();

    assert(n_rhs == 1);
    assert(n_row == n_col);

    std::complex<double>* a = A.GetData();
    std::complex<double>* b = B.GetData();
    std::unique_ptr<int[]> ipiv(new int[n_row]);

    if(useExpertVersion) {
        std::unique_ptr<double[]> rpivot(new double [n_row]);

        char equed = 'A';

        std::complex<double> af[n_row][n_col];
        std::unique_ptr<double[]> r(new double[n_row]);    // row scale factors
        std::unique_ptr<double[]> c(new double[n_row]);    // col scale factors
        Matrix<std::complex<double>> X(n_row, n_rhs);
        std::complex<double>* x = X.GetData();

        double rcond;
        std::unique_ptr<double[]> ferr(new double[n_rhs]);
        std::unique_ptr<double[]> berr(new double[n_rhs]);

        lapack_int info = LAPACKE_zgesvx(LAPACK_ROW_MAJOR, 'E', 'N',
                               n_row, n_rhs, a, n_col, *af, n_row,
                               ipiv.get(), &equed, r.get(), c.get(),
                               b, n_rhs, x, n_rhs,
                               &rcond, ferr.get(), berr.get(),
                               rpivot.get() );


        std::cout << " info : " << info << std::endl;
        std::cout << " condition : " << rcond << std::endl;
        std::cout << " ferr : " << ferr[0] << std::endl;
        std::cout << " berr : " << berr[0] << std::endl;
        std::cout << " equed : " << equed << std::endl;

        /*
        std::cout << " ipiv : " << std::endl;
        int* ipiv_p = ipiv.get();
        for(int i = 0; i < n_row; ++i) {
            std::cout << "(" << i + 1 << ", " << ipiv_p[i] << ")  ";
        }
        std::cout << std::endl;*/

        return X;
    } else {
        lapack_int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n_row, n_rhs, a, n_col, ipiv.get(), b, n_rhs);

        Matrix<std::complex<double>> X(n_row, n_rhs);
        std::complex<double>* x = X.GetData();
        for(int i = 0; i < n_row*n_rhs; ++i) {
            x[i] = b[i];
        }

        return X;
    }
}

Matrix<std::complex<double>> SolveLinear_Ei(Matrix<std::complex<double>>& A, Matrix<std::complex<double>>& B) {
    const int n_row = A.GetNRow();
    const int n_col = A.GetNCol();
    const int n_rhs = B.GetNCol();

    assert(n_rhs == 1);
    assert(n_row == n_col);

    Eigen::MatrixXcd A_ei(n_row, n_col);
    Eigen::MatrixXcd B_ei(n_row, n_rhs);

    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_col; ++j) {
            A_ei(i, j) = A(i, j);
        }
        for(int j = 0; j < n_rhs; ++j) {
            B_ei(i, j) = B(i, j);
        }
    }

    //std::cout << A_ei << std::endl;
    Matrix<std::complex<double>> X(n_row, n_rhs);
    Eigen::MatrixXcd X_ei = A_ei.fullPivLu().solve(B_ei);
    std::cout << X_ei << std::endl;

    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_rhs; ++j) {
            X(i, j) = X_ei(i, j);
        }
    }

    double relative_error = (A_ei*X_ei - B_ei).norm() / B_ei.norm(); // norm() is L2 norm
    std::cout << "The relative error is:" << relative_error << std::endl;
    relative_error = (A*X - B).GetNorm() / B.GetNorm(); // norm() is L2 norm
    std::cout << "The relative error is:" << relative_error << std::endl;

    return X;
}

Matrix<std::complex<double>> SolveLinear_np_file(Matrix<std::complex<double>>& A, Matrix<std::complex<double>>& B) {
    remove("out/x.data");
    A.WriteToFile("out/A.data");
    B.WriteToFile("out/b.data");
    system("python3 PyMatSolver.py");

    auto x = ReadMatrixFromFile<std::complex<double>>("out/x.data");
    int m = x.GetNRow();
    int n = x.GetNCol();
    assert(m == 1);
    x.Reshape(m*n, 1);
    //x.Print();
    return x;
}


#include "suitesparse/umfpack.h"

Matrix<std::complex<double>> SolveLinear_umfpack(Matrix<std::complex<double>>& A, Matrix<std::complex<double>>& B) {
    const int n = A.GetNRow();
    assert(A.GetNCol() == n && B.GetNRow() == n && B.GetNCol() == 1);

    int nnz = 0;
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            if( A(i, j) != 0.0 ) {
                nnz += 1;
            }
        }
    }

    std::cout << "n: " << n << "  nnz: " << nnz << std::endl;

    auto a_cols = std::make_unique<int[]>(n + 1);
    auto a_rows = std::make_unique<int[]>(nnz);
    auto a_elems = std::make_unique<double[]>(2*nnz);
    auto b_elems = std::make_unique<double[]>(2*n);
    auto x_elems = std::make_unique<double[]>(2*n);

    int* Ap = a_cols.get();
    int* Ai = a_rows.get();
    double* Ax = a_elems.get();
    double* b = b_elems.get();
    double* x = x_elems.get();

    nnz = 0;
    Ap[0] = 0;
    for(int j = 0; j < n; ++j) {
        for(int i = 0; i < n; ++i) {
            if( A(i, j) != 0.0 ) {
                Ai[nnz] = i;
                Ax[2*nnz]     = A(i, j).real();
                Ax[2*nnz + 1] = A(i, j).imag();
                nnz += 1;
            }
        }
        Ap[j + 1] = nnz;
    }

    for(int i = 0; i < n; ++i) {
        b[2*i]     = B(i, 0).real();
        b[2*i + 1] = B(i, 0).imag();
    }

    double *null = (double *) NULL;
    void *Symbolic, *Numeric;
    (void) umfpack_zi_symbolic(n, n, Ap, Ai, Ax, NULL, &Symbolic, null, null);
    (void) umfpack_zi_numeric(Ap, Ai, Ax, NULL, Symbolic, &Numeric, null, null);
    umfpack_zi_free_symbolic(&Symbolic);
    (void) umfpack_zi_solve (UMFPACK_A, Ap, Ai, Ax, NULL, x, NULL, b, NULL, Numeric, null, null);
    umfpack_zi_free_numeric (&Numeric);

    Matrix<std::complex<double>> X(n, 1);
    //for (int i = 0 ; i < n ; i++) std::cout << i << " " <<  x_complex[i] << std::endl;

    for(int i = 0; i < n; ++i) {
        X(i, 0) = std::complex<double>(x[2*i], x[2*i + 1]);
    }

    return X;
}

#endif

