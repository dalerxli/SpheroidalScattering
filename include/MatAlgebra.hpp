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

#endif

