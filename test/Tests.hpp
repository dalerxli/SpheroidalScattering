
#ifndef __TESTS__
#define __TESTS__

#include <utility>
#include <string>


#include "boost/lexical_cast.hpp"

#include <sys/types.h>  // directory check
#include <sys/stat.h>

#include <cstdlib>

#include "MatAlgebra.hpp"


void TestExpansionCoefficients(int m, int n, double c) {
    double dmn[200];
    int dmn_inds[200];
    int len_dmn = 0;
    double thresh = 1.0e-15;
    GetExpansionCoefficientsDmn(m, n, c, dmn, dmn_inds, &len_dmn, thresh);

    for(int i = 0; i < len_dmn; ++i) {
        std::cout << "i:" << i << "  ind: " << dmn_inds[i] <<  "  d_mn[i]: " << dmn[i] << std::endl;
    }
}

void TestRadialProlate(int m, int n, double c, double x) {
    double cv;
    double eg[200];

    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 3;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d);

    std::cout << "Prolate functions and derivatives: " << std::endl;
    std::cout << "r1f: " << r1f << "  r1d: " << r1d << "  r2f: " << r2f << "  r2d: " << r2d << std::endl;
}

void TestRadialOblate(int m, int n, double c, double x) {
    double cv;
    double eg[200];

    // get characteristic value
    int kd = -1; // oblate
    segv_(&m, &n, &c, &kd, &cv, eg);

    // get prolate radials and derivatives
    int kf = 3;
    double r1f, r1d, r2f, r2d;
    rswfo_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d);

    std::cout << "Oblate functions and derivatives: " << std::endl;
    std::cout << "r1f: " << r1f << "  r1d: " << r1d << "  r2f: " << r2f << "  r2d: " << r2d << std::endl;
}

void TestAngularFirstKind(int m, int n, double c, double x) {
    double cv;
    double eg[200];

    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);

    // get prolate radial and derivatives
    double s1f, s1d;
    aswfa_(&m, &n, &c, &x, &kd, &cv, &s1f, &s1d);

    std::cout << "Prolate angular function and derivative: " << std::endl;
    std::cout << "s1f: " << s1f << "  s1d: " << s1d << std::endl;
}

lapack_int TestLapack_dgels() {
   double a[5][3] = {1,1,1,2,3,4,3,5,2,4,2,5,5,4,3};
   double b[5][2] = {-10,-3,12,14,14,12,16,16,18,16};
   lapack_int info,m,n,lda,ldb,nrhs;
   int i,j;

   m = 5;
   n = 3;
   nrhs = 2;
   lda = 3;
   ldb = 2;

   info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,*a,lda,*b,ldb);


   for(i=0;i<n;i++)
   {
      for(j=0;j<nrhs;j++)
      {
         printf("%lf ",b[i][j]);
      }
      printf("\n");
   }
   return(info);
}

lapack_int TestLapack_dgesv() {
    srand (time(NULL));

    const int n_row = 5;
    const int n_col = n_row;
    const int n_rhs = 1;

    double a[n_row][n_col];
    double a0[n_row][n_col];
    double b[n_row][n_rhs];
    double b0[n_row][n_rhs];
    int ipiv[n_row];

    std::cout << "A : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_col; ++j) {
            a[i][j] = (double)std::rand() / RAND_MAX;
            a0[i][j] = a[i][j];
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_rhs; ++j) {
            b[i][j] = (double)std::rand() / RAND_MAX;
            b0[i][j] = b[i][j];
            std::cout << b[i][j] << " ";
        }
        std::cout << std::endl;
    }

    lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n_row, n_rhs, *a, n_col, ipiv, *b, n_rhs);

    std::cout << "A : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_col; ++j) {
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_rhs; ++j) {
            std::cout << b[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Ax - B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        double ax_b_i = 0.0;
        for(int j = 0; j < n_col; ++j) {
            ax_b_i += a0[i][j] * b[j][0];
        }
        std::cout << ax_b_i - b0[i][0] << std::endl;
    }

    return info;
}

lapack_int TestLapack_dgesvx() {
    srand (time(NULL));

    const int n_row = 5;
    const int n_col = n_row;
    const int n_rhs = 1;

    double a[n_row][n_col];
    double a0[n_row][n_col];
    double b[n_row][n_rhs];
    double b0[n_row][n_rhs];
    int ipiv[n_row];
    double rpivot[n_row];

    char equed = 'N';

    double af[n_row][n_col];
    double r[n_row];    // row scale factors
    double c[n_row];    // col scale factors
    double x[n_row][n_rhs];

    double rcond;
    double ferr[n_rhs];
    double berr[n_rhs];

    std::cout << "A : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_col; ++j) {
            a[i][j] = (double)std::rand() / RAND_MAX;
            a0[i][j] = a[i][j];
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_rhs; ++j) {
            b[i][j] = (double)std::rand() / RAND_MAX;
            b0[i][j] = b[i][j];
            std::cout << b[i][j] << " ";
        }
        std::cout << std::endl;
    }

    lapack_int info = LAPACKE_dgesvx(LAPACK_ROW_MAJOR, 'E', 'N',
                           n_row, n_rhs, *a, n_col, *af, n_row,
                           ipiv, &equed, r, c,
                           *b, n_rhs, *x, n_rhs,
                           &rcond, ferr, berr,
                           rpivot );

    std::cout << "A : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_col; ++j) {
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_rhs; ++j) {
            std::cout << b[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Ax - B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        double ax_b_i = 0.0;
        for(int j = 0; j < n_col; ++j) {
            ax_b_i += a0[i][j] * x[j][0];
        }
        std::cout << ax_b_i - b0[i][0] << std::endl;
    }

    std::cout << " condition : " << rcond << std::endl;
    std::cout << " ferr : " << ferr[0] << std::endl;
    std::cout << " berr : " << berr[0] << std::endl;

    return info;
}


lapack_int TestLapack_zgesvx() {
    srand (time(NULL));

    const int n_row = 5;
    const int n_col = n_row;
    const int n_rhs = 1;

    lapack_complex_double a[n_row][n_col];
    lapack_complex_double a0[n_row][n_col];
    lapack_complex_double b[n_row][n_rhs];
    lapack_complex_double b0[n_row][n_rhs];
    int ipiv[n_row];
    double rpivot[n_row];

    char equed = 'N';

    lapack_complex_double af[n_row][n_col];
    double r[n_row];    // row scale factors
    double c[n_row];    // col scale factors
    lapack_complex_double x[n_row][n_rhs];

    double rcond;
    double ferr[n_rhs];
    double berr[n_rhs];

    std::cout << "A : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_col; ++j) {
            a[i][j] = lapack_complex_double((double)std::rand() / RAND_MAX,
                                            (double)std::rand() / RAND_MAX);
            a0[i][j] = a[i][j];
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_rhs; ++j) {
            b[i][j] = lapack_complex_double((double)std::rand() / RAND_MAX,
                                            (double)std::rand() / RAND_MAX);
            b0[i][j] = b[i][j];
            std::cout << b[i][j] << " ";
        }
        std::cout << std::endl;
    }

    lapack_int info = LAPACKE_zgesvx(LAPACK_ROW_MAJOR, 'E', 'N',
                           n_row, n_rhs, *a, n_col, *af, n_row,
                           ipiv, &equed, r, c,
                           *b, n_rhs, *x, n_rhs,
                           &rcond, ferr, berr,
                           rpivot );

    std::cout << "A : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_col; ++j) {
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        for(int j = 0; j < n_rhs; ++j) {
            std::cout << b[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Ax - B : " << std::endl;
    for(int i = 0; i < n_row; ++i) {
        lapack_complex_double ax_b_i = 0.0;
        for(int j = 0; j < n_col; ++j) {
            ax_b_i += a0[i][j] * x[j][0];
        }
        std::cout << ax_b_i - b0[i][0] << std::endl;
    }

    std::cout << " condition : " << rcond << std::endl;
    std::cout << " ferr : " << ferr[0] << std::endl;
    std::cout << " berr : " << berr[0] << std::endl;

    return info;
}

void TestMatSolver() {
    srand (time(NULL));

    int N = 10;
    Matrix<std::complex<double>> A(N, N);
    Matrix<std::complex<double>> B(N, 1);

    std::cout << "A0 : " << std::endl;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            A[{i, j}] = std::complex<double>((double)std::rand() / RAND_MAX,
                                             (double)std::rand() / RAND_MAX);
            std::cout << A[{i, j}] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "B0 : " << std::endl;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < 1; ++j) {
            B[{i, j}] = std::complex<double>((double)std::rand() / RAND_MAX,
                                             (double)std::rand() / RAND_MAX);
            std::cout << B[{i, j}] << " ";
        }
        std::cout << std::endl;
    }

    Matrix<std::complex<double>> A0(A);
    Matrix<std::complex<double>> B0(B);

    Matrix<std::complex<double>> X = SolveLinear(A, B, false);


    std::cout << "A : " << std::endl;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            std::cout << A[{i, j}] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "B : " << std::endl;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < 1; ++j) {
            std::cout << B[{i, j}] << " ";
        }
        std::cout << std::endl;
    }


    std::cout << "X : " << std::endl;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < 1; ++j) {
            std::cout << X[{i, j}] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Ax - B : " << std::endl;
    for(int i = 0; i < N; ++i) {
        std::complex<double> ax_b_i = 0.0;
        for(int j = 0; j < N; ++j) {
            ax_b_i += A0[{i, j}] * X[{j, 0}];
        }
        std::cout << ax_b_i - B0(i, 0) << std::endl;
    }
}

void TestMatrixRead() {
    Matrix<std::complex<double>> A = ReadMatrixFromFile<std::complex<double>>("out/E_ksi.data");
    int n =  A.GetNCol();
    int n_row = (int)std::sqrt(n);
    assert(n_row * n_row == n);
    A.Reshape(n_row, n_row);

    Matrix<std::complex<double>> B = ReadMatrixFromFile<std::complex<double>>("out/A.data");
    n =  B.GetNCol();
    n_row = (int)std::sqrt(n);
    assert(n_row * n_row == n);
    B.Reshape(n_row, n_row);

    (A - B).Print();
}


#include "suitesparse/umfpack.h"

void TestUmfpack() {
    int n = 5;
    int Ap[] = {0, 2, 5, 9, 10, 12};
    int Ai[] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
    double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
    double b[] = {8., 45., -3., 3., 19.};
    double x[5];

    double *null = (double *) NULL;
    int i;
    void *Symbolic, *Numeric;
    (void) umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, null, null);
    (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic (&Symbolic);
    (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
    umfpack_di_free_numeric (&Numeric);
    for (i = 0 ; i < n ; i++) printf ("x [%d] = %g\n", i, x [i]) ;
}

void TestUmfpackComplex() {
    int n = 5;
    int Ap[] = {0, 2, 5, 9, 10, 12};
    int Ai[] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
    double Ax[] = {0.0, 2., 0.0, 3., 0.0, 3., 0.0, -1., 0.0, 4., 0.0, 4., 0.0, -3., 0.0, 1., 0.0, 2., 0.0, 2., 0.0, 6., 0.0, 1.};
    double b[] = {0.0, 8., 0.0, 45., 0.0, -3., 0.0, 3., 0.0, 19.};
    double x[10];

    double *null = (double *) NULL;
    int i;
    void *Symbolic, *Numeric;
    (void) umfpack_zi_symbolic(n, n, Ap, Ai, Ax, NULL, &Symbolic, null, null);
    (void) umfpack_zi_numeric(Ap, Ai, Ax, NULL, Symbolic, &Numeric, null, null);
    umfpack_zi_free_symbolic(&Symbolic);
    (void) umfpack_zi_solve (UMFPACK_A, Ap, Ai, Ax, NULL, x, NULL, b, NULL, Numeric, null, null);
    umfpack_zi_free_numeric (&Numeric);
    for (i = 0 ; i < n ; i++) printf ("x [%d] = %g, %g\n", i, x[2*i], x[2*i + 1]) ;
}

/*
#include <omp.h>

void TestOMP() {
    std::cout << "Limit: " << omp_get_max_threads() << std::endl;
    omp_set_num_threads(4);

    double a[100];
    double b[100];
    int i;
    #pragma omp parallel for
    for (i=1; i<100; i++) {
        std::cout << omp_get_thread_num() << "/" << omp_get_num_threads() << std::endl;
        b[i] = (a[i] + a[i-1]) / 2.0;
    }
}
*/

#endif


