#ifndef __SPHEROIDAL_FUNC__
#define __SPHEROIDAL_FUNC__

#include <complex>
#include <cassert>

#ifdef __cplusplus
extern "C" {
#endif
    // expansion coefficients d_mn for spheroidal functions
    void sdmn_(int* m, int* n, double* c, double* cv, int* kd, double* df);         // kd:  1  prolate
                                                                                    //      2  oblate
    
    // characteristic values for spheroidal functions
    void segv_(int* m, int* n, double* c, int* kd, double* cv, double* eg);         // kd:  1  prolate
                                                                                    //      2  oblate
    
    // prolate radial spheroidal functions of the first and second kind and derivatives
    void rswfp_(int* m, int* n, double* c, double* x, double* cv, int* kf,          // kf : 1  first kind and derivative
                double* r1f, double* r1d, double* r2f, double* r2d);                //      2  second kind and derivative
                                                                                    //      3  first and second kind and derivatives
    // oblate radial spheroidal functions of the first and second kind and derivatives
    void rswfo_(int* m, int* n, double* c, double* x, double* cv, int* kf,          // kf : 1  first kind and derivative
                double* r1f, double* r1d, double* r2f, double* r2d);                //      2  second kind and derivative

    // angular spheroidal function of the first kind and derivative
    void aswfa_(int* m, int* n, double* c, double* x, int* kd, double* cv, double* s1f, double* s1d); // kd :  1  prolate
                                                                                                      //      -1  oblate
                                                                                                      
    // angular spheroidal function of the first kind and derivative
    // based on associated legendre functions expansion
    void aswfb_(int* m, int* n, double* c, double* x, int* kd, double* cv, double* s1f, double* s1d); // kd :  1  prolate
                                                                                                      //      -1  oblate
                                                                                                      
#ifdef __cplusplus
}
#endif

void GetExpansionCoefficientsDmn(int m, int n, double c, double* dmn, int* dmn_inds, int* len_dmn, double thresh = 1.0e-15) {
    // dmn : preallocated array of length 200 ---> double[200]
    // d_inds: int[200]
    // the effective length will be output to len_mn 
    double cv;
    double eg[200];
    int kd = 1;
    double df[200];
    for(int i = 0; i < 200; ++i) {
        df[i] = 0.0;
    }
    
    // get characteristic value
    segv_(&m, &n, &c, &kd, &cv, eg);
        
    // get expansion coefficient
    sdmn_(&m, &n, &c, &cv, &kd, df);
    
    // get num of elems
    int n_elem = 0;
    for(int i = 199; i > 1; --i) {
        if(df[i] != 0.0) {
            n_elem = i + 1;
            break;
        }
    }
    
    // get max value
    double d_mn_max = 0.0;
    for(int i = 0; i < n_elem; ++i) {
        if(std::abs(df[i]) > d_mn_max) {
            d_mn_max = std::abs(df[i]);
        }
    }
    
    // filter based on threshold
    double d_mn_lim = d_mn_max * thresh;
    bool keep[200];
    int n_keep = 0;
    for(int i = 0; i < n_elem; ++i) {
        if(std::abs(df[i]) > d_mn_lim) {
            keep[i] = true;
            n_keep++;
        } else {
            keep[i] = false;
        }
    }
    
    int inds[200];
    if( (n - m) % 2 == 0 ) {
        for(int i = 0; i < n_elem; ++i) {
            inds[i] = 2 * i;
        }
    } else {
        for(int i = 0; i < n_elem; ++i) {
            inds[i] = 2 * i + 1;
        }
    }
    
    int i_kept = 0;
    for(int i = 0; i < n_elem; ++i) {
        if(keep[i]) {
            dmn[i_kept] = df[i];
            dmn_inds[i_kept] = inds[i];
            i_kept++; 
        }
    }
    
    if( m % 2 == 1 ) {
        for(int i = 0; i < n_keep; ++i) {
            dmn[i] *= -1.0;
        }
    }
    
    assert( i_kept == n_keep );
    *len_dmn = n_keep;
}

double GetProlateRadial1(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 1;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d); 
    
    return r1f;
}

double GetProlateRadial2(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 2;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d); 
    
    return r2f;
}

std::complex<double> GetProlateRadial3(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 3;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d); 
    
    return std::complex<double>(r1f, r2f);
}

std::complex<double> GetProlateRadial4(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 3;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d); 
    
    return std::complex<double>(r1f, -r2f);
}

double GetProlateRadialDerivative1(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 1;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d); 
    
    return r1d;
}

double GetProlateRadialDerivative2(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 2;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d); 
    
    return r2d;
}

std::complex<double> GetProlateRadialDerivative3(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 3;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d); 
    
    return std::complex<double>(r1d, r2d);
}

std::complex<double> GetProlateRadialDerivative4(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);


    // get prolate radials and derivatives
    int kf = 3;
    double r1f, r1d, r2f, r2d;
    rswfp_(&m, &n, &c, &x, &cv, &kf, &r1f, &r1d, &r2f, &r2d); 
    
    return std::complex<double>(r1d, -r2d);
}

double GetProlateAngular1(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);
    
    // get prolate radial and derivatives
    double s1f, s1d;
    aswfa_(&m, &n, &c, &x, &kd, &cv, &s1f, &s1d);

    return s1f;    
}

double GetProlateAngularDerivative1(int m, int n, double c, double x) {
    double cv;
    double eg[200];
    
    // get characteristic value
    int kd = 1; // prolate
    segv_(&m, &n, &c, &kd, &cv, eg);
    
    // get prolate radial and derivatives
    double s1f, s1d;
    aswfa_(&m, &n, &c, &x, &kd, &cv, &s1f, &s1d);

    return s1d;    
}


std::complex<double> GetM_mplus1n_o_plus_eta(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return  ((ksi*ksi - 1.0)*GetProlateRadialDerivative4(m, n, c, ksi) - ksi*m*GetProlateRadial4(m, n, c, ksi)) \
                            *GetProlateAngular1(m, n, c, eta)*std::cos(phi*(m + 1)) \
                            /(d*std::sqrt(-eta*eta + ksi*ksi)*std::sqrt(ksi*ksi - 1.0));
}

std::complex<double> GetM_mplus1n_o_plus_ksi(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return ((eta*eta - 1)*GetProlateAngularDerivative1(m, n, c, eta) - eta*m*GetProlateAngular1(m, n, c, eta))
                         *GetProlateRadial4(m, n, c, ksi)*std::cos(phi*(m + 1))
                         /(d*std::sqrt(-eta*eta + 1)*std::sqrt(-eta*eta + ksi*ksi));
}

std::complex<double> GetM_mplus1n_o_plus_phi(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return -(eta*(ksi*ksi - 1)*GetProlateAngular1(m, n, c, eta)*GetProlateRadialDerivative4(m, n, c, ksi)
             -ksi*(eta*eta - 1)*GetProlateRadial4(m, n, c, ksi)*GetProlateAngularDerivative1(m, n, c, eta))
            *std::sin(phi*m + phi)/(d*(eta*eta - ksi*ksi));
}

std::complex<double> GetM_mminus1n_o_minus_eta(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return  ((1.0 - ksi*ksi)*GetProlateRadialDerivative4(m, n, c, ksi) - ksi*m*GetProlateRadial4(m, n, c, ksi))
            *GetProlateAngular1(m, n, c, eta)
            *std::cos(phi*(m - 1))
            /(d*std::sqrt(-eta*eta + ksi*ksi)*std::sqrt(ksi*ksi - 1));
}

std::complex<double> GetM_mminus1n_o_minus_ksi(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return  ((1.0-eta*eta)*GetProlateAngularDerivative1(m, n, c, eta) - eta*m*GetProlateAngular1(m, n, c, eta))
            *GetProlateRadial4(m, n, c, ksi)
            *std::cos(phi*(m - 1))/(d*std::sqrt(-eta*eta + 1)*std::sqrt(-eta*eta + ksi*ksi));
}

std::complex<double> GetM_mminus1n_o_minus_phi(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return  -(eta*(ksi*ksi - 1)*GetProlateAngular1(m, n, c, eta)*GetProlateRadialDerivative4(m, n, c, ksi) 
              - ksi*(eta*eta - 1)*GetProlateRadial4(m, n, c, ksi)*GetProlateAngularDerivative1(m, n, c, eta))
              *std::sin(phi*(m - 1))/(d*(eta*eta - ksi*ksi));
}

std::complex<double> GetM_mn_o_z_eta(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return  -2.0*eta*m*GetProlateRadial4(m, n, c, ksi)*GetProlateAngular1(m, n, c, eta)
                *std::cos(phi*m)/(d*std::sqrt(-eta*eta + 1)*std::sqrt(-eta*eta + ksi*ksi));
}

std::complex<double> GetM_mn_o_z_ksi(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return  2.0*ksi*m*GetProlateRadial4(m, n, c, ksi)*GetProlateAngular1(m, n, c, eta)
               *std::cos(phi*m)/(d*std::sqrt(-eta*eta + ksi*ksi)*std::sqrt(ksi*ksi - 1));
}

std::complex<double> GetM_mn_o_z_phi(double eta, double ksi, double phi, int m, int n, double c, double d) {
    return  2.0*std::sqrt(-eta*eta + 1)*std::sqrt(ksi*ksi - 1)*(-eta*GetProlateRadial4(m, n, c, ksi)
               *GetProlateAngularDerivative1(m, n, c, eta) 
             + ksi*GetProlateAngular1(m, n, c, eta)*GetProlateRadialDerivative4(m, n, c, ksi))
                  *std::sin(phi*m)/(d*(eta*eta - ksi*ksi));
}



#endif

