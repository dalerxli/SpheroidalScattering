#ifndef __SPHEROID_SCATTERING__
#define __SPHEROID_SCATTERING__

#include <cstddef>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include <map>
#include <utility>

#include "SpheroidalFunc.hpp"
#include "SpheroidalIntegrals.hpp"
#include "MatAlgebra.hpp"
#include "SpheroidalInterpolator.hpp"


constexpr double speedOfLight = 2.99792458e8;    // m/s

class SpheroidScattering {
    public:
    SpheroidScattering(double tip_radius, double ellipsod_length) {
        tipRadius = tip_radius;
        length = ellipsod_length;

        GetProlateSpheroidParameters(tipRadius, length,
                                      ellipse_a, ellipse_b, spheroid_ksi, spheroid_d);

        numOfHarmonics = 8;
    }

    void GetProlateSpheroidParameters(double tipRadius, double length,
                                      double& a, double& b, double& ksi, double& d) {
        double b2_div_a = tipRadius;
        a = length/2.0;
        b = std::sqrt(b2_div_a * a);

        // d*ksi = a     d*sqrt(ksi**2 - 1) = b
        // ksi**2 * (1 - (b/a)**2) = 1
        ksi = 1.0/(1.0 - (b/a)*(b/a));
        d = a / ksi;

        return;
    }

    void SetFrequency(double f) {
        frequency = f;

        wavelength = speedOfLight / frequency;
        wavenumber = 2.0 * M_PI / wavelength;
        spheroid_c = wavenumber * spheroid_d / 2.0;
    }

    void SetIncidentAngle(double angle_rad) {
        incidenceAngle = angle_rad;
    }

    void SetFieldAmp(std::complex<double> e0) {
        e0_incidence = e0;
    }

    void SetNumberOfHarmonics(int n_max) {
        numOfHarmonics = n_max;
    }

    void Map2DIndexTo1D(const int m_0, const int m_1, const int ind_start,
                        std::map<std::pair<int, int>, int>& map2DTo1D,
                        std::map<int, std::pair<int, int>>& map1DTo2D
                        ) {
        // m = m_0 .. m_1-1   n = m ... m_1-1
        int ind = ind_start;
        for(int m = m_0; m < m_1; ++m) {
            for(int n = m; n < m_1; ++n) {
                map2DTo1D[{m,n}] = ind;
                map1DTo2D[ind] = {m, n};
                ind += 1;
            }
        }
    }

    std::complex<double> GetIncExpansionCoeffs_Amn(const int m, const int n) {
        const std::complex<double>& E0 = e0_incidence;
        const double theta_0 = incidenceAngle;
        const double k = wavenumber;
        const double c = spheroid_c;

        double eps_m = 2.0;
        if( m == 0) {
            eps_m = 1.0;
        }
        double N_mn = GetInt_Sm_mpn_Sm_mpN(c, m, n-m, n-m);
        double A_mn = 2.0 * eps_m * GetProlateAngular1(m, n, c, std::cos(theta_0)) / N_mn;

        std::complex<double> j_nm1 = 1.0;
        if( n == 0 ) {
            j_nm1 = std::complex<double>(0.0, -1.0);
        } else {
            int nm1_4 = (n - 1) % 4;
            for(int i = 0; i < nm1_4; ++i) {
                j_nm1 *= std::complex<double>(0.0, 1.0);
            }
        }

        return E0 / k * j_nm1 * A_mn;
    }

    std::array<Matrix<std::complex<double>>, 2> ConstructMatrix() {
        const double k = wavenumber;
        const double theta_0 = incidenceAngle;
        const double ksi_0 = spheroid_ksi;
        const double c_0 = spheroid_c;

        const std::complex<double>& E0 = e0_incidence;
        const int N_t = numOfHarmonics;

        std::map<std::pair<int, int>, int> alphaInd_2DTo1D;
        std::map<int, std::pair<int, int>> alphaInd_1DTo2D;
        Map2DIndexTo1D(0, N_t, 0, alphaInd_2DTo1D, alphaInd_1DTo2D);
        int n_total = alphaInd_2DTo1D.size();

        std::map<std::pair<int, int>, int> betaInd_2DTo1D;
        std::map<int, std::pair<int, int>> betaInd_1DTo2D;
        Map2DIndexTo1D(1, N_t+1, n_total, betaInd_2DTo1D, betaInd_1DTo2D);
        n_total += betaInd_2DTo1D.size();

        assert(N_t > 2);
        std::vector<int> gammaInd(N_t - 1);
        for(std::size_t i = 0; i < N_t - 1; ++i) {
            gammaInd[i] = n_total + i;
        }
        n_total += gammaInd.size();

        // construct coefficient marix
        Matrix<std::complex<double>> A(n_total, n_total);
        Matrix<std::complex<double>> b(n_total, 1);

        // eta: cos(m*phi) m=1..Nt
        for(int m = 0; m < N_t; ++m) {
            for(int N = m; N < N_t; ++N) {
                int ind_row = alphaInd_2DTo1D[{m, N}];
                for(int n = m; n < N_t; ++n) {
                    int ind_col = alphaInd_2DTo1D[{m, n}];
                    std::complex<double> elem =
                            ((ksi_0*ksi_0 - 1)*GetProlateRadialDerivative4(m, n, c_0, ksi_0) \
                            - ksi_0*m*GetProlateRadial4(m, n, c_0, ksi_0)) \
                            * GetInt_Sm_mpn_Sm_mpN(c_0, m, n-m, N-m);

                    A(ind_row, ind_col) += elem;

                    ind_col = betaInd_2DTo1D[{m+1, n+1}];
                    elem = -2.0*std::sqrt(ksi_0*ksi_0 - 1)*(m+1)*GetProlateRadial4(m+1, n+1, c_0, ksi_0) \
                            * GetInt_Smp1_mpnp1_Sm_mpN_x_div_sqrt_1mx2(c_0, m, n-m, N-m);

                    A(ind_row, ind_col) += elem;

                    // ---- rhs
                    std::complex<double> A_mn = GetIncExpansionCoeffs_Amn(m, n);
                    b(ind_row, 0) -= A_mn * \
                                  ( \
                                      -ksi_0*m*GetProlateRadial1(m, n, c_0, ksi_0) \
                                      + (ksi_0*ksi_0 - 1)*GetProlateRadialDerivative1(m, n, c_0, ksi_0) \
                                  ) * GetInt_Sm_mpn_Sm_mpN(c_0, m, n-m, N-m);

                    std::complex<double> A_mp2np2 = GetIncExpansionCoeffs_Amn(m+2, n+2);
                    b(ind_row, 0) -= A_mp2np2 * \
                                  ( \
                                      ksi_0*(m+2)*GetProlateRadial1(m+2, n+2, c_0, ksi_0) \
                                      + (ksi_0*ksi_0 - 1)*GetProlateRadialDerivative1(m+2, n+2, c_0, ksi_0) \
                                  ) * GetInt_Smp2_mpnp2_Sm_mpN(c_0, m, n-m, N-m);
                    if( m==0 ) {
                        std::complex<double> A_0n = GetIncExpansionCoeffs_Amn(0, n);
                        b(ind_row, 0) -= A_0n * (ksi_0*ksi_0 - 1)*GetProlateRadialDerivative1(m, n, c_0, ksi_0) \
                                           * GetInt_Sm_mpn_Sm_mpN(c_0, m, n-m, N-m);
                    }
                }
            }
        }

        // eta: cos(0*phi)
        for(int N = 0; N < N_t - 1; ++N) {
            int ind_row = gammaInd[N];
            for(int n = 0; n < N_t - 1; ++n) {
                int ind_col = gammaInd[n];

                std::complex<double> elem =
                        (-(ksi_0*ksi_0 - 1)*GetProlateRadialDerivative4(1, n+1, c_0, ksi_0) \
                        - ksi_0*1*GetProlateRadial4(1, n+1, c_0, ksi_0)) \
                        * GetInt_Sm_mpn_Sm_mpN(c_0, 1, n, N);

                A(ind_row, ind_col) += elem;

                // rhs
                std::complex<double> A_1np1 = GetIncExpansionCoeffs_Amn(1, n+1);
                b(ind_row, 0) -= A_1np1 * \
                                  ( \
                                      ksi_0*GetProlateRadial1(1, n+1, c_0, ksi_0) \
                                      + (ksi_0*ksi_0 - 1)*GetProlateRadialDerivative1(1, n+1, c_0, ksi_0) \
                                  ) * GetInt_Sm_mpn_Sm_mpN(c_0, 1, n, N);
            }
        }

        // phi: sin(m*phi), m=1...Nt-2
        for(int m = 0; m < N_t; ++m) {
            for(int N = m; N < N_t; ++N) {
                int ind_row = betaInd_2DTo1D[{m+1, N+1}];
                for(int n = m; n < N_t; ++n) {

                    int ind_col = alphaInd_2DTo1D[{m, n}];
                    std::complex<double> elem =
                        (ksi_0*ksi_0 - 1.0)*GetProlateRadialDerivative4(m, n, c_0, ksi_0) \
                            * GetInt_Sm_mpn_Sm_mpN_x(c_0, m, n-m, N-m) \
                        +  \
                            ksi_0 * GetProlateRadial4(m, n, c_0, ksi_0) \
                            * GetInt_dxSm_mpn_Sm_mpN_1mx2(c_0, m, n-m, N-m);

                    A(ind_row, ind_col) += elem;

                    ind_col = betaInd_2DTo1D[{m+1, n+1}];
                    elem = 2.0*std::sqrt(ksi_0*ksi_0 - 1.0) * \
                            ( \
                                GetProlateRadial4(m+1, n+1, c_0, ksi_0) \
                                * GetInt_dxSmp1_mpnp1_Sm_mpN_x_sqrt_1mx2(c_0, m, n-m, N-m) \
                            - \
                                ksi_0 * GetProlateRadialDerivative4(m+1, n+1, c_0, ksi_0) \
                                * GetInt_Smp1_mpnp1_Sm_mpN_sqrt_1mx2(c_0, m, n-m, N-m) \
                            );

                    A(ind_row, ind_col) += elem;

                    //---- rhs
                    std::complex<double> A_mn = GetIncExpansionCoeffs_Amn(m, n);
                    b(ind_row, 0) -= A_mn * \
                        ( \
                           ksi_0*GetProlateRadial1(m, n, c_0, ksi_0)*GetInt_dxSm_mpn_Sm_mpN_1mx2(c_0, m, n-m, N-m) \
                           + (ksi_0*ksi_0 - 1)*GetProlateRadialDerivative1(m, n, c_0, ksi_0) \
                                           *GetInt_Sm_mpn_Sm_mpN_x(c_0, m, n-m, N-m) \
                        );

                    std::complex<double> A_mp2np2 = GetIncExpansionCoeffs_Amn(m+2, n+2);
                    b(ind_row, 0) += A_mp2np2 * \
                        ( \
                           ksi_0*GetProlateRadial1(m+2, n+2, c_0, ksi_0)*GetInt_dxSmp2_mpnp2_Sm_mpN_1mx2(c_0, m, n-m, N-m) \
                           + (ksi_0*ksi_0 - 1)*GetProlateRadialDerivative1(m+2, n+2, c_0, ksi_0) \
                                           *GetInt_Smp2_mpnp2_Sm_mpN_x(c_0, m, n-m, N-m) \
                        );
                    if( m==0 ) {
                        std::complex<double> A_0n = GetIncExpansionCoeffs_Amn(0, n);
                        b(ind_row, 0) -= A_0n * \
                        ( \
                           ksi_0*GetProlateRadial1(m, n, c_0, ksi_0)*GetInt_dxSm_mpn_Sm_mpN_1mx2(c_0, m, n-m, N-m) \
                           + (ksi_0*ksi_0 - 1)*GetProlateRadialDerivative1(m, n, c_0, ksi_0) \
                                           *GetInt_Sm_mpn_Sm_mpN_x(c_0, m, n-m, N-m) \
                        );
                    }
                }
            }
        }

        return std::array<Matrix<std::complex<double>>, 2>{A, b};
    }

    std::array<Matrix<std::complex<double>>, 3> GetAlphaBetaGamma_from_X(Matrix<std::complex<double>>& x) {
        const int N_t = numOfHarmonics;

        std::map<std::pair<int, int>, int> alphaInd_2DTo1D;
        std::map<int, std::pair<int, int>> alphaInd_1DTo2D;
        Map2DIndexTo1D(0, N_t, 0, alphaInd_2DTo1D, alphaInd_1DTo2D);
        int n_total = alphaInd_2DTo1D.size();
        int n_end_alpha = n_total;

        std::map<std::pair<int, int>, int> betaInd_2DTo1D;
        std::map<int, std::pair<int, int>> betaInd_1DTo2D;
        Map2DIndexTo1D(1, N_t+1, n_total, betaInd_2DTo1D, betaInd_1DTo2D);
        n_total += betaInd_2DTo1D.size();
        int n_end_beta = n_total;

        assert(N_t > 2);
        std::vector<int> gammaInd(N_t - 1);
        for(std::size_t i = 0; i < N_t - 1; ++i) {
            gammaInd[i] = n_total + i;
        }
        n_total += gammaInd.size();

        Matrix<std::complex<double>> alpha(N_t, N_t);
        Matrix<std::complex<double>> beta(N_t+1, N_t+1);
        Matrix<std::complex<double>> gamma(N_t, 1);

        for(std::size_t i = 0; i < n_end_alpha; ++i) {
            alpha[alphaInd_1DTo2D[i]] = x(i, 0);
        }
        for(std::size_t i = n_end_alpha; i < n_end_beta; ++i) {
            beta[betaInd_1DTo2D[i]] = x(i, 0);
        }
        for(std::size_t i = n_end_beta; i < x.GetNRow(); ++i) {
            gamma(i - n_end_beta + 1, 0) = x(i, 0);
        }

        return std::array<Matrix<std::complex<double>>, 3>{alpha, beta, gamma};
    }

    void GetETMonSurface_direct(std::vector<double>& etas, double ksi_0, double phi_0,
            std::vector<std::complex<double>>& E_eta,
            std::vector<std::complex<double>>& E_ksi
            ) {
        const std::complex<double>& E_0 = e0_incidence;
        const double k = wavenumber;
        const double d = spheroid_d;

        assert(phi_0 == 0);

        std::size_t n_eta = etas.size();
        E_eta.resize(n_eta);
        E_ksi.resize(n_eta);

        const std::complex<double> _1j(0.0, 1.0);
        for(std::size_t i = 0; i < n_eta; ++i) {
            double eta = etas[i];
            double z_hat_eta = ksi_0*std::sqrt((1.0 - eta*eta)/(ksi_0*ksi_0 - eta*eta));
            double z_hat_ksi = eta * std::sqrt((ksi_0*ksi_0 - 1)/(ksi_0*ksi_0 - eta*eta));
            double x = d/2*std::sqrt(1.0 - eta*eta)*std::sqrt(ksi_0*ksi_0 - 1)*std::cos(phi_0);
            E_eta[i] = E_0*std::exp(_1j*k*x)*z_hat_eta;
            E_ksi[i] = E_0*std::exp(_1j*k*x)*z_hat_ksi;
        }
    }

    void GetETMonSurface_expansion(std::vector<double>& etas, double ksi_0, double phi_0,
            std::vector<std::complex<double>>& E_eta,
            std::vector<std::complex<double>>& E_ksi
            ) {
        const std::complex<double>& E_0 = e0_incidence;
        const double k = wavenumber;
        const double d = spheroid_d;
        const double c_0 = spheroid_c;

        assert( phi_0 == 0 );

        const double theta_0 = incidenceAngle; //np.pi/2;
        assert( theta_0 == M_PI/2);

        std::size_t n_eta = etas.size();
        E_eta.resize(n_eta);
        E_ksi.resize(n_eta);
        const int N = numOfHarmonics;
        for(int i = 0; i < n_eta; ++i) {
            double eta = etas[i];
            for(int m = 0; m < N; ++m) {
                for(int n = m; n < N; ++n) {
                    std::complex<double> A_mn = GetIncExpansionCoeffs_Amn(m, n);
                    E_eta[i] += A_mn * 2.0*(ksi_0*ksi_0 - 1)*GetProlateRadialDerivative1(m, n, c_0, ksi_0) \
                                         *GetProlateAngular1(m, n, c_0, eta) \
                                        /(d*std::sqrt(ksi_0*ksi_0 - eta*eta)*std::sqrt(ksi_0*ksi_0 - 1));
                    E_ksi[i] += A_mn * (-2.0)*(1.0 - eta*eta)*GetProlateAngularDerivative1(m, n, c_0, eta) \
                                        *GetProlateRadial1(m, n, c_0, ksi_0) \
                                        /(d*std::sqrt(ksi_0*ksi_0 - eta*eta)*std::sqrt(1.0 - eta*eta));
                }
            }
        }
    }

    void GetFieldOnSurface(Matrix<std::complex<double>>& alpha,
                           Matrix<std::complex<double>>& beta,
                           Matrix<std::complex<double>>& gamma,
                           std::vector<double>& etas, double ksi, double phi,
                           std::vector<std::complex<double>>& E_eta,
                           std::vector<std::complex<double>>& E_ksi,
                           std::vector<std::complex<double>>& E_phi
                           ) {
        const double d = spheroid_d;
        const double c = spheroid_c;

        int n_pts = etas.size();
        E_eta.resize(n_pts);
        E_ksi.resize(n_pts);
        E_phi.resize(n_pts);
        for(int i = 0; i < n_pts; ++i) {
            double eta = etas[i];
            int M = alpha.GetNRow();
            int N = alpha.GetNCol();
            for(int m = 0; m < M; ++m) {
                for(int n = m; n < N; ++n) {
                    E_eta[i] += alpha(m, n)*GetM_mplus1n_o_plus_eta(eta, ksi, phi, m, n, c, d);
                    E_ksi[i] += alpha(m, n)*GetM_mplus1n_o_plus_ksi(eta, ksi, phi, m, n, c, d);
                    E_phi[i] += alpha(m, n)*GetM_mplus1n_o_plus_phi(eta, ksi, phi, m, n, c, d);
                }
            }

            M = beta.GetNRow();
            N = beta.GetNCol();
            for(int m = 0; m < M; ++m) {
                for(int n = m; n < N; ++n) {
                    E_eta[i] += beta(m, n)*GetM_mn_o_z_eta(eta, ksi, phi, m, n, c, d);
                    E_ksi[i] += beta(m, n)*GetM_mn_o_z_ksi(eta, ksi, phi, m, n, c, d);
                    E_phi[i] += beta(m, n)*GetM_mn_o_z_phi(eta, ksi, phi, m, n, c, d);
                }
            }

            assert( gamma.GetNCol() == 1 );
            N = gamma.GetNRow();
            for(int n = 1; n < N; ++n) {
                E_eta[i] += gamma(n, 0)*GetM_mminus1n_o_minus_eta(eta, ksi, phi, 1, n, c, d);
                E_ksi[i] += gamma(n, 0)*GetM_mminus1n_o_minus_ksi(eta, ksi, phi, 1, n, c, d);
                E_phi[i] += gamma(n, 0)*GetM_mminus1n_o_minus_phi(eta, ksi, phi, 1, n, c, d);
            }
        }
    }


    void GetFieldAtSpheroidalPoints(
                           Matrix<std::complex<double>>& alpha,
                           Matrix<std::complex<double>>& beta,
                           Matrix<std::complex<double>>& gamma,
                           std::vector<std::array<double, 3>>& r_pts,   //(eta, ksi, phi)
                           std::vector<std::complex<double>>& E_eta,
                           std::vector<std::complex<double>>& E_ksi,
                           std::vector<std::complex<double>>& E_phi
                           ) {
        const double d = spheroid_d;
        const double c = spheroid_c;

        int n_pts = r_pts.size();
        E_eta.resize(n_pts);
        E_ksi.resize(n_pts);
        E_phi.resize(n_pts);
        for(int i = 0; i < n_pts; ++i) {
            double eta = r_pts[i][0];
            double ksi = r_pts[i][1];
            double phi = r_pts[i][2];
            int M = alpha.GetNRow();
            int N = alpha.GetNCol();
            for(int m = 0; m < M; ++m) {
                for(int n = m; n < N; ++n) {
                    E_eta[i] += alpha(m, n)*GetM_mplus1n_o_plus_eta(eta, ksi, phi, m, n, c, d);
                    E_ksi[i] += alpha(m, n)*GetM_mplus1n_o_plus_ksi(eta, ksi, phi, m, n, c, d);
                    E_phi[i] += alpha(m, n)*GetM_mplus1n_o_plus_phi(eta, ksi, phi, m, n, c, d);
                }
            }

            M = beta.GetNRow();
            N = beta.GetNCol();
            for(int m = 0; m < M; ++m) {
                for(int n = m; n < N; ++n) {
                    E_eta[i] += beta(m, n)*GetM_mn_o_z_eta(eta, ksi, phi, m, n, c, d);
                    E_ksi[i] += beta(m, n)*GetM_mn_o_z_ksi(eta, ksi, phi, m, n, c, d);
                    E_phi[i] += beta(m, n)*GetM_mn_o_z_phi(eta, ksi, phi, m, n, c, d);
                }
            }

            assert( gamma.GetNCol() == 1 );
            N = gamma.GetNRow();
            for(int n = 1; n < N; ++n) {
                E_eta[i] += gamma(n, 0)*GetM_mminus1n_o_minus_eta(eta, ksi, phi, 1, n, c, d);
                E_ksi[i] += gamma(n, 0)*GetM_mminus1n_o_minus_ksi(eta, ksi, phi, 1, n, c, d);
                E_phi[i] += gamma(n, 0)*GetM_mminus1n_o_minus_phi(eta, ksi, phi, 1, n, c, d);
            }
        }
    }
    
    SpheroidalInterpolator GetFieldInterpolator (
        Matrix<std::complex<double>>& alpha,
        Matrix<std::complex<double>>& beta,
        Matrix<std::complex<double>>& gamma,
        std::array<std::pair<double, double>, 3>& coord_limits,   //((eta0, eta1), (ksi0, ksi1), (phi0, phi1))
        std::array<std::size_t, 3>& numOfSamples
    ) {
        SpheroidalInterpolator sphInterpolator;
        sphInterpolator.SetEtaLimits(coord_limits[0]);
        sphInterpolator.SetKsiLimits(coord_limits[1]);
        sphInterpolator.SetPhiLimits(coord_limits[2]);
        
        sphInterpolator.SetNumOfSamples(numOfSamples[0], numOfSamples[1], numOfSamples[2]);
        
        sphInterpolator.SetupMesh();
        
        GetFieldAtSpheroidalPoints(alpha, beta, gamma, sphInterpolator.GetMesh(), 
                                                       sphInterpolator.GetE_eta(),
                                                       sphInterpolator.GetE_ksi(),
                                                       sphInterpolator.GetE_phi()
                                                       );
        int degree = 3;
        sphInterpolator.SetupSamples();
        sphInterpolator.SetupBsplines(degree);
        
        return sphInterpolator;        
    }
    
    void GetFieldAtCartesianPoints(Matrix<std::complex<double>>& alpha,
                                   Matrix<std::complex<double>>& beta,
                                   Matrix<std::complex<double>>& gamma,
                                   std::vector<std::array<double, 3>>& r_pts,
                                   std::vector<std::complex<double>>& E_eta,
                                   std::vector<std::complex<double>>& E_ksi,
                                   std::vector<std::complex<double>>& E_phi,
                                   bool totalField = true) {
        // r = [x, y, z]
        const auto& a = ellipse_a;
        const auto& c = spheroid_c;
        const auto& d = spheroid_d;
        const auto& k = wavenumber;
        const auto& E0 = e0_incidence;

        std::size_t n_pts = r_pts.size();
        E_eta.resize(n_pts);
        E_ksi.resize(n_pts);
        E_phi.resize(n_pts);
        double eta, ksi, phi;
        int M, N;

        std::complex<double> _1j(0.0, 1.0);

        for(std::size_t i = 0; i < n_pts; ++i) {
            auto& r = r_pts[i];
            double x = r[0];
            double y = r[1];
            double z = r[2];
            CoordinatePointTransformRectToSpheroid(x, y, z, eta, ksi, phi);
            M = alpha.GetNRow();
            N = alpha.GetNCol();
            for(int m = 0; m < M; ++m) {
                for(int n = m; n < N; ++n) {
                    E_eta[i] += alpha(m, n)*GetM_mplus1n_o_plus_eta(eta, ksi, phi, m, n, c, d);
                    E_ksi[i] += alpha(m, n)*GetM_mplus1n_o_plus_ksi(eta, ksi, phi, m, n, c, d);
                    E_phi[i] += alpha(m, n)*GetM_mplus1n_o_plus_phi(eta, ksi, phi, m, n, c, d);
                }
            }

            M = beta.GetNRow();
            N = beta.GetNCol();
            for(int m = 0; m < M; ++m) {
                for(int n = m; n < N; ++n) {
                    E_eta[i] += beta(m, n)*GetM_mn_o_z_eta(eta, ksi, phi, m, n, c, d);
                    E_ksi[i] += beta(m, n)*GetM_mn_o_z_ksi(eta, ksi, phi, m, n, c, d);
                    E_phi[i] += beta(m, n)*GetM_mn_o_z_phi(eta, ksi, phi, m, n, c, d);
                }
            }

            N = gamma.GetNRow();
            assert(gamma.GetNCol() == 1);
            for(int n = 1; n < N; ++n) {
                E_eta[i] += gamma(n, 0)*GetM_mminus1n_o_minus_eta(eta, ksi, phi, 1, n, c, d);
                E_ksi[i] += gamma(n, 0)*GetM_mminus1n_o_minus_ksi(eta, ksi, phi, 1, n, c, d);
                E_phi[i] += gamma(n, 0)*GetM_mminus1n_o_minus_phi(eta, ksi, phi, 1, n, c, d);
            }

            if(totalField) {
                auto z_hat_eta = ksi * std::sqrt((1.0 - eta*eta)/(ksi*ksi - eta*eta));
                auto z_hat_ksi = eta * std::sqrt((ksi*ksi - 1.0)/(ksi*ksi - eta*eta));
                //x = d/2*np.sqrt(1 - eta**2)*np.sqrt(ksi**2 - 1)*np.cos(phi);
                assert(incidenceAngle == M_PI/2.0);
                E_eta[i] += E0*std::exp(_1j*k*x)*z_hat_eta;
                E_ksi[i] += E0*std::exp(_1j*k*x)*z_hat_ksi;
            }
        }
    }

    void CoordinatePointTransformSpheroidToRect(double eta, double ksi, double phi,
                                                double& x, double& y, double& z) {
        const auto& d = spheroid_d;
        x = d/2*std::sqrt((1.0 - eta*eta))*std::sqrt((ksi*ksi - 1.0))*std::cos(phi);
        y = d/2*std::sqrt((1.0 - eta*eta))*std::sqrt((ksi*ksi - 1.0))*std::sin(phi);
        z = d/2*eta*ksi;
    }

    void CoordinatePointTransformRectToSpheroid(double x, double y, double z,
                                                double& eta, double& ksi, double& phi) {
        const auto& d = spheroid_d;
        ksi = (std::sqrt(x*x + y*y + (z + d/2)*(z + d/2)) + std::sqrt(x*x + y*y + (z - d/2)*(z - d/2)))/d;
        eta = (std::sqrt(x*x + y*y + (z + d/2)*(z + d/2)) - std::sqrt(x*x + y*y + (z - d/2)*(z - d/2)))/d;
        phi = std::atan2(y, x);
    }

    auto GetFieldAroundTipAtXZPlane(double Dx, double Dz, int nx, int nz,
                                    Matrix<std::complex<double>>& alpha,
                                    Matrix<std::complex<double>>& beta,
                                    Matrix<std::complex<double>>& gamma,
                                    bool totalField = true) {
        const auto& ksi_0 = spheroid_ksi;
        const auto& a = ellipse_a;
        const auto& c = spheroid_c;
        const auto& d = spheroid_d;
        const auto& E0 = e0_incidence;
        const auto& k = wavenumber;

        std::vector<double> x(nx);
        for(int i = 0; i < nx; ++i) {
            x[i] = -Dx/2.0 + (double)i * Dx / (nx - 1);
        }

        std::vector<double> z(nz);
        for(int i = 0; i < nz; ++i) {
            z[i] = a/2.0 - Dz/2.0 + (double)i * Dz / (nz - 1);
        }

        std::vector<std::array<double, 3>> r_pts;
        std::vector<std::array<int, 2>> r_inds;
        double eta, ksi, phi;
        for(int i = 0; i < nx; ++i) {
            double x_i = x[i];
            for(int j = 0; j < nz; ++j) {
                double z_j = z[j];
                CoordinatePointTransformRectToSpheroid(x_i, 0.0, z_j, eta, ksi, phi);
                if( ksi > ksi_0 ) {
                    r_pts.emplace_back(std::array<double, 3>{x_i, 0.0, z_j});
                    r_inds.emplace_back(std::array<int, 2>{i, j});
                }
            }
        }

        std::vector<std::complex<double>> E_eta;
        std::vector<std::complex<double>> E_ksi;
        std::vector<std::complex<double>> E_phi;
        GetFieldAtCartesianPoints(alpha, beta, gamma, r_pts, E_eta, E_ksi, E_phi, true);

        Matrix<std::complex<double>> E_ksi_mesh(nx, nz);
        for(int i = 0; i < r_inds.size(); ++i) {
            E_ksi_mesh(r_inds[i][0], r_inds[i][1]) = E_ksi[i];
        }

        return E_ksi_mesh;
    }

    std::array<Matrix<std::complex<double>>, 3> VerifySurfaceField() {
        std::cout << "Number of harmonics : " << numOfHarmonics << std::endl;

        auto A_b = ConstructMatrix();
        auto& A = A_b[0];
        auto& b = A_b[1];

        auto x = SolveLinear(A, b, true);
        //std::cout << "x : " << std::endl;
        //x.Print();

        std::cout << "error : " << (A*x - b).GetNorm() << std::endl;

        auto alpha_beta_gamma = GetAlphaBetaGamma_from_X(x);
        auto& alpha = alpha_beta_gamma[0];
        auto& beta = alpha_beta_gamma[1];
        auto& gamma = alpha_beta_gamma[2];

        double eps = 1.0e-5;
        int n_pts = 20;
        double eta_0 = 0.0;
        double eta_1 = 1.0 - eps;
        std::vector<double> etas(n_pts);
        for(int i = 0; i < n_pts; ++i) {
            etas[i] = eta_0 + (double)i * (eta_1 - eta_0) / (n_pts - 1);
        }

        double phi_0 = 0.0;
        double ksi_0 = spheroid_ksi;
        double d = spheroid_d;

        std::vector<std::complex<double>> E_eta0, E_ksi0;
        GetETMonSurface_direct(etas, ksi_0, phi_0, E_eta0, E_ksi0);

        std::vector<std::complex<double>> E_eta1, E_ksi1;
        GetETMonSurface_expansion(etas, ksi_0, phi_0, E_eta1, E_ksi1);

        std::vector<std::complex<double>> E_eta2, E_ksi2, E_phi2;
        GetFieldOnSurface(alpha, beta, gamma, etas, ksi_0, phi_0, E_eta2, E_ksi2, E_phi2);

        for(int i = 0; i < n_pts; ++i) {
            std::cout <<  E_eta0[i] << "  "
                      <<  E_eta1[i] << "  "
                      <<  E_eta2[i] << "  "
                      <<  std::abs(E_eta2[i]/E_eta0[i])
                      << std::endl;
        }

        std::cout << "------------------------" << std::endl;

        for(int i = 0; i < n_pts; ++i) {
            printf("(%.4e, %.4e)  (%.4e, %.4e)  (%.4e, %.4e)   %.5f \n", E_ksi0[i].real(), E_ksi0[i].imag(),
                                                                         E_ksi1[i].real(), E_ksi1[i].imag(),
                                                                         E_ksi2[i].real(), E_ksi2[i].imag(),
                                                                         std::abs(E_ksi2[i]));
        }

        return alpha_beta_gamma;
    }

    private:
    double tipRadius;
    double length;
    double ellipse_a;
    double ellipse_b;
    double spheroid_ksi;
    double spheroid_d;
    double spheroid_c;

    double frequency;
    double wavelength;
    double wavenumber;

    std::complex<double> e0_incidence;
    double incidenceAngle;

    int numOfHarmonics;
};

#endif

