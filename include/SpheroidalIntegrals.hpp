#ifndef __SPHEROIDAL_INTEGRALS__
#define __SPHEROIDAL_INTEGRALS__

#include <complex>
#include <cassert>
#include <map>
#include <utility>

#include "SpheroidalFunc.hpp"

double Get_Int_PmnPmn(int m, int n) {
    if(m > n) {
        return 0.0;
    }
    double fact_mpn_div_fact_mmn = 1.0;
    for(int i = n + m; i > std::max(n - m, 0); --i) {
        fact_mpn_div_fact_mmn *= i;
    }
    return 2.0/(2*n + 1)*fact_mpn_div_fact_mmn;
}

double Get_Int_Pm2muPmnu(int m, int mu, int nu) {
    if( m + 2 > mu || m > nu) {
        return 0.0;
    }
    if( mu < nu ) {
        return 0.0;
    } else if( mu == nu ) {
        double fact_div = 1.0;
        for( int i = nu + m; i > std::max(nu - m - 2, 0); --i) {
            fact_div *= i;    
        }
        return -2.0/(2.0*nu + 1)*fact_div;
    } else if( (mu + nu) % 2 == 0) {
        double fact_div = 1.0;
        for(int i = nu + m; i > std::max(nu - m, 0); --i) {
            fact_div *= i;
        }               
        return 4.0*(m + 1)*fact_div;
    } else {
        return 0.0;
    }
}

double Get_Int_PmnPmn_div_1mx2(int m, int n) {
    if( m > n ) {
        return 0.0;
    }
    if( m != 0) {
        double fact_div = 1.0;
        for( int i = n + m; i > std::max(n - m, 0); --i) {
            fact_div *= i;
        }
        return fact_div / m;
    } else {
        std::cout << "Invalid value. m = " << m << std::endl; 
        assert(false);
    }
}

std::map<std::pair<int, int>, double> 
    Expand_P_minusm_n(int m, int n) {
    double coeff = 0.0;
    if( std::abs(m) <= n ) {
        double fact_mpn_div_fact_mmn = 1.0;
        for(int i = n + m; i > std::max(n - m, 0); --i) {
            fact_mpn_div_fact_mmn *= i;
        }
        coeff = 1.0 / fact_mpn_div_fact_mmn;
        if(m % 2 == 1) {
            coeff *= -1.0;
        }
    }
    std::pair<int, int> pair(m, n);
    std::map<std::pair<int, int>, double> map;
    map[pair] = coeff;
    return map;
}

std::map<std::pair<int, int>, double> 
    Expand_xPmn(int m, int n) {
    return std::map<std::pair<int, int>, double>{{std::pair<int, int>(m, n-1), (double)(n + m)/(2*n + 1)},
                                                 {std::pair<int, int>(m, n+1), (double)(n - m + 1)/(2*n + 1)}
                                                };
}

std::map<std::pair<int, int>, double> 
    Expand_x_div_sqrt_1mx2_Pmn(int m, int n) {
    return std::map<std::pair<int, int>, double>{{std::pair<int, int>(m-1, n), -(double)(n + m)*(n - m + 1)/(2*m)}, 
                                                 {std::pair<int, int>(m+1, n), -1.0/(2*m)}
                                                };
}
    
std::map<std::pair<int, int>, double> 
    Expand_1_div_sqrt_1mx2_Pmn__nm1(int m, int n) {
    return std::map<std::pair<int, int>, double>{{std::pair<int, int>(m-1, n-1), -(double)(n + m)*(n + m - 1)/(2*m)}, 
                                                 {std::pair<int, int>(m+1, n-1), -1.0/(2*m)}
                                                };
}

std::map<std::pair<int, int>, double> 
    Expand_1_div_sqrt_1mx2_Pmn__np1(int m, int n) {
    return std::map<std::pair<int, int>, double>{{std::pair<int, int>(m-1, n+1), -(double)(n - m + 1)*(n - m + 2)/(2*m)}, 
                                                 {std::pair<int, int>(m+1, n+1), -1.0/(2*m)}
                                                };
}

std::map<std::pair<int, int>, double> 
    Expand_1mx2_dxPmn(int m, int n) {
    return std::map<std::pair<int, int>, double>{{std::pair<int, int>(m, n-1), (double)(n + 1)*(n + m)/(2*n + 1)}, 
                                                 {std::pair<int, int>(m, n+1), -(double)n*(n - m + 1)/(2*n + 1)}
                                                };
}

std::map<std::pair<int, int>, double> 
    Expand_sqrt_1mx2_Pmn__mm1(int m, int n) {
    return std::map<std::pair<int, int>, double>{{std::pair<int, int>(m-1, n-1), -(double)(n + m)*(n + m - 1)/(2*n + 1)}, 
                                                 {std::pair<int, int>(m-1, n+1), (double)(n - m + 1)*(n - m + 2)/(2*n + 1)}
                                                };
}

std::map<std::pair<int, int>, double> 
    Expand_sqrt_1mx2_Pmn__mp1(int m, int n) {
    return std::map<std::pair<int, int>, double>{{std::pair<int, int>(m+1, n-1), 1.0/(2*n + 1)}, 
                                                 {std::pair<int, int>(m+1, n+1), -1.0/(2*n + 1)}
                                                };
}

std::map<std::pair<int, int>, double> 
    Expand_sqrt_1mx2_dxPmn(int m, int n) {
    return std::map<std::pair<int, int>, double>{{std::pair<int, int>(m-1, n), (double)(n + m)*(n - m + 1)/2}, 
                                                 {std::pair<int, int>(m+1, n), -1.0/2}
                                                };
}


// I_2mNn
double GetInt_Sm_mpn_Sm_mpN(double c, int m, int n, int N) {
    if( n == N ) {
        double d_mn[200];
        int d_inds[200];
        int len_d = 0;
        GetExpansionCoefficientsDmn(m, m+n, c, d_mn, d_inds, &len_d);
        double I = 0.0;
        for(int i = 0; i < len_d; ++i) {
            I += d_mn[i]*d_mn[i] * Get_Int_PmnPmn(m, m+d_inds[i]);
        }
        return I;
    } else {
        return 0.0;
    }
}

// I_2mNn
double GetInt_Smp1_mpnp1_Sm_mpN_x_div_sqrt_1mx2(double c, int m, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(m+1, m+n+1, c, d_q, d_q_inds, &len_q);

    int m_q = m + 1;

    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(m, m+N, c, d_r, d_r_inds, &len_r);

    int m_r = m;
        
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        auto q_expansion = Expand_x_div_sqrt_1mx2_Pmn(m_q, m_q + d_q_inds[q]);
        for( auto& q_exp_elem : q_expansion ) {
            auto& mn_pair = q_exp_elem.first;
            int m_q_expan = mn_pair.first;
            int n_q_expan = mn_pair.second;
            double coeff = q_exp_elem.second; //q_expansion[mn_pair];
            if(m_q_expan == m_r) {
                for(int r = 0; r < len_r; ++r) {
                    if( n_q_expan == m_r + d_r_inds[r] ) {
                        I += coeff*d_q[q]*d_r[r]*Get_Int_PmnPmn(m_q_expan, n_q_expan);
                    }
                }
            } else if( m_q_expan == m_r + 2 ) {
                for(int r = 0; r < len_r; ++r) {
                    I += coeff*d_q[q]*d_r[r]*Get_Int_Pm2muPmnu(m_r, n_q_expan, m_r + d_r_inds[r]);
                }
            } else {
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}

// I_3mNn
double GetInt_Sm_mpn_Sm_mpN_x(double c, int m, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(m, m+n, c, d_q, d_q_inds, &len_q);

    int m_q = m;

    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(m, m+N, c, d_r, d_r_inds, &len_r);

    int m_r = m;
    
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        auto q_expansion = Expand_xPmn(m_q, m_q + d_q_inds[q]);
        for(auto& q_exp_elem : q_expansion ) {
            auto& mn_pair = q_exp_elem.first;
            int m_q_expan = mn_pair.first;
            int n_q_expan = mn_pair.second;
            double coeff = q_exp_elem.second; //q_expansion[mn_pair]
            if(m_q_expan == m_r) {
                for(int r = 0; r < len_r; ++r) {
                    if(n_q_expan == m_r + d_r_inds[r]) {
                        I += coeff*d_q[q]*d_r[r]*Get_Int_PmnPmn(m_q_expan, n_q_expan);
                    }
                }
            } else {
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}

// I_4mNn
double GetInt_dxSm_mpn_Sm_mpN_1mx2(double c, int m, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(m, m+n, c, d_q, d_q_inds, &len_q);

    int m_q = m;

    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(m, m+N, c, d_r, d_r_inds, &len_r);

    int m_r = m;
        
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        auto q_expansion = Expand_1mx2_dxPmn(m_q, m_q + d_q_inds[q]);
        for(auto& q_exp_elem : q_expansion ) {
            auto& mn_pair = q_exp_elem.first;
            int m_q_expan = mn_pair.first;
            int n_q_expan = mn_pair.second;
            double coeff = q_exp_elem.second; //q_expansion[mn_pair]
            if(m_q_expan == m_r) {
                for(int r = 0; r < len_r; ++r) {
                    if(n_q_expan == m_r + d_r_inds[r]) {
                        I += coeff*d_q[q]*d_r[r]*Get_Int_PmnPmn(m_q_expan, n_q_expan);
                    }
                }
            } else {
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}

// I_5mNn
double GetInt_Smp1_mpnp1_Sm_mpN_sqrt_1mx2(double c, int m, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(m+1, m+n+1, c, d_q, d_q_inds, &len_q);

    int m_q = m + 1;

    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(m, m+N, c, d_r, d_r_inds, &len_r);

    int m_r = m;
    
    double I = 0.0;
    for(int q = 0;  q < len_q; ++q) {
        auto q_expansion = Expand_sqrt_1mx2_Pmn__mm1(m_q, m_q + d_q_inds[q]);
        for(auto& q_exp_elem : q_expansion ) {
            auto& mn_pair = q_exp_elem.first;
            int m_q_expan = mn_pair.first;
            int n_q_expan = mn_pair.second;
            double coeff = q_exp_elem.second; //q_expansion[mn_pair]
            if(m_q_expan == m_r) {
                for(int r = 0; r < len_r; ++r) {
                    if(n_q_expan == m_r + d_r_inds[r]) {
                        I += coeff*d_q[q]*d_r[r]*Get_Int_PmnPmn(m_q_expan, n_q_expan);
                    }
                }
            } else {
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}

// I_6mNn
double GetInt_dxSmp1_mpnp1_Sm_mpN_x_sqrt_1mx2(double c, int m, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(m+1, m+n+1, c, d_q, d_q_inds, &len_q);

    int m_q = m + 1;

    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(m, m+N, c, d_r, d_r_inds, &len_r);

    int m_r = m;
    
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        auto q_expansion = Expand_sqrt_1mx2_dxPmn(m_q, m_q + d_q_inds[q]);
        for(auto& q_exp_elem : q_expansion ) {
            auto& mn_q_pair = q_exp_elem.first;
            int m_q_expan = mn_q_pair.first; 
            int n_q_expan = mn_q_pair.second;
            double coeff_q = q_exp_elem.second; // q_expansion[mn_q_pair]
            if( m_q_expan == m_r ) {
                for(int r = 0; r < len_r; ++r) {
                    auto r_expansion = Expand_xPmn(m_r, m_r + d_r_inds[r]);
                    for(auto& r_exp_elem : r_expansion ) {
                        auto& mn_r_pair = r_exp_elem.first;
                        int m_r_expan = mn_r_pair.first;
                        int n_r_expan = mn_r_pair.second;
                        double coeff_r = r_exp_elem.second; //r_expansion[mn_r_pair]
                        assert(m_r_expan == m_q_expan);
                        if(n_q_expan == n_r_expan) {
                            I += coeff_q*coeff_r*d_q[q]*d_r[r]*Get_Int_PmnPmn(m_q_expan, n_q_expan);
                        }
                    }
                }
            } else if(m_q_expan == m_r + 2) {
                for(int r = 0; r < len_r; ++r) {
                    auto r_expansion = Expand_xPmn(m_r, m_r + d_r_inds[r]);
                    for(auto& r_exp_elem : r_expansion) {
                        auto& mn_r_pair = r_exp_elem.first;
                        int m_r_expan = mn_r_pair.first;
                        int n_r_expan = mn_r_pair.second;
                        double coeff_r = r_exp_elem.second; //r_expansion[mn_r_pair]
                        assert( m_r_expan + 2 == m_q_expan );
                        I += coeff_q*coeff_r*d_q[q]*d_r[r]*Get_Int_Pm2muPmnu(m_r_expan, n_q_expan, n_r_expan);
                    }
                }
            } else {
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}

// I_7mNn
double GetInt_Smp2_mpnp2_Sm_mpN(double c, int m, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(m+2, m+n+2, c, d_q, d_q_inds, &len_q);

    double m_q = m + 2;
    
    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(m, m+N, c, d_r, d_r_inds, &len_r);

    int m_r = m;
    
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        for(int r = 0; r < len_r; ++r) {
            I += d_q[q]*d_r[r]*Get_Int_Pm2muPmnu(m_r, m_q + d_q_inds[q], m_r + d_r_inds[r]);
        }
    }
    return I;
}

// I_8mNn
double GetInt_Smp2_mpnp2_Sm_mpN_x(double c, int m, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(m+2, m+n+2, c, d_q, d_q_inds, &len_q);

    int m_q = m + 2;

    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(m, m+N, c, d_r, d_r_inds, &len_r);

    int m_r = m;
        
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        auto q_expansion = Expand_xPmn(m_q, m_q + d_q_inds[q]);
        for(auto& q_exp_elem : q_expansion) {
            auto& mn_pair = q_exp_elem.first;
            int m_q_expan = mn_pair.first;
            int n_q_expan = mn_pair.second;
            double coeff = q_exp_elem.second; //q_expansion[mn_pair]
            if( m_q_expan == m_r + 2 ) {
                for(int r = 0; r < len_r; ++ r) {
                    I += coeff*d_q[q]*d_r[r]*Get_Int_Pm2muPmnu(m_r, n_q_expan, m_r + d_r_inds[r]);
                }
            } else {
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}


// I_9mNn
double GetInt_dxSmp2_mpnp2_Sm_mpN_1mx2(double c, int m, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(m+2, m+n+2, c, d_q, d_q_inds, &len_q);

    int m_q = m + 2;
    
    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(m, m+N, c, d_r, d_r_inds, &len_r);

    int m_r = m;
    
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        auto q_expansion = Expand_1mx2_dxPmn(m_q, m_q + d_q_inds[q]);
        for(auto& q_exp_elem : q_expansion) {
            auto& mn_pair = q_exp_elem.first;
            int m_q_expan = mn_pair.first;
            int n_q_expan = mn_pair.second;
            double coeff = q_exp_elem.second; //q_expansion[mn_pair]
            if( m_q_expan == m_r + 2 ) {
                for(int r = 0;  r < len_r; ++r) {
                    I += coeff*d_q[q]*d_r[r]*Get_Int_Pm2muPmnu(m_r, n_q_expan, m_r + d_r_inds[r]);
                }
            } else {
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}


// I_10Nn
double GetInt_S0_n_S1_1pN_sqrt_1mx2(double c, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(0, n, c, d_q, d_q_inds, &len_q);

    int m_q = 0;
    
    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(1, 1+N, c, d_r, d_r_inds, &len_r);

    int m_r = 1;
    
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        auto q_expansion = Expand_sqrt_1mx2_Pmn__mp1(m_q, m_q + d_q_inds[q]);
        for(auto& q_exp_elem : q_expansion) {
            auto& mn_pair = q_exp_elem.first;
            int m_q_expan = mn_pair.first;
            int n_q_expan = mn_pair.second;
            double coeff = q_exp_elem.second; //q_expansion[mn_pair]
            if( m_q_expan == m_r ) {
                for(int r = 0; r < len_r; ++r) {
                    if( n_q_expan == m_r + d_r_inds[r] ) {
                        I += coeff*d_q[q]*d_r[r]*Get_Int_PmnPmn(m_q_expan, n_q_expan);
                    }
                }
            } else {
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}


// I_11Nn
double GetInt_dxS0_n_S1_1pN_sqrt_x_sqrt_1mx2(double c, int n, int N) {
    double d_q[200];
    int d_q_inds[200];
    int len_q = 0;
    GetExpansionCoefficientsDmn(0, n, c, d_q, d_q_inds, &len_q);

    int m_q = 0;

    double d_r[200];
    int d_r_inds[200];
    int len_r = 0;
    GetExpansionCoefficientsDmn(1, 1+N, c, d_r, d_r_inds, &len_r);

    int m_r = 1;
    
    double I = 0.0;
    for(int q = 0; q < len_q; ++q) {
        auto q_expansion = Expand_sqrt_1mx2_dxPmn(m_q, m_q + d_q_inds[q]);
        for(auto& q_exp_elem : q_expansion) {
            auto& mn_q_pair = q_exp_elem.first;
            int m_q_expan = mn_q_pair.first;
            int n_q_expan = mn_q_pair.second;
            double coeff_q = q_exp_elem.second; //q_expansion[mn_q_pair]
            if( m_q_expan < 0 ) {
                auto minusm_expansion = Expand_P_minusm_n(abs(m_q_expan), n_q_expan);
                assert( minusm_expansion.size() == 1);
                int count = 0;
                for(auto& minusm_exp_elem : minusm_expansion) {
                    auto& mn_minusm_pair = minusm_exp_elem.first;
                    m_q_expan = mn_minusm_pair.first;
                    n_q_expan = mn_minusm_pair.second;
                    coeff_q *= minusm_exp_elem.second;
                    count++;
                }
                assert(count == 1);
            }
            if( m_q_expan == m_r ) {
                for(int r = 0; r < len_r; ++r) {
                    auto r_expansion = Expand_xPmn(m_r, m_r + d_r_inds[r]);
                    for(auto& r_exp_elem : r_expansion) {
                        auto& mn_r_pair = r_exp_elem.first;
                        int m_r_expan = mn_r_pair.first;
                        int n_r_expan = mn_r_pair.second;
                        double coeff_r = r_exp_elem.second; //r_expansion[mn_r_pair]
                        assert( m_r_expan == m_q_expan );
                        if( n_q_expan == n_r_expan ) {
                            I += coeff_q*coeff_r*d_q[q]*d_r[r]*Get_Int_PmnPmn(m_q_expan, n_q_expan);
                        }
                    }
                }
            } else{
                std::cout << "error: invalid input .. " << std::endl;
                assert(false);
            }
        }
    }
    return I;
}



#endif

