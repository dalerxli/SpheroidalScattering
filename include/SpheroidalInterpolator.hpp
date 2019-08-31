
#ifndef __SPHEROIDAL_INTERPOLATOR__
#define __SPHEROIDAL_INTERPOLATOR__


#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <complex>
#include <fstream>

#include <splinter/datatable.h>
#include <splinter/bspline.h>
#include <splinter/bsplinebuilder.h>

class SpheroidalInterpolator {
    public:
    SpheroidalInterpolator(int numVariables = 3) : bspline_Eeta_real(numVariables),
                                                   bspline_Eeta_imag(numVariables),
                                                   bspline_Eksi_real(numVariables),
                                                   bspline_Eksi_imag(numVariables),
                                                   bspline_Ephi_real(numVariables),
                                                   bspline_Ephi_imag(numVariables)
                                                   { assert(numVariables = 3); }
    ~SpheroidalInterpolator() { }

    void SetEtaLimits(double eta_0, double eta_1) {
        eta_limits = {eta_0, eta_1};
    }

    void SetEtaLimits(std::pair<double, double> limits) {
        eta_limits = limits;
    }

    void SetKsiLimits(double ksi_0, double ksi_1) {
        ksi_limits = {ksi_0, ksi_1};
    }

    void SetKsiLimits(std::pair<double, double> limits) {
        ksi_limits = limits;
    }

    void SetPhiLimits(double phi_0, double phi_1) {
        phi_limits = {phi_0, phi_1};
    }

    void SetPhiLimits(std::pair<double, double> limits) {
        phi_limits = limits;
    }

    void SetNumOfSamples(std::size_t n_eta, std::size_t n_ksi, std::size_t n_phi) {
        numOfSamples = {n_eta, n_ksi, n_phi};
    }

    void SetKsiNonuniformScale(int pow_scale) {
        assert(pow_scale >= 1);
        ksi_nonuniform_scale = pow_scale;
    }

    auto& GetMesh() {return mesh;}
    auto& GetE_eta() {return E_eta;}
    auto& GetE_ksi() {return E_ksi;}
    auto& GetE_phi() {return E_phi;}
    auto& GetNumOfSamples() {return numOfSamples;}
    auto& GetEtaLimits() {return eta_limits;}
    auto& GetKsiLimits() {return ksi_limits;}
    auto& GetPhiLimits() {return phi_limits;}
    auto GetKsiNonuniformScale() {return ksi_nonuniform_scale;}

    void SetupMesh(int _ksi_nonuniform_scale = -1) {
        if (_ksi_nonuniform_scale > 0) {
            ksi_nonuniform_scale = ksi_nonuniform_scale;
        }

        double& eta_0 = eta_limits.first;
        double& eta_1 = eta_limits.second;
        double& ksi_0 = ksi_limits.first;
        double& ksi_1 = ksi_limits.second;
        double& phi_0 = phi_limits.first;
        double& phi_1 = phi_limits.second;

        std::size_t& n_eta = numOfSamples[0];
        std::size_t& n_ksi = numOfSamples[1];
        std::size_t& n_phi = numOfSamples[2];

        double eta_i, ksi_j, phi_k;
        for(std::size_t i = 0; i < n_eta; ++i) {
            eta_i = eta_0 + (double)i/(n_eta - 1) * (eta_1 - eta_0);
            for(std::size_t j = 0; j < n_ksi; ++j) {
                if(ksi_nonuniform_scale == 1) {
                    ksi_j = ksi_0 + (double)j/(n_ksi - 1) * (ksi_1 - ksi_0);
                } else {
                    double d_ksi_scale = 1.0;
                    for(int j_2 = 0; j_2 < ksi_nonuniform_scale; ++j_2) {
                        d_ksi_scale *= (double)j/(n_ksi - 1);
                    }
                    ksi_j = ksi_0 +  d_ksi_scale * (ksi_1 - ksi_0);
                }
                for(std::size_t k = 0; k < n_phi; ++k) {
                    phi_k = phi_0 + (double)k/(n_phi - 1) * (phi_1 - phi_0);

                    mesh.emplace_back(std::array<double, 3>{eta_i, ksi_j, phi_k});
                }
            }
        }
    }

    void SetupSamples() {
        std::size_t n_pts = mesh.size();
        assert(E_eta.size() == n_pts && E_ksi.size() == n_pts && E_phi.size() == n_pts);

        SPLINTER::DenseVector x(3);
        for(std::size_t i = 0; i < n_pts; ++i) {
            auto& e_eta_i = E_eta[i];
            auto& e_ksi_i = E_ksi[i];
            auto& e_phi_i = E_phi[i];
            auto& r_i = mesh[i];


            x(0) = r_i[0];
            x(1) = r_i[1];
            x(2) = r_i[2];

            //std::cout << "r: " << x(0) << " " << x (1) << " " << x(2) << std::endl;
            //std::cout << "e: " << E_eta[i] << E_ksi[i] << E_phi[i] << std::endl;

            e_eta_real_Samples.addSample(x, e_eta_i.real());
            e_eta_imag_Samples.addSample(x, e_eta_i.imag());
            e_ksi_real_Samples.addSample(x, e_ksi_i.real());
            e_ksi_imag_Samples.addSample(x, e_ksi_i.imag());
            e_phi_real_Samples.addSample(x, e_phi_i.real());
            e_phi_imag_Samples.addSample(x, e_phi_i.imag());
        }

        //std::cout << "=========================" << std::endl;
    }

    void SetupBsplines(int degree = 3) {
        bspline_Eeta_real = SPLINTER::BSpline::Builder(e_eta_real_Samples).degree(degree).build();
        bspline_Eeta_imag = SPLINTER::BSpline::Builder(e_eta_imag_Samples).degree(degree).build();
        bspline_Eksi_real = SPLINTER::BSpline::Builder(e_ksi_real_Samples).degree(degree).build();
        bspline_Eksi_imag = SPLINTER::BSpline::Builder(e_ksi_imag_Samples).degree(degree).build();
        bspline_Ephi_real = SPLINTER::BSpline::Builder(e_phi_real_Samples).degree(degree).build();
        bspline_Ephi_imag = SPLINTER::BSpline::Builder(e_phi_imag_Samples).degree(degree).build();
    }

    void InterpolateFieldAtSpheroidalPoints(std::vector<std::array<double, 3>>& r_pts,
                                            std::vector<std::complex<double>>& e_eta,
                                            std::vector<std::complex<double>>& e_ksi,
                                            std::vector<std::complex<double>>& e_phi,
                                            std::vector<bool>& mask
                                            ) {
        // r_pts: [eta, ksi, phi]
        std::size_t n_pts = r_pts.size();
        SPLINTER::DenseVector x(3);

        e_eta.resize(n_pts);
        e_ksi.resize(n_pts);
        e_phi.resize(n_pts);

        //std::cout << "etalimits: " << eta_limits.first << ", " << eta_limits.second << std::endl;
        //std::cout << "ksilimits: " << ksi_limits.first << ", " << ksi_limits.second << std::endl;
        //std::cout << "philimits: " << phi_limits.first << ", " << phi_limits.second << std::endl;

        //double e_ksi_max = 0.0;

        for(std::size_t i = 0; i < n_pts; ++i) {
            if (mask[i]) {
                auto& r_i = r_pts[i];

                x(0) = r_i[0];
                x(1) = r_i[1];
                x(2) = r_i[2];

                //std::cout << "r: " << x(0) << " " << x (1) << " " << x(2) << std::endl;

                e_eta[i] = std::complex<double>(bspline_Eeta_real.eval(x), bspline_Eeta_imag.eval(x));
                e_ksi[i] = std::complex<double>(bspline_Eksi_real.eval(x), bspline_Eksi_imag.eval(x));
                e_phi[i] = std::complex<double>(bspline_Ephi_real.eval(x), bspline_Ephi_imag.eval(x));

                //std::cout << "e: " << e_eta[i] << e_ksi[i] << e_phi[i] << std::endl;

                //if (std::abs(e_ksi[i]) > e_ksi_max) {
                //    e_ksi_max = std::abs(e_ksi[i]);
                //}
            }
        }
        //std::cout << e_ksi_max << "  =================================" << std::endl;
    }

    void WriteMeshToFile(std::string fileName) {
        std::ofstream fileOut(fileName.c_str(), std::ios::out | std::ios::binary);
        assert(fileOut.is_open());

        fileOut.write((char*)(numOfSamples.data()), numOfSamples.size() * sizeof(std::size_t));
        fileOut.write((char*)(&eta_limits.first), sizeof(double));
        fileOut.write((char*)(&eta_limits.second), sizeof(double));
        fileOut.write((char*)(&ksi_limits.first), sizeof(double));
        fileOut.write((char*)(&ksi_limits.second), sizeof(double));
        fileOut.write((char*)(&phi_limits.first), sizeof(double));
        fileOut.write((char*)(&phi_limits.second), sizeof(double));

        fileOut.write((char*)(&ksi_nonuniform_scale), sizeof(int));

        fileOut.write((char*)(E_eta.data()), E_eta.size() * sizeof(std::complex<double>));
        fileOut.write((char*)(E_ksi.data()), E_ksi.size() * sizeof(std::complex<double>));
        fileOut.write((char*)(E_phi.data()), E_phi.size() * sizeof(std::complex<double>));
        fileOut.close();
    }


    void ReadMeshFromFile(std::string fileName) {
        std::ifstream fileIn(fileName.c_str(), std::ios::in | std::ios::binary);
        assert(fileIn.is_open());


        fileIn.read((char*)(numOfSamples.data()), numOfSamples.size() * sizeof(std::size_t));
        fileIn.read((char*)(&eta_limits.first), sizeof(double));
        fileIn.read((char*)(&eta_limits.second), sizeof(double));
        fileIn.read((char*)(&ksi_limits.first), sizeof(double));
        fileIn.read((char*)(&ksi_limits.second), sizeof(double));
        fileIn.read((char*)(&phi_limits.first), sizeof(double));
        fileIn.read((char*)(&phi_limits.second), sizeof(double));

        fileIn.read((char*)(&ksi_nonuniform_scale), sizeof(int));

        int n_sample_3d = numOfSamples[0]*numOfSamples[1]*numOfSamples[2];
        E_eta.resize(n_sample_3d);
        E_ksi.resize(n_sample_3d);
        E_phi.resize(n_sample_3d);
        fileIn.read((char*)(E_eta.data()), E_eta.size() * sizeof(std::complex<double>));
        fileIn.read((char*)(E_ksi.data()), E_ksi.size() * sizeof(std::complex<double>));
        fileIn.read((char*)(E_phi.data()), E_phi.size() * sizeof(std::complex<double>));
        fileIn.close();
    }

    private:
    std::pair<double, double> eta_limits;
    std::pair<double, double> ksi_limits;
    std::pair<double, double> phi_limits;

    int ksi_nonuniform_scale = 1;

    std::array<std::size_t, 3> numOfSamples;

    std::vector<std::array<double, 3>> mesh;
    std::vector<std::complex<double>> E_eta;
    std::vector<std::complex<double>> E_ksi;
    std::vector<std::complex<double>> E_phi;

    SPLINTER::DataTable e_eta_real_Samples;
    SPLINTER::DataTable e_eta_imag_Samples;
    SPLINTER::DataTable e_ksi_real_Samples;
    SPLINTER::DataTable e_ksi_imag_Samples;
    SPLINTER::DataTable e_phi_real_Samples;
    SPLINTER::DataTable e_phi_imag_Samples;

    SPLINTER::BSpline bspline_Eeta_real;
    SPLINTER::BSpline bspline_Eeta_imag;
    SPLINTER::BSpline bspline_Eksi_real;
    SPLINTER::BSpline bspline_Eksi_imag;
    SPLINTER::BSpline bspline_Ephi_real;
    SPLINTER::BSpline bspline_Ephi_imag;
};


#endif  // __SPHEROIDAL_INTERPOLATOR__



