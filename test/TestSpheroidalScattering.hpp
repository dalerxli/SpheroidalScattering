
#include <utility>
#include <string>
#include <vector>
#include <string>

#include "boost/lexical_cast.hpp"

#include <sys/types.h>  // directory check
#include <sys/stat.h>

#include "SpheroidScattering.hpp"

void TestSpheroidVectorTransforms() {

    double eta = 0.9;
    double ksi = 1.3;
    double phi = M_PI/6.0;

    std::complex<double> e_x(2.0, 1.0);
    std::complex<double> e_y(1.0, 2.0);
    std::complex<double> e_z(3.0, 4.0);

    std::complex<double> e_eta;
    std::complex<double> e_ksi;
    std::complex<double> e_phi;

    std::cout << "e_x: " << e_x << std::endl;
    std::cout << "e_y: " << e_y << std::endl;
    std::cout << "e_z: " << e_z << std::endl;

    SpheroidScattering::VectorTransformFromRectToSpheroid(eta, ksi, phi, e_x, e_y, e_z, e_eta, e_ksi, e_phi);

    std::cout << "e_eta: " << e_eta << std::endl;
    std::cout << "e_ksi: " << e_ksi << std::endl;
    std::cout << "e_phi: " << e_phi << std::endl;

    e_x = 0.0;
    e_y = 0.0;
    e_z = 0.0;

    SpheroidScattering::VectorTransformFromSpheroidToRect(eta, ksi, phi, e_eta, e_ksi, e_phi, e_x, e_y, e_z);

    std::cout << "e_x: " << e_x << std::endl;
    std::cout << "e_y: " << e_y << std::endl;
    std::cout << "e_z: " << e_z << std::endl;

}

void TestSpheroidalScattering() {
    double tipRadius= 50.0 * 1.0e-9;
    double length = 900.0 * 1.0e-6;
    double freq = 1.0 * 1.0e12;

    int N_t = (int)(freq/speedOfLight*length*3.0 + 4);

    SpheroidScattering spheroid(tipRadius, length);
    spheroid.SetFrequency(freq);
    spheroid.SetIncidentAngle(M_PI/2.0);
    spheroid.SetFieldAmp(1.0);
    spheroid.SetNumberOfHarmonics(N_t);

    std::complex<double> tip_field;
    auto alpha_beta_gamma = spheroid.VerifySurfaceField(&tip_field);
    auto& alpha = alpha_beta_gamma[0];
    auto& beta = alpha_beta_gamma[1];
    auto& gamma = alpha_beta_gamma[2];

    std::cout << "tip_field: " << tip_field << std::endl;

    bool plotFields = true;
    if(plotFields) {
        double Dx = 4.0*tipRadius;
        double Dz = 4.0*tipRadius;
        int nx = 100;
        int nz = 100;
        auto E_ksi_mesh = spheroid.GetFieldAroundTipAtXZPlane(Dx, Dz, nx, nz, alpha, beta, gamma);

        E_ksi_mesh.WriteToFile("out/E_ksi.data");
    }
}

void TestSpheroidalScatteringWideband() {
    double tipRadius= 50.0 * 1.0e-9;
    double length = 300.0 * 1.0e-6;
    double f0 = 0.01 * 1.0e12;
    double f1 = 3.0 * 1.0e12;
    int n_f = 200;

    std::string folder("out/");
    folder += std::string("R=") + boost::lexical_cast<std::string>(tipRadius/1.0e-9)
                        + "_L=" + boost::lexical_cast<std::string>(length/1.0e-6)
                        + "_f0=" + boost::lexical_cast<std::string>(f0/1.0e12)
                        + "_f1=" + boost::lexical_cast<std::string>(f1/1.0e12)
                        + "_nf=" + boost::lexical_cast<std::string>(n_f)
                        ;


    int dir_err = mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        printf("Error creating directory %s!\n", folder.c_str());
    }

    folder += "/";

    SpheroidScattering spheroid(tipRadius, length);

    for(int i = 0; i < n_f; ++i) {
        double freq = f0 + (double)i / (n_f - 1) * (f1 - f0);
        int N_t = (int)(freq/speedOfLight*length*3.0 + 4);

        std::cout << "======================================================================" << std::endl;
        std::cout << " i: " << i << "  , f: " << freq/1.0e12 << "  , N_t: " << N_t << std::endl;

        spheroid.SetFrequency(freq);
        spheroid.SetIncidentAngle(M_PI/2.0);
        spheroid.SetFieldAmp(1.0);
        spheroid.SetNumberOfHarmonics(N_t);

        std::complex<double> tip_field;
        auto alpha_beta_gamma = spheroid.VerifySurfaceField(&tip_field);
        auto& alpha = alpha_beta_gamma[0];
        auto& beta = alpha_beta_gamma[1];
        auto& gamma = alpha_beta_gamma[2];


        std::string fileSuffix = std::string("_")                                               \
                                        + "_i=" + boost::lexical_cast<std::string>(i)           \
                                        + "_f=" + boost::lexical_cast<std::string>(freq)        \
                                        + "_nt=" + boost::lexical_cast<std::string>(N_t);

        alpha.WriteToFile(folder + std::string("A") + fileSuffix);
        beta.WriteToFile(folder + std::string("B") + fileSuffix);
        gamma.WriteToFile(folder + std::string("G") + fileSuffix);

        Matrix<std::complex<double>> tipMat(1, 1);
        tipMat(0, 0) = tip_field;
        tipMat.WriteToFile(folder + std::string("T") + fileSuffix);

        bool plotFields = true;
        if(plotFields) {
            double Dx = 4.0*tipRadius;
            double Dz = 4.0*tipRadius;
            int nx = 40;
            int nz = 40;
            auto E_ksi_mesh = spheroid.GetFieldAroundTipAtXZPlane(Dx, Dz, nx, nz, alpha, beta, gamma);

            std::string eksi_fileSuffix = std::string("_")
                                          + "_nx=" + boost::lexical_cast<std::string>(nx)
                                          + "_nz=" + boost::lexical_cast<std::string>(nz)
                                          + "_Dx=" + boost::lexical_cast<std::string>(Dx)
                                          + "_Dz=" + boost::lexical_cast<std::string>(Dz);


            E_ksi_mesh.WriteToFile(folder + "E_ksi" + fileSuffix + eksi_fileSuffix + ".data");
        }

    }
}


void GetThreadStartAndEndFrequencyIndices(int n_f, int n_thread, int ind_thread, std::vector<int>& inds_f) {
    inds_f.clear();
    int i = 0;
    while(i < n_f) {
        for(int j = 0; j < n_thread; ++j) {
            if(j == ind_thread) {
                inds_f.push_back(i);
            }
            ++i;

            if(i >= n_f) {
                break;
            }
        }
    }
}


void SpheroidalScatteringWideband_ThreadJob(int ind_thread, int n_thread) {
    assert(n_thread >= 1 && ind_thread < n_thread);
    double tipRadius= 50.0 * 1.0e-9;
    double length = 600.0 * 1.0e-6;
    double f0 = 0.01 * 1.0e12;
    double f1 = 3.0 * 1.0e12;
    int n_f = 200;

    std::string folder("out/");
    folder += std::string("R=") + boost::lexical_cast<std::string>(tipRadius/1.0e-9)
                        + "_L=" + boost::lexical_cast<std::string>(length/1.0e-6)
                        + "_f0=" + boost::lexical_cast<std::string>(f0/1.0e12)
                        + "_f1=" + boost::lexical_cast<std::string>(f1/1.0e12)
                        + "_nf=" + boost::lexical_cast<std::string>(n_f)
                        ;


    if( ind_thread == 0 ) {
        int dir_err = mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err)
        {
            printf("Error creating directory %s!\n", folder.c_str());
        }
    }

    folder += "/";

    std::vector<int> inds_f;
    GetThreadStartAndEndFrequencyIndices(n_f, n_thread, ind_thread, inds_f);

    for(int i : inds_f) {
        SpheroidScattering spheroid(tipRadius, length);

        double freq = f0 + (double)i / (n_f - 1) * (f1 - f0);
        int N_t = (int)(freq/speedOfLight*length*3.0 + 4);

        std::cout << "======================================================================" << std::endl;
        std::cout << " i: " << i << "  , f: " << freq/1.0e12 << "  , N_t: " << N_t << std::endl;

        spheroid.SetFrequency(freq);
        spheroid.SetIncidentAngle(M_PI/2.0);
        spheroid.SetFieldAmp(1.0);
        spheroid.SetNumberOfHarmonics(N_t);

        std::complex<double> tip_field;
        auto alpha_beta_gamma = spheroid.VerifySurfaceField(&tip_field);
        auto& alpha = alpha_beta_gamma[0];
        auto& beta = alpha_beta_gamma[1];
        auto& gamma = alpha_beta_gamma[2];


        std::string fileSuffix = std::string("_")                                               \
                                        + "_i=" + boost::lexical_cast<std::string>(i)           \
                                        + "_f=" + boost::lexical_cast<std::string>(freq)        \
                                        + "_nt=" + boost::lexical_cast<std::string>(N_t);

        alpha.WriteToFile(folder + std::string("A") + fileSuffix);
        beta.WriteToFile(folder + std::string("B") + fileSuffix);
        gamma.WriteToFile(folder + std::string("G") + fileSuffix);

        Matrix<std::complex<double>> tipMat(1, 1);
        tipMat(0, 0) = tip_field;
        tipMat.WriteToFile(folder + std::string("T") + fileSuffix);

        bool plotFields = true;
        if(plotFields) {
            double Dx = 4.0*tipRadius;
            double Dz = 4.0*tipRadius;
            int nx = 100;
            int nz = 100;
            auto E_ksi_mesh = spheroid.GetFieldAroundTipAtXZPlane(Dx, Dz, nx, nz, alpha, beta, gamma);

            std::string eksi_fileSuffix = std::string("_")
                                          + "_nx=" + boost::lexical_cast<std::string>(nx)
                                          + "_nz=" + boost::lexical_cast<std::string>(nz)
                                          + "_Dx=" + boost::lexical_cast<std::string>(Dx)
                                          + "_Dz=" + boost::lexical_cast<std::string>(Dz);


            E_ksi_mesh.WriteToFile(folder + "E_ksi" + fileSuffix + eksi_fileSuffix + ".data");
        }
    }
}

#include <thread>

void TestSpheroidalScatteringWideband_Threaded(int n_thread = 4) {
    std::vector<std::thread> thread_pool;
    for(int i = 0; i < n_thread; ++i) {
        thread_pool.emplace_back(SpheroidalScatteringWideband_ThreadJob, i, n_thread);
    }

    for(int i = 0; i < n_thread; ++i) {
        thread_pool[i].join();
    }
}

