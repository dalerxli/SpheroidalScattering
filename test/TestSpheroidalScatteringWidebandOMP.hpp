

#include <utility>
#include <string>


#include "boost/lexical_cast.hpp"

#include <sys/types.h>  // directory check
#include <sys/stat.h>

#include "SpheroidScattering.hpp"

#include <omp.h>


void TestSpheroidalScatteringWideband_OMP(int n_thread = 4) {
    std::cout << "thread limit: " << omp_get_max_threads() << std::endl;
    omp_set_num_threads(n_thread);

    double tipRadius= 50.0 * 1.0e-9;
    double length = 300.0 * 1.0e-6;
    double f0 = 0.01 * 1.0e12;
    double f1 = 3.0 * 1.0e12;
    int n_f = 100;

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

    #pragma omp parallel for
    for(int i = 0; i < n_f; ++i) {
        SpheroidScattering spheroid(tipRadius, length);

        double freq = f0 + (double)i / (n_f - 1) * (f1 - f0);
        int N_t = (int)(freq/speedOfLight*length*2.0 + 4);

        std::cout << omp_get_thread_num() << "/" << omp_get_num_threads();
        std::cout << "======================================================================" << std::endl;
        std::cout << " i: " << i << "  , f: " << freq/1.0e12 << "  , N_t: " << N_t << std::endl;

        spheroid.SetFrequency(freq);
        spheroid.SetIncidentAngle(M_PI/2.0);
        spheroid.SetFieldAmp(1.0);
        spheroid.SetNumberOfHarmonics(N_t);

        auto alpha_beta_gamma = spheroid.VerifySurfaceField();
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

