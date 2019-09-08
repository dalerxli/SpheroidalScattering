
#include <utility>
#include <string>
#include <vector>
#include <string>

#include "boost/lexical_cast.hpp"

#include <sys/types.h>  // directory check
#include <sys/stat.h>
#include <filesystem>

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

    /*
    std::filesystem::path folder_path = folder;
    if ( std::filesystem::exists(folder_path) ) {
        std::filesystem::remove_all(folder_path);
        //std::filesystem::remove(folder.c_str());
    }
    std::filesystem::create_directories(folder_path);*/

    struct stat info;

    if( stat( folder.c_str(), &info ) != 0 ) {
        printf( "cannot access %s\n", folder.c_str() );
    } else if( info.st_mode & S_IFDIR ) {
        printf( "%s is a directory\n", folder.c_str() );
        //delete
        std::string syscommand = std::string("rm ") + folder + "/*";
        std::system( syscommand.c_str() );
        std::cout << "content delete!" << std::endl;
    }

    int dir_err = mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        printf("Error creating directory %s!\n", folder.c_str());
    }

    folder += "/";

    SpheroidScattering spheroid(tipRadius, length);

    for(int i = 0; i < n_f; ++i) {
        if (i == 3) break;
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


void SpheroidalScatteringWideband_ThreadJob(int ind_thread, int n_thread,
                                            double tipRadius, double length, double f0, double f1, int n_f,
                                            std::string folder) {
    assert(n_thread >= 1 && ind_thread < n_thread);
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

        std::string interpolator_fileSuffix = std::string("_")                                               \
                                        + "_i=" + boost::lexical_cast<std::string>(i)           \
                                        + "_f=" + boost::lexical_cast<std::string>(freq);
        auto& sphInterpolator = spheroid.GetSpheroidalInterpolator();
        sphInterpolator.WriteMeshToFile(folder + "interpMesh" + interpolator_fileSuffix + ".data");

    }
}

#include <thread>

void TestSpheroidalScatteringWideband_Threaded(int n_thread = 4) {
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


    /*std::filesystem::path folder_path = folder;
    if ( std::filesystem::exists(folder_path) ) {
        std::filesystem::remove_all(folder_path);
        //std::filesystem::remove(folder.c_str());
    }
    std::filesystem::create_directories(folder_path);*/

    struct stat info;

    if( stat( folder.c_str(), &info ) != 0 ) {
        printf( "cannot access %s\n", folder.c_str() );
    } else if( info.st_mode & S_IFDIR ) {
        printf( "%s is a directory\n", folder.c_str() );
        //delete
        std::string syscommand = std::string("rm ") + folder + "/*";
        std::system( syscommand.c_str() );
        std::cout << "content delete!" << std::endl;
    }

    int dir_err = mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        printf("Error creating directory %s!\n", folder.c_str());
    }

    std::vector<std::thread> thread_pool;
    for(int i = 0; i < n_thread; ++i) {
        thread_pool.emplace_back(SpheroidalScatteringWideband_ThreadJob, i, n_thread, tipRadius, length, f0, f1, n_f, folder);
    }

    for(int i = 0; i < n_thread; ++i) {
        thread_pool[i].join();
    }
}

#include <dirent.h>
//#include <algorithm>
#include "Utility.hpp"

void GenerateTemporalInterpolatorFromFrequencyInterpolator(std::string folder) {
    int ind_nf = folder.find("_nf=");
    int ind_f0 = folder.find("_f0=");
    int ind_f1 = folder.find("_f1=");
    assert(ind_nf != std::string::npos);
    assert(ind_f0 != std::string::npos);
    assert(ind_f1 != std::string::npos);

    const int n_f = std::stoi(folder.substr(ind_nf + 4, folder.size() - (ind_nf + 4)));
    double f_0 = std::stod(folder.substr(ind_f0 + 4, ind_f1 - (ind_f0 + 4))) * 1.0e12;
    double f_1 = std::stod(folder.substr(ind_f1 + 4, ind_nf - (ind_f1 + 4))) * 1.0e12;

    std::cout << "f_0: " << f_0 << "  f_1: " << f_1 << "  n_f: " << n_f << std::endl;

    // directory content
    DIR* dirp = opendir(folder.c_str());
    struct dirent * dp;
    std::vector<std::string> filenames;
    while ((dp = readdir(dirp)) != NULL) {
        std::string filename(dp->d_name);
        if (filename.find("interpMesh") != std::string::npos) {
            std::cout << filename << std::endl;
            filenames.push_back(filename);
        }
    }
    closedir(dirp);

    // get indices
    std::cout << "========== sorted =======" << std::endl;
    std::vector<int> filenames_inds(filenames.size(), -1);
    for(int i = 0; i < filenames.size(); ++i) {
        std::string& name = filenames[i];
        int ind = std::stoi(name.substr(name.find("_i=") + 3, name.find("_f=") - (name.find("_i=") + 3)));
        assert(ind >= 0 && ind < n_f);
        filenames_inds[ind] = i;
    }
    for(int i = 0; i < n_f; ++i) {
        std::cout << filenames[filenames_inds[i]] << std::endl;
    }

    //Fourier incident
    std::vector<double> t;
    std::vector<double> e_inc_t;
    GetIncidentField_Experiment(t, e_inc_t, true);
    std::size_t n_t = t.size();
    for(std::size_t i = 0; i < n_t; ++i) {
        //std::cout << "t: " << t[i] << " , E: " << e_inc_t[i] << std::endl;
    }

    std::vector<double> f;
    std::vector<std::complex<double>> e_inc_f;
    GetFourierTransform(t, e_inc_t, f_0, f_1, n_f, f, e_inc_f);
    for(std::size_t i = 0; i < n_f; ++i) {
        //std::cout << "f: " << f[i] << " , E(f): " << e_inc_f[i] << std::endl;
    }

    //Fourier inverse
    std::size_t n_pts_3d = 0;
    std::array<std::size_t, 3> numOfSamples;
    std::pair<double, double> eta_limits;
    std::pair<double, double> ksi_limits;
    std::pair<double, double> phi_limits;
    int ksiNonuniformPowScale;

    std::vector<std::vector<std::complex<double>>> E_eta_f_vec(n_f);
    std::vector<std::vector<std::complex<double>>> E_ksi_f_vec(n_f);
    std::vector<std::vector<std::complex<double>>> E_phi_f_vec(n_f);
    for(std::size_t i = 0; i < n_f; ++i) {
        SpheroidalInterpolator sphInterp;
        sphInterp.ReadMeshFromFile(folder + "/" + filenames[filenames_inds[i]]);
        E_eta_f_vec[i] = sphInterp.GetE_eta();
        E_ksi_f_vec[i] = sphInterp.GetE_ksi();
        E_phi_f_vec[i] = sphInterp.GetE_phi();
        if(i == 0) {
            n_pts_3d = E_eta_f_vec[i].size();
            auto& nSamples = sphInterp.GetNumOfSamples();
            numOfSamples[0] = nSamples[0];
            numOfSamples[1] = nSamples[1];
            numOfSamples[2] = nSamples[2];
            eta_limits = sphInterp.GetEtaLimits();
            ksi_limits = sphInterp.GetKsiLimits();
            phi_limits = sphInterp.GetPhiLimits();
            ksiNonuniformPowScale = sphInterp.GetKsiNonuniformScale();
            std::cout << "eta_limits: " << eta_limits.first << " " << eta_limits.second << std::endl;
            std::cout << "ksi_limits: " << ksi_limits.first << " " << ksi_limits.second << std::endl;
            std::cout << "phi_limits: " << phi_limits.first << " " << phi_limits.second << std::endl;
            assert(n_pts_3d > 0);
            assert(E_eta_f_vec[i].size() == n_pts_3d);
            assert(E_ksi_f_vec[i].size() == n_pts_3d);
            assert(E_phi_f_vec[i].size() == n_pts_3d);
        } else {
            assert(E_eta_f_vec[i].size() == n_pts_3d);
            assert(E_ksi_f_vec[i].size() == n_pts_3d);
            assert(E_phi_f_vec[i].size() == n_pts_3d);
        }
    }

    std::vector<std::vector<std::complex<double>>> E_eta_t_vec(n_t);
    std::vector<std::vector<std::complex<double>>> E_ksi_t_vec(n_t);
    std::vector<std::vector<std::complex<double>>> E_phi_t_vec(n_t);

    std::complex<double> _j(0.0, 1.0);
    for(std::size_t i = 0; i < n_t; ++i) {
        E_eta_t_vec[i].resize(n_pts_3d, 0.0);
        E_ksi_t_vec[i].resize(n_pts_3d, 0.0);
        E_phi_t_vec[i].resize(n_pts_3d, 0.0);
        for(std::size_t ind_pt = 0; ind_pt < n_pts_3d; ++ind_pt) {
            for(std::size_t j = 0; j < n_f; ++j) {
                double w_j = 2.0*M_PI*f[j];

                E_eta_t_vec[i][ind_pt] += 2.0*std::real(E_eta_f_vec[j][ind_pt] * e_inc_f[j] * std::exp(_j*w_j*t[i]));
                E_ksi_t_vec[i][ind_pt] += 2.0*std::real(E_ksi_f_vec[j][ind_pt] * e_inc_f[j] * std::exp(_j*w_j*t[i]));
                E_phi_t_vec[i][ind_pt] += 2.0*std::real(E_phi_f_vec[j][ind_pt] * e_inc_f[j] * std::exp(_j*w_j*t[i]));
            }
            E_eta_t_vec[i][ind_pt] *= (f[1] - f[0]);
            E_ksi_t_vec[i][ind_pt] *= (f[1] - f[0]);
            E_phi_t_vec[i][ind_pt] *= (f[1] - f[0]);

            //std::cout << E_ksi_t_vec[i][ind_pt] << std::endl;
        }
        //std::cout << t[i] << "=========================" << std::endl;
    }

    //write to disk
    for(std::size_t i = 0; i < n_t; ++i) {
        SpheroidalInterpolator sphInterp;
        sphInterp.SetEtaLimits(eta_limits);
        sphInterp.SetKsiLimits(ksi_limits);
        sphInterp.SetPhiLimits(phi_limits);
        sphInterp.SetKsiNonuniformScale(ksiNonuniformPowScale);
        sphInterp.SetNumOfSamples(numOfSamples[0], numOfSamples[1], numOfSamples[2]);

        auto& sphInt_e_eta = sphInterp.GetE_eta();
        auto& sphInt_e_ksi = sphInterp.GetE_ksi();
        auto& sphInt_e_phi = sphInterp.GetE_phi();

        sphInt_e_eta.resize(n_pts_3d);
        sphInt_e_ksi.resize(n_pts_3d);
        sphInt_e_phi.resize(n_pts_3d);

        auto& E_eta_t_vec_i = E_eta_t_vec[i];
        auto& E_ksi_t_vec_i = E_ksi_t_vec[i];
        auto& E_phi_t_vec_i = E_phi_t_vec[i];
        for(std::size_t j = 0; j < n_pts_3d; ++j) {
            sphInt_e_eta[j] = E_eta_t_vec_i[j];
            sphInt_e_ksi[j] = E_ksi_t_vec_i[j];
            sphInt_e_phi[j] = E_phi_t_vec_i[j];
        }


        std::string interpolator_fileSuffix = std::string("_")                                  \
                                        + "_i=" + boost::lexical_cast<std::string>(i)           \
                                        + "_t=" + boost::lexical_cast<std::string>(t[i]);
        sphInterp.WriteMeshToFile(folder + "/" + "ipt" + interpolator_fileSuffix + ".data");
    }
}

void PlotTemporllySavedFields(std::string folder) {
    double tipRadius= std::stod(folder.substr(folder.find("R=") + 2, folder.find("L=_") - (folder.find("R=") + 2))) * 1.0e-9;
    double length = std::stod(folder.substr(folder.find("_L=") + 3, folder.find("f0=_") - (folder.find("_L=") + 3))) * 1.0e-6;
    SpheroidScattering spheroid(tipRadius, length);
    spheroid.SetTemporalFieldInterpolators(folder);

    auto& t_arr = spheroid.GetTemporalSamples();
    int n_t = t_arr.size();

    double Dx = 4.0*tipRadius;
    double Dz = 4.0*tipRadius;
    int nx = 100;
    int nz = 100;

    std::vector<std::array<double, 3>> r_pts;   // spheroidal coords
    std::vector<std::array<int, 2>> r_inds;
    spheroid.GetXZGridPointsInSpheroidalCoords(Dx, Dz, nx, nz, r_pts, r_inds);

    std::vector<std::complex<double>> e_eta;
    std::vector<std::complex<double>> e_ksi;
    std::vector<std::complex<double>> e_phi;
    std::vector<std::complex<double>> e_x;
    std::vector<std::complex<double>> e_y;
    std::vector<std::complex<double>> e_z;
    for(int i = 0; i < n_t; ++i) {
        spheroid.GetTemporalFieldAtGridPoints_SpatialInterpolation(r_pts, e_eta, e_ksi, e_phi, i);
        spheroid.VectorTransformFromSpheroidToRect(r_pts, e_eta, e_ksi, e_phi, e_x, e_y, e_z);

        Matrix<std::complex<double>> e_x_mat(nx, nz);
        Matrix<std::complex<double>> e_y_mat(nx, nz);
        Matrix<std::complex<double>> e_z_mat(nx, nz);
        for(int i = 0; i < r_inds.size(); ++i) {
            e_x_mat(r_inds[i][0], r_inds[i][1]) = e_x[i];
            e_y_mat(r_inds[i][0], r_inds[i][1]) = e_y[i];
            e_z_mat(r_inds[i][0], r_inds[i][1]) = e_z[i];
        }

        std::string mat_fileSuffix = std::string("_")                          \
                                + "_i=" + boost::lexical_cast<std::string>(i)           \
                                + "_t=" + boost::lexical_cast<std::string>(t_arr[i])
                                + "_nx=" + boost::lexical_cast<std::string>(nx)
                                + "_nz=" + boost::lexical_cast<std::string>(nz)
                                + "_Dx=" + boost::lexical_cast<std::string>(Dx)
                                + "_Dz=" + boost::lexical_cast<std::string>(Dz)
                                + ".data";

        e_x_mat.WriteToFile(folder + "/" + "e_x" + mat_fileSuffix);
        e_y_mat.WriteToFile(folder + "/" + "e_y" + mat_fileSuffix);
        e_z_mat.WriteToFile(folder + "/" + "e_z" + mat_fileSuffix);
    }
}

#include "FowlerNordheimEmission.hpp"
#include "TipEmission.hpp"
#include "ChargedParticleTracer.hpp"

void TestSpheroidalScattering_Emission(std::string folder) {
    TipEmission emitter(folder);
    emitter.SetElectricFieldAmplitude(-3.0e7);
    emitter.SetMetalWorkFunction(4.5);
    double max_patch_area = (5.0*1.0e-9) * (5.0*1.0e-9);
    emitter.SubdevideSurface( 1.5 * emitter.GetSpheroid().GetTipRadius(), max_patch_area);
    auto n_particles = emitter.GetTotalNumberOfEmittedParticles();

    int n_total = 0;
    for(int j = 0; j < n_particles.size(); ++j) {
        n_total += n_particles[j];
        std::cout << j << "  - n_e: " << n_particles[j] << std::endl;
    }
    std::cout << "n_total : " << n_total << std::endl;

    int n_t = emitter.GetTimeSamples().size();
    for(int i = 0; i < n_t; ++i) {
        auto n_e = emitter.GetNumberOfEmittedParticles(i);
        std::cout << i << " - n_e(t) " << n_e[0] << std::endl;
    }

    auto& emissionPoints = emitter.GetEmissionPoints();
    auto& emissionPtNormals = emitter.GetEmissionPointNormals();
    auto& eFieldNormals = emitter.GetNormalEFields();
    auto& t_arr = emitter.GetTimeSamples();
    double vf = 1.0e4;
    double elec_charge = -PhysicalConstants_SI::electronCharge;

    std::string particle_folder = folder + "/particles";
    CreateFolderIfItDoesNotExists(particle_folder);

    ChargedParticleTracer etracer;
    etracer.ReserveMemory(n_total);
    std::vector<double> n_e(emissionPoints.size());
    for(int i = 0; i < n_t; ++i) {
        emitter.AddNumberOfEmittedParticles(i, n_e);
        std::cout << i << std::endl;
        for(int j = 0; j < n_e.size(); ++j) {
            if(n_e[j] > 1.0) {
                auto& a_n = emissionPtNormals[j];
                auto ej = eFieldNormals[j];
                std::array<double, 3> v0_e{vf*a_n[0], vf*a_n[1], vf*a_n[2]};
                std::array<double, 3> f0_e{elec_charge*ej*a_n[0], elec_charge*ej*a_n[1], elec_charge*ej*a_n[2]};
                etracer.AddParticle(n_e[j]*elec_charge, n_e[j]*PhysicalConstants_SI::electronMass, emissionPoints[j],
                                    v0_e, f0_e);
                n_e[j] = 0.0;
            }
        }

        if(i > 0) {
            double dt = t_arr[i] - t_arr[i - 1];
            emitter.GetElectricForce(etracer.GetCharges(), etracer.GetPositions(), etracer.GetForces(), i);
            etracer.UpdateParticles(dt);
        }
        etracer.SaveData(particle_folder, i, t_arr[i]);
    }
}








