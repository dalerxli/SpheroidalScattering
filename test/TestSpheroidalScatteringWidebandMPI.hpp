
#include <mpi.h>

#include <utility>
#include <string>
#include <vector>

#include "boost/lexical_cast.hpp"

#include <sys/types.h>  // directory check
#include <sys/stat.h>
//#include <filesystem>

#include "SpheroidScattering.hpp"


void GetStartAndEndFrequencyIndices(int n_f, int numOfProcesses, int processRank, int& ind_f_start, int& ind_f_end) {
    ind_f_start = processRank * n_f / numOfProcesses;
    ind_f_end =  (processRank + 1) * n_f / numOfProcesses;
}

void GetStartAndEndFrequencyIndices(int n_f, int numOfProcesses, int processRank, std::vector<int>& inds_f) {
    inds_f.clear();
    int i = 0;
    while(i < n_f) {
        for(int j = 0; j < numOfProcesses; ++j) {
            if(j == processRank) {
                inds_f.push_back(i);
            }
            ++i;

            if(i >= n_f) {
                break;
            }
        }
    }
}

void TestSpheroidalScatteringWidebandMPI() {

    int numOfProcesses;
    int processRank;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int processorNameLen;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Get_processor_name(processorName, &processorNameLen);

    double tipRadius= 50.0 * 1.0e-9;
    double length = 1500.0 * 1.0e-6;
    double f0 = 0.01 * 1.0e12;
    double f1 = 3.0 * 1.0e12;
    int n_f = 600;

    std::string folder("out/");
    folder += std::string("R=") + boost::lexical_cast<std::string>(tipRadius/1.0e-9)
                        + "_L=" + boost::lexical_cast<std::string>(length/1.0e-6)
                        + "_f0=" + boost::lexical_cast<std::string>(f0/1.0e12)
                        + "_f1=" + boost::lexical_cast<std::string>(f1/1.0e12)
                        + "_nf=" + boost::lexical_cast<std::string>(n_f)
                        ;


    if(processRank == 0) {
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
    }

    MPI_Barrier(MPI_COMM_WORLD);

    folder += "/";

    SpheroidScattering spheroid(tipRadius, length);

    std::vector<int> inds_f;
    GetStartAndEndFrequencyIndices(n_f, numOfProcesses, processRank, inds_f);

    for(int i : inds_f) {

        double freq =  f0 + (double)i / (n_f - 1) * (f1 - f0);
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

        bool plotFields = false;
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
        } else{
            spheroid.SetupFieldInterpolator(alpha, beta, gamma);
        }

        std::string interpolator_fileSuffix = std::string("_")                                  \
                                        + "_i=" + boost::lexical_cast<std::string>(i)           \
                                        + "_f=" + boost::lexical_cast<std::string>(freq);
        auto& sphInterpolator = spheroid.GetSpheroidalInterpolator();
        sphInterpolator.WriteMeshToFile(folder + "interpMesh" + interpolator_fileSuffix + ".data");

    }

    MPI_Finalize();
}


#include "FowlerNordheimEmission.hpp"
#include "TipEmission.hpp"
#include "ChargedParticleTracer.hpp"

void TestSpheroidalScattering_Emission_MPI(std::string folder) {
    int numOfProcesses;
    int processRank;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int processorNameLen;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Get_processor_name(processorName, &processorNameLen);

    TipEmission emitter(folder);
    emitter.SetElectricFieldAmplitude(-2.0e7);
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
    if(processRank == 0) {
        CreateFolderIfItDoesNotExists(particle_folder, true /*delete_content*/);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    ChargedParticleTracer etracer;
    etracer.ReserveMemory(n_total);
    std::vector<double> n_e(emissionPoints.size());

    int n_time_subdivisions = 10;

    std::string filename_mc = std::string("mpi_masscharge_") + boost::lexical_cast<std::string>(processRank) + "_";
    std::string filename_pvm = std::string("mpi_posvelmom_") + boost::lexical_cast<std::string>(processRank) + "_";

    int emitted_particle_index = 0;

    for(int i = 0; i < n_t; ++i) {
        for(int i_sub = 0; i_sub < n_time_subdivisions; ++i_sub) {
            emitter.AddNumberOfEmittedParticles(i, n_e, i_sub, n_time_subdivisions);
            std::cout << i << std::endl;
            for(int j = 0; j < n_e.size(); ++j) {
                if(n_e[j] > 1.0) {
                    auto& a_n = emissionPtNormals[j];
                    auto ej = eFieldNormals[j];
                    std::array<double, 3> v0_e{vf*a_n[0], vf*a_n[1], vf*a_n[2]};
                    std::array<double, 3> f0_e{elec_charge*ej*a_n[0], elec_charge*ej*a_n[1], elec_charge*ej*a_n[2]};

                    //assign the particle to only one of the processes
                    if(emitted_particle_index % numOfProcesses  == processRank) {
                        etracer.AddParticle(n_e[j]*elec_charge, n_e[j]*PhysicalConstants_SI::electronMass, emissionPoints[j],
                                            v0_e, f0_e);
                    }

                    n_e[j] = 0.0;
                    emitted_particle_index++;
                }
            }

            if(i > 0) {
                double dt = (t_arr[i] - t_arr[i - 1])/n_time_subdivisions;
                emitter.GetElectricForce(etracer.GetCharges(), etracer.GetPositions(), etracer.GetForces(), i, i_sub, n_time_subdivisions);
                etracer.UpdateParticles(dt);
            }
        }
        etracer.SaveData(particle_folder, i, t_arr[i], filename_mc, filename_pvm);
    }

    MPI_Finalize();
}





