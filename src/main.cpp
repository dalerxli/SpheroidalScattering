
#include <iostream>
#include "stdlib.h"
#include "time.h"

#include "SpheroidalFunc.hpp"
#include "SpheroidalIntegrals.hpp"
#include "SpheroidScattering.hpp"
#include "Tests.hpp"
#include "TestInterpolation.hpp"
#include "TestSpheroidalScattering.hpp"
#include "TestSpheroidalScatteringWidebandMPI.hpp"

int main() {

    //TestExpansionCoefficients(3, 5, 1.0);

    //TestRadialProlate(3, 5, 2.0, 1.5);
    //TestRadialOblate(3, 5, 2.0, 1.5);

    //TestAngularFirstKind(3, 5, 2.0, 0.5);


    //std::cout << GetInt_dxSmp2_mpnp2_Sm_mpN_1mx2(1.1, 3, 1, 2) << std::endl;

    //std::cout << GetInt_dxS0_n_S1_1pN_sqrt_x_sqrt_1mx2(1.2, 7, 3) << std::endl;

    //lapack_int info = TestLapack_dgesv();
    //std::cout << "info: " << info << std::endl;

    //info = TestLapack_zgesvx();
    //std::cout << "info: " << info << std::endl;

    //std::cout << "Testing the Matsolver: " << std::endl;
    //TestMatSolver();
    //TestMatrixRead();
    //TestUmfpack();
    //TestUmfpackComplex();
    //TestOMP();

    //std::cout << "=================================================== " << std::endl;
    //std::cout << "Testing spheroidal scattering: " << std::endl;
    //TestSpheroidalScattering();
    //TestSpheroidalScatteringWideband();
    //TestSpheroidalScatteringWidebandMPI();    // bug: wrong results when optimized
    //TestSpheroidalScatteringWideband_Threaded(4);
    //GenerateTemporalInterpolatorFromFrequencyInterpolator("out/R=50_L=300_f0=0.01_f1=3_nf=200");
    PlotTemporllySavedFields("out/R=50_L=300_f0=0.01_f1=3_nf=200");

    //testInterpolation3D();
    //TestSpheroidVectorTransforms();

    return 0;
}

