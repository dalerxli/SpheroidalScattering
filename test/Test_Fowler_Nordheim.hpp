

#include "FowlerNordheimEmission.hpp"

void test_Fowler_Nordheim() {
    double eField_si = 3.0e9; // v/m
    double workFunction_eV = 4.5;
    FowlerNordheimEmission fn(eField_si, workFunction_eV);

    int ny = 20;
    double dy = 1.0/ny;
    for(int i = 0; i < ny; ++i) {
        double y = (double)i * dy;
        double vy = fn.GetParameterVy(y);

        std::cout << "y : " << y << " , v(y) : " << vy <<std::endl;
    }

    double numOfElectrons = fn.GetNumberOfEmittedElectrons(1.0e-7, 1.0e-12);
    std::cout << "numOfElectrons : " << numOfElectrons << std::endl;

    std::cout << "-----------------------------------" << std::endl;
    double eField_si_start = 1.0e7; // v/m
    double eField_si_end = 1.0e10; // v/m
    int n_field = 100;
    for(int i = 0; i < n_field; ++i) {
        double eField_si_i = eField_si_start + (double)i/n_field * (eField_si_end - eField_si_start);
        fn.SetEfield(eField_si_i);
        numOfElectrons = fn.GetNumberOfEmittedElectrons(1.0e-7, 1.0e-12);
        std::cout << "eField_si_i: " << eField_si_i << " , numOfElectrons: " << numOfElectrons << std::endl;
    }

};
