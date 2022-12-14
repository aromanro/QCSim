#pragma once

#include <complex>
#include <Eigen/eigen>
#include <random>

bool approxEqual(double val1, double val2, double err = 1E-10);
bool approxEqual(std::complex<double> val1, std::complex<double> val2, double err = 1E-10);
bool checkUnitary(const Eigen::MatrixXcd& m);

bool registerMeasurementsTests();
bool SimulationTests();

bool ParadoxesTests();
bool quantumAdderTests();

bool ErrorCorrectionTests();

bool tests();

#ifndef TESTS_CPP_

extern std::random_device rd;
extern std::mt19937 gen;

extern std::uniform_int_distribution<int> dist_bool;
extern std::uniform_real_distribution<double> dist_ampl;

#endif
