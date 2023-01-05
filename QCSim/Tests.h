#pragma once

#include <complex>
#include <Eigen/eigen>

bool approxEqual(double val1, double val2, double err = 1E-10);
bool approxEqual(std::complex<double> val1, std::complex<double> val2, double err = 1E-10);
bool checkUnitary(const Eigen::MatrixXcd& m);

bool registerMeasurementsTests();
bool SimulationTests();

bool ParadoxesTests();
bool quantumAdderTests();

bool tests();
