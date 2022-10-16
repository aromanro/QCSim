#pragma once

#include <complex>

bool approxEqual(double val1, double val2, double err = 1E-10);
bool approxEqual(std::complex<double> val1, std::complex<double> val2, double err = 1E-10);

bool registerMeasurementsTests();

bool tests();
