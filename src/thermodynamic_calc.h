#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

#include <cmath>

double const EPSILON = 0.62194572947593730; //ratio of g/mol over g/mol; non SI, but all formulae are adapted to it
double const G = 9.80665; //gravitational acceleration in m*s^-2

double calcMixingRatio(double temperature, double pressure);

double calcVirtualTemperature(double temperature, double mixRatio);

double calcTemperatureInAdiabat(double pressure, double gamma, double lambda);

double calcBouyancyForce(double parcelTv, double envTv);

double calcWBPotentialTemperature(double temperature, double mixingRatio, double satMixingRatio, double pressure);

double calcGamma(double mixingRatio);

double calcLambda(double temperature, double pressure, double gamma);

#endif
