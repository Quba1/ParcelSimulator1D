#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

#include <cmath>

double const EPSILON = 0.62194572947593730; //ratio of g/mol over g/mol; non SI, but all formulae are adapted to it
double const G = 9.80665; //gravitational acceleration in m*s^-2
double const C_P = 1005.7; //J kg^-1 K^-1 specific heat capacity of dry air at constant pressure
double const C_V = 718; //J kg^-1 K^-1 specific heat capacity of dry air at constant volume
double const C_PV = 1870; //J kg^-1 K^-1 specific heat capacity of water vapour at constant pressure
double const C_VV = 1410; //J kg^-1 K^-1 specific heat capacity of water vapour at constant volume

double calcMixingRatio(double temperature, double pressure);

double calcVirtualTemperature(double temperature, double mixRatio);

double calcTemperatureInAdiabat(double pressure, double gamma, double lambda);

double calcBouyancyForce(double parcelTv, double envTv);

double calcWBPotentialTemperature(double temperature, double mixingRatio, double satMixingRatio, double pressure);

double calcGamma(double mixingRatio);

double calcLambda(double temperature, double pressure, double gamma);

#endif
