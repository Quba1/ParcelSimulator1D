#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

//externaly defined constants
double const G = 9.80665; //gravitational acceleration in m*s^-2
double const R = 8.31446261815324; //J K^-1 mol^-1 universal gas constant

double const M_D = 0.0289644; // kg mol^-1 molar mass of dry air (ECMWF, 2020)
double const M_V = 0.0180152833; //kg mol^-1 molar mass of water vapour

double const C_P = 1004.709; //J kg^-1 K^-1 specific heat capacity of dry air at constant pressure (ECMWF, 2020)
double const C_V = 717.6493; //J kg^-1 K^-1 specific heat capacity of dry air at constant volume (ECMWF, 2020)
double const C_PV = 1846.1; //J kg^-1 K^-1 specific heat capacity of water vapour at constant pressure (ECMWF, 2020)
double const C_VV = 1384.575; //J kg^-1 K^-1 specific heat capacity of water vapour at constant volume (ECMWF, 2020)

double const L_V = 2500800.0; //J kg^1 latent heat of vapourization of water (ECMWF, 2020)

//pre-calculated constants
double const EPSILON = M_V / M_D; //ratio of molar masses of dry air and water vapour
double const R_D = R / M_D; //specific gas constant for dry air

double calcVapourPressure(double temperature, double pressure);

double calcMixingRatio(double temperature, double pressure);

double calcVirtualTemperature(double temperature, double mixRatio);

double calcTemperatureInAdiabat(double pressure, double gamma, double lambda);

double calcBouyancyForce(double parcelTv, double envTv);

double calcWBPotentialTemperature(double temperature, double mixingRatio, double satMixingRatio, double pressure);

double calcGamma(double mixingRatio);

double calcLambda(double temperature, double pressure, double gamma);

#endif
