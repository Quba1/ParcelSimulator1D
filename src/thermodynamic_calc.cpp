#include "thermodynamic_calc.h"
#include <cmath>

double calcVapourPressure(double temperature, double pressure)
{
   //input in K & Pa; output in ratio of Pa
   temperature -= 273.15; //convert to C
   pressure /= 100.0; //convert to hPa
   
   const double a = 6.1121;
   const double b = 18.729;
   const double c = 257.87;
   const double d = 227.3;
   
   const double A = 0.00072; 
   const double B = 0.0000032;
   const double C = 0.00000000059;
   
   const double e = a * exp(((b - (temperature / d)) * temperature)/(temperature + c));
   const double f = 1.0 + A + (pressure * (B + (C * temperature * temperature)));
   
   return (e * f) * 100.0; //return in Pa  
}

double calcMixingRatio(double temperature, double pressure)
{
    //input in K & Pa; output in ratio of kg/kg
    //function for caluclating both mixing ratio and saturation mixing ratio

    //first calculate (saturation) vapour pressure
    double wvpres = calcVapourPressure(temperature, pressure);

    //second return calculated mixing ratio
    return EPSILON * (wvpres / (pressure - wvpres));
}

double calcVirtualTemperature(double temperature, double mixRatio)
{
    //input in K & Pa/Pa, output in K
    return temperature * ((1.0 + (mixRatio / EPSILON)) / (1.0 + mixRatio));
}

double calcTemperatureInAdiabat(double pressure, double gamma, double lambda)
{
    //input in K; output in K
    return pow(lambda / pow(pressure, 1.0 - gamma), 1.0 / gamma);
}

double calcBouyancyForce(double parcelTv, double envTv)
{
    return G * ((parcelTv - envTv) / envTv);
}

double calcWBPotentialTemperature(double temperature, double mixingRatio, double satMixingRatio, double pressure)
{
    //input in K & kg/kg & Pa; output in K
    double vapourPressure = calcVapourPressure(temperature, pressure);
    double dryTemperature = temperature * pow(100000.0 / ((pressure - vapourPressure)), 0.2854);
    double equiTemperature = dryTemperature * pow(mixingRatio / satMixingRatio, -(0.2854 * (mixingRatio / EPSILON))) * exp((2555000.0 * mixingRatio) / (C_P * temperature)); //Bryan (2008)
    double wetTemperature = (45.114 - (51.489 * pow(273.15 / equiTemperature, 3.504))) + 273.15;

    return wetTemperature;
}

double calcGamma(double mixingRatio)
{
    double gamma = (C_P * ((1.0 + (mixingRatio * (C_PV / C_P))) / (1.0 + mixingRatio))) / (C_V * ((1.0 + (mixingRatio * (C_VV / C_V))) / (1.0 + mixingRatio))); //Bailyn (1994)
    return gamma;
}

double calcLambda(double temperature, double pressure, double gamma)
{
    double lambda = pow(pressure, 1.0 - gamma) * pow(temperature, gamma);
    return lambda;
}
