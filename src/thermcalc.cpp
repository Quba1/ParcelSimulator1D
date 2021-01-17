#include "thermcalc.h"
#include <cmath>

double calcMixingRatio(double temperature, double pressure)
{
    //input in K & Pa; output in ratio of kg/kg
    //function for caluclating both mixing ratio and saturation mixing ratio

    //first calculate (saturation) vapour pressure
    temperature = temperature - 273.15;
    double wvpres = (6.1121 * exp((17.502 * temperature) / (240.97 + temperature))) * 100; //Buck (1981)

    //second return calculated mixing ratio
    return EPSILON * (wvpres / (pressure - wvpres));
}

double calcVirtualTemperature(double temperature, double mixRatio)
{
    //input in K & Pa/Pa, output in K
    return temperature * ((1 + (mixRatio / EPSILON)) / (1 + mixRatio));
}

double calcTemperatureInAdiabat(double pressure, double gamma, double lambda)
{
    //input in K; output in K
    return pow(lambda / pow(pressure, 1 - gamma), 1 / gamma);
}

double calcBouyancyForce(double parcelTv, double envTv)
{
    return G * ((parcelTv - envTv) / envTv);
}

double calcVelocity(double position1, double position2, double dt)
{
    return (position2 - position1) / dt;
}

double calcWBPotentialTemperature(double temperature, double mixingRatio, double satMixingRatio, double pressure)
{
    //input in K & kg/kg & Pa; output in K
    double vapourPressure = (6.1121 * exp((17.502 * (temperature - 273.15)) / (240.97 + (temperature - 273.15)))) * 100.0;
    double dryTemperature = temperature * pow(100000.0 / ((pressure - vapourPressure)), 0.2854);
    double equiTemperature = dryTemperature * pow(mixingRatio / satMixingRatio, -(0.2854 * (mixingRatio / EPSILON))) * exp((2555000.0 * mixingRatio) / (1005.7 * temperature)); //Bryan (2008)
    double wetTemperature = (45.114 - (51.489 * pow(273.15 / equiTemperature, 3.504))) + 273.15;

    return wetTemperature;
}
