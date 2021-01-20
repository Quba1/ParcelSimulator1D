#include "thermodynamic_calc.h"
#include "parcel.h"
#include "pseudoadiabatic_scheme.h"

double NumericalPseudoadiabat::calculateCurrentPseudoadiabaticTemperature(Parcel::Slice currentParcelSlice, double deltaPressure, double WetBulbTheta)
{
    //input in Pa & K; output in K
    //Bakhshaii & Stull (2013)

    //necesary unit conversions and declarations
    double pressure = (currentParcelSlice.pressure + deltaPressure) / 1000.0;
    double wbTemp = WetBulbTheta - 273.15;

    //compute temperature
    if (wbTemp > -30.0 && wbTemp <= 4.0)
    {
        double g1 = -20.3313 - (0.0253 * pressure);
        double g2 = sin(pow(wbTemp + pressure, 0.5)) + (wbTemp / pressure) + pressure - 2.8565;
        double g3 = cos(19.6836 + pow(1.0 + exp(-wbTemp), -(1.0 / 3.0)) + (pressure / 15.0252));
        double g4 = (4.4653 * sin(pow(pressure, 0.5))) - 71.9358;
        double g5 = pow(exp(wbTemp - (2.71828 * cos(pressure / 18.5219))), (1.0 / 6.0));
        double g6 = wbTemp - sin(pow(pressure + wbTemp + atan(wbTemp) + 6.6165, 0.5));

        return g1 + g2 + g3 + g4 + g5 + g6 + 273.15;

    }
    else if (wbTemp > 4.0 && wbTemp <= 21.0)
    {
        double g1 = -9.6285 + cos(log(atan(atan(exp((-9.2121 * wbTemp) / pressure)))));
        double g2 = wbTemp - ((19.9563 / pressure) * atan(wbTemp)) + (pow(wbTemp, 2.0) / (5.47162 * pressure));
        double g3 = sin(log(8.0 * pow(pressure, 3.0))) * log(2.0 * pow(pressure, 1.5));
        double g4 = wbTemp + (((pressure * wbTemp) - pressure + wbTemp) / (pressure - 190.2578));
        double g5 = pressure - ((pressure - 383.0292) / ((15.4014 * pressure) - pow(pressure, 2.0)));
        double g6 = ((1.0 / 3.0) * log(339.0316 - pressure)) + atan(wbTemp - pressure + 95.9839);
        double g7 = -(log(pressure) * ((298.2909 + (16.5109 * pressure)) / (pressure - 2.2183)));

        return g1 + g2 + g3 + g4 + g5 + g6 + g7 + 273.15;
    }
    else if (wbTemp > 21.0 && wbTemp < 45.0)
    {
        double g1 = 0.3919 * pow(wbTemp, 7 / 3) * pow(pressure * (pressure + 15.8148), -1.0);
        double g2 = (19.9724 + (797.7921 / pressure)) * sin(-19.9724 / wbTemp);
        double g3 = pow(log(-3.927765 + wbTemp + pressure) * cos(log(wbTemp + pressure)), 3.0);
        double g4 = pow(exp(pow(wbTemp + pow(1.0 + exp(-pressure), -1.0), 0.5) - 1.5603), 0.5);
        double g5 = pow(pressure + wbTemp, 0.5) * exp(atan((pressure + wbTemp) / 7.9081));
        double g6 = ((pressure / pow(wbTemp, 2)) * std::min(9.6112, pressure - wbTemp)) - 13.73;
        double g7 = sin(pow(sin(std::min(pressure, 17.3170)), 3.0) - pow(pressure, 0.5) + (25.5113 / wbTemp));

        return g1 + g2 + g3 + g4 + g5 + g6 + g7 + 273.15;
    }

    //if theta out of range
    return -999.0;
}

double FiniteDifferencePseudoadiabat::calculateCurrentPseudoadiabaticTemperature(Parcel::Slice currentParcelSlice, double deltaPressure, double WetBulbTheta)
{
    return 0.0;
}

double RungeKuttaPseudoadiabat::calculateCurrentPseudoadiabaticTemperature(Parcel::Slice currentParcelSlice, double deltaPressure, double WetBulbTheta)
{
    return 0.0;
}
