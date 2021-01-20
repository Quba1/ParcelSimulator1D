#include "environment.h"
#include "parcel.h"
#include "thermodynamic_calc.h"
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

Parcel::Parcel()
{
    timeDeltaSquared = 0;
    timeDelta = 0;
    currentTimeStep = 0;
    ascentSteps = 0;
}

Parcel::Parcel(const std::map<std::string, std::string>& parcelConfiguration) :
    parcelConfiguration(parcelConfiguration),
    outputFileName(parcelConfiguration.at("output_filename"))
{
    calculateConstants();
    setupVariableFields();
    setInitialConditionsAndLocation();
}

void Parcel::calculateConstants()
{
    double period = std::stod(parcelConfiguration.at("period"));

    timeDelta = std::stod(parcelConfiguration.at("timestep"));
    timeDeltaSquared = timeDelta * timeDelta;

    ascentSteps = static_cast<size_t>(floor((period * 3600) / timeDelta) + 1); //including step zero
}

void Parcel::setupVariableFields()
{
    std::vector<double> holder(ascentSteps, -999.0);
    position = holder;
    velocity = holder;
    pressure = holder;
    temperature = holder;
    temperatureVirtual = holder;
    mixingRatio = holder;
    mixingRatioSaturated = holder;
}

void Parcel::setInitialConditionsAndLocation()
{
    //write initial conditions into parcel and convert to SI units

    //initial conditions from configuration
    position[0] = std::stod(parcelConfiguration.at("init_height"));
    velocity[0] = std::stod(parcelConfiguration.at("init_velocity"));
    temperature[0] = std::stod(parcelConfiguration.at("init_temp")) + 273.15;

    currentTimeStep = 0;

    currentLocation.position = position[0];
    currentLocation.updateSector();

    //intermediate variables initial conditions
    pressure[0] = Environment::getPressureAtLocation(currentLocation);

    mixingRatio[0] = calcMixingRatio((std::stod(parcelConfiguration.at("init_dewpoint")) + 273.15), pressure[0]);
    temperatureVirtual[0] = calcVirtualTemperature(temperature[0], mixingRatio[0]);
    mixingRatioSaturated[0] = calcMixingRatio(temperature[0], pressure[0]);
}

void Parcel::updateCurrentDynamicsAndPressure()
{
    currentLocation.position = position[currentTimeStep];
    currentLocation.updateSector();
    pressure[currentTimeStep] = Environment::getPressureAtLocation(currentLocation); //pressure of parcel always equalises with atmosphere
}

void Parcel::updateCurrentThermodynamicsAdiabatically(double lambda, double gamma)
{
    temperature[currentTimeStep] = calcTemperatureInAdiabat(pressure[currentTimeStep], gamma, lambda);
    mixingRatio[currentTimeStep] = mixingRatio[currentTimeStep - 1]; //mixing ratio is conservative during adiabatic ascent
    mixingRatioSaturated[currentTimeStep] = calcMixingRatio(temperature[currentTimeStep], pressure[currentTimeStep]);
    temperatureVirtual[currentTimeStep] = calcVirtualTemperature(temperature[currentTimeStep], mixingRatio[currentTimeStep]);
}

void Parcel::updateCurrentThermodynamicsPseudoadiabatically()
{
    mixingRatioSaturated[currentTimeStep] = calcMixingRatio(temperature[currentTimeStep], pressure[currentTimeStep]);
    mixingRatio[currentTimeStep] = mixingRatioSaturated[currentTimeStep];
    temperatureVirtual[currentTimeStep] = calcVirtualTemperature(temperature[currentTimeStep], mixingRatio[currentTimeStep]);
}

Parcel::Slice Parcel::getSlice(size_t stepsForwardFromCurrent)
{
    size_t timestep = currentTimeStep + stepsForwardFromCurrent;
    Parcel::Slice slice;

    slice.position = position[timestep];
    slice.pressure = pressure[timestep];
    slice.temperature = temperature[timestep];
    slice.mixingRatio = mixingRatio[timestep];
    slice.mixingRatioSaturated = mixingRatioSaturated[timestep];
    slice.temperatureVirtual = temperatureVirtual[timestep];

    return slice;
}