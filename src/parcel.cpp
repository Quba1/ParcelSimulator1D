#include "environment.h"
#include "parcel.h"
#include "thermodynamic_calc.h"
#include "dynamic_scheme.h"
#include "pseudoadiabatic_scheme.h"
#include <cmath>
#include <map>
#include <memory>
#include <string>

Parcel::Parcel(const std::map<std::string, std::string>& parcelConfiguration, DynamicScheme*& dynamicScheme, PseudoAdiabaticScheme*& pseudadiabaticScheme) :
    parcelConfiguration(parcelConfiguration),
    dynamicScheme(dynamicScheme),
    pseudadiabaticScheme(pseudadiabaticScheme)
{
    calculateParcelSize();
    setupVariableFields();    
    setInitialConditionsAndLocation();
}

void Parcel::calculateParcelSize()
{
    double period = std::stod(parcelConfiguration.at("period"));
    double timeDelta = std::stod(parcelConfiguration.at("timestep"));
    
    ascentSteps = static_cast<size_t>(floor((period * 3600) / timeDelta) + 1); //including step zero
}

void Parcel::setupVariableFields()
{
    position = std::make_unique<double[]>(ascentSteps);
    velocity = std::make_unique<double[]>(ascentSteps);
    pressure = std::make_unique<double[]>(ascentSteps);
    temperature = std::make_unique<double[]>(ascentSteps);
    temperatureVirtual = std::make_unique<double[]>(ascentSteps);
    mixingRatio = std::make_unique<double[]>(ascentSteps);
    mixingRatioSaturated = std::make_unique<double[]>(ascentSteps);
}

void Parcel::setInitialConditionsAndLocation()
{
    //write initial conditions into parcel and convert to SI units

    //initial conditions from configuration
    position[0] = std::stod(parcelConfiguration.at("init_height"));
    velocity[0] = std::stod(parcelConfiguration.at("init_velocity"));
    temperature[0] = std::stod(parcelConfiguration.at("init_temp")) + 273.15;

    currentLocation.position = position[0];
    currentLocation.updateSector();

    currentTimeStep = 0;

    //intermediate variables initial conditions
    pressure[0] = Environment::getPressureAtLocation(currentLocation);

    mixingRatio[0] = calcMixingRatio((std::stod(parcelConfiguration.at("init_dewpoint")) + 273.15), pressure[0]);
    temperatureVirtual[0] = calcVirtualTemperature(temperature[0], mixingRatio[0]);
    mixingRatioSaturated[0] = calcMixingRatio(temperature[0], pressure[0]);
}

void Parcel::ascentAlongMoistAdiabat() 
{
    //calculate ascent constants
    double gamma = calcGamma(mixingRatio[currentTimeStep]);
    double lambda = calcLambda(temperature[currentTimeStep], pressure[currentTimeStep], gamma);
    
    //setup DynamicPair to be updated
    DynamicPair currentDynamicPair;
    currentDynamicPair.velocity = velocity[currentTimeStep];
    currentDynamicPair.location = currentLocation;

    //compute first timestep
    double bouyancyForce = getCurrentBouyancyForce();
    currentDynamicPair = dynamicScheme->computeFirstTimeStep(bouyancyForce, currentDynamicPair);
    
    currentTimeStep++;
    updateCurrentDynamics(currentDynamicPair);
    updateCurrentThermodynamicAdiabatically(lambda, gamma);

}

double Parcel::getCurrentBouyancyForce()
{
    double envPressure = Environment::getPressureAtLocation(currentLocation);
    double envTemperature = Environment::getTemperatureAtLocation(currentLocation);
    double envDewpoint = Environment::getDewpointAtLocation(currentLocation);

    double envMixingRatio = calcMixingRatio(envDewpoint, envPressure);
    double envVirtualTemperature = calcVirtualTemperature(envTemperature, envMixingRatio);

    return calcBouyancyForce(temperatureVirtual[currentTimeStep], envVirtualTemperature);
}

void Parcel::updateCurrentDynamics(DynamicPair dynamics)
{
    currentLocation = dynamics.location;
    position[currentTimeStep] = dynamics.location.position;
    velocity[currentTimeStep] = dynamics.velocity;
}

void Parcel::updateCurrentThermodynamicAdiabatically(double lambda, double gamma)
{
    pressure[currentTimeStep] = Environment::getPressureAtLocation(currentLocation); //pressure of parcel always equalises with atmosphere
    mixingRatio[currentTimeStep] = mixingRatio[currentTimeStep - 1]; //mixing ratio is conservative during adiabatic ascent
    temperature[currentTimeStep] = calcTemperatureInAdiabat(pressure[currentTimeStep], gamma, lambda);
    mixingRatioSaturated[currentTimeStep] = calcMixingRatio(temperature[currentTimeStep], pressure[currentTimeStep]);
    temperatureVirtual[currentTimeStep] = calcVirtualTemperature(temperature[currentTimeStep], mixingRatio[currentTimeStep]);
}

void Parcel::ascentAlongDryAdiabat()
{

}

void Parcel::ascentAlongPseudoAdiabat()
{

}