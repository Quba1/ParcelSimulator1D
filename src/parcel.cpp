#include "environment.h"
#include "parcel.h"
#include "thermodynamic_calc.h"
#include <cmath>
#include <map>
#include <memory>
#include <string>

Parcel::Parcel(const std::map<std::string, std::string>& parcelConfiguration, const Environment& parcelEnvironment):
    parcelConfiguration(parcelConfiguration),
    parcelEnvironment(parcelEnvironment)
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
    currentLocation.sector = updateSector(currentLocation);

    //intermediate variables initial conditions
    pressure[0] = parcelEnvironment.getPressureAtLocation(currentLocation);

    mixingRatio[0] = calcMixingRatio((std::stod(parcelConfiguration.at("init_dewpoint")) + 273.15), pressure[0]);
    temperatureVirtual[0] = calcVirtualTemperature(temperature[0], mixingRatio[0]);
    mixingRatioSaturated[0] = calcMixingRatio(temperature[0], pressure[0]);
}

Sector Parcel::updateSector(Location& location)
{
    //assuming sorted array of ascending values in heightField and location within bounds of heightField

    size_t nearestPoint;

    //check which sector boundary is closest
    double distToUpper = abs(location.position - parcelEnvironment.height[location.sector.upperBoundary]);
    double distToLower = abs(location.position - parcelEnvironment.height[location.sector.lowerBoundary]);

    if (distToLower == distToUpper)
    {
        //no need for sector correction
        return;
    }
    else if (distToLower < distToUpper)
    {
        if (location.sector.lowerBoundary == 0)
        {
            //location is in the lowest sector
            location.sector.upperBoundary = 1;
            return;
        }

        //move sector gradually down
        nearestPoint = location.sector.lowerBoundary;
        double dist = distToLower;
        double newDist = abs(location.position - parcelEnvironment.height[nearestPoint - 1]);

        while (newDist < dist)
        {
            dist = newDist;
            nearestPoint--;

            if (nearestPoint > 0)
            {
                newDist = abs(location.position - parcelEnvironment.height[nearestPoint - 1]);
            }
        }
    }
    else if (distToLower > distToUpper)
    {
        if (location.sector.upperBoundary == (parcelEnvironment.height.size() - 1))
        {
            //location is in the highest sector
            location.sector.lowerBoundary = parcelEnvironment.height.size() - 2;
            return;
        }

        //move sector gradually up
        nearestPoint = location.sector.upperBoundary;
        double dist = distToUpper;
        double newDist = abs(location.position - parcelEnvironment.height[nearestPoint + 1]);

        while (newDist < dist)
        {
            dist = newDist;
            nearestPoint++;

            if (nearestPoint < (parcelEnvironment.height.size() - 1))
            {
                newDist = abs(location.position - parcelEnvironment.height[nearestPoint + 1]);
            }
        }
    }

    double distToNearest = location.position - parcelEnvironment.height[nearestPoint];

    if (distToNearest >= 0)
    {
        location.sector.lowerBoundary = nearestPoint;
        location.sector.upperBoundary = nearestPoint + 1;
        return;
    }
    else
    {
        location.sector.lowerBoundary = nearestPoint - 1;
        location.sector.upperBoundary = nearestPoint;
        return;
    }

}
