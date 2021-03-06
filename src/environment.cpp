#include "thermodynamic_calc.h"
#include "environment.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

std::vector<double> Environment::height;
std::vector<double> Environment::pressure;
std::vector<double> Environment::temperature;
std::vector<double> Environment::dewpoint;
double Environment::highestPoint;

Sector::Sector()
{
    lowerBoundary = 0;
    upperBoundary = 1;
}

Environment::Location::Location()
{
    sector = Sector();
    position = 0.0;
}

Environment::Environment(std::string configurationFileName)
{
    std::ifstream configurationFile(configurationFileName);
    importDataFrom(configurationFile);
    configurationFile.close();

    highestPoint = height[height.size() - 1];
}

void Environment::importDataFrom(std::ifstream& file)
{
    std::string line;

    //dummy getlines for skipping the header lines
    getline(file, line);
    getline(file, line);

    //read variables line by line and push into separate vectors
    while (getline(file, line))
    {
        std::stringstream lineStream(line);
        std::string var;

        getline(lineStream, var, ';');
        height.push_back(stod(var));
        getline(lineStream, var, ';');
        pressure.push_back(stod(var));
        getline(lineStream, var, ';');
        temperature.push_back(stod(var));
        getline(lineStream, var, ';');
        dewpoint.push_back(stod(var));
    }
}

double Environment::getInterpolatedValueofFieldAtLocation(const std::vector<double>& variableField, const Location& location)
{
    //do linear interpolation of the field within the sector
    double b = (variableField[location.sector.upperBoundary] - variableField[location.sector.lowerBoundary]) / (height[location.sector.upperBoundary] - height[location.sector.lowerBoundary]);
    double a = variableField[location.sector.lowerBoundary] - (height[location.sector.lowerBoundary] * b);

    //caluclate and return pressure
    return (a + (b * location.position));
}

double Environment::getPressureAtLocation(const Location& location)
{
    //input in m; output in Pa
    double value = getInterpolatedValueofFieldAtLocation(pressure, location);
    
    return value * 100.0;

}

double Environment::getTemperatureAtLocation(const Location& location)
{
    //input in m; output in K
    double value = getInterpolatedValueofFieldAtLocation(temperature, location);

    return value + 273.15;

}

double Environment::getDewpointAtLocation(const Location& location)
{
    //input in m; output in K
    double value = getInterpolatedValueofFieldAtLocation(dewpoint, location);

    return value + 273.15;

}

double Environment::getVirtualTemperatureAtLocation(const Location& location)
{
    double press = getPressureAtLocation(location);
    double temp = getTemperatureAtLocation(location);
    double dwpt = getDewpointAtLocation(location);

    double mixr = calcMixingRatio(dwpt, press);
    return calcVirtualTemperature(temp, mixr);
}

void Environment::Location::updateSector()
{
    //assuming sorted array of ascending values in heightField and location within bounds of heightField

    size_t nearestPoint;

    //check which sector boundary is closest
    double distToUpper = abs(position - height[sector.upperBoundary]);
    double distToLower = abs(position - height[sector.lowerBoundary]);

    if (distToLower == distToUpper)
    {
        //no need for sector correction
        return;
    }
    else if (distToLower < distToUpper)
    {
        if (sector.lowerBoundary == 0)
        {
            //location is in the lowest sector
            sector.upperBoundary = 1;
            return;
        }

        //move sector gradually down
        nearestPoint = sector.lowerBoundary;
        double dist = distToLower;
        double newDist = abs(position - height[nearestPoint - 1]);

        while (newDist < dist)
        {
            dist = newDist;
            nearestPoint--;

            if (nearestPoint > 0)
            {
                newDist = abs(position - height[nearestPoint - 1]);
            }
            else
            {
                sector.lowerBoundary = 0;
                sector.upperBoundary = 1;
                return;
            }
        }
    }
    else if (distToLower > distToUpper)
    {
        if (sector.upperBoundary == (height.size() - 1))
        {
            //location is in the highest sector
            sector.lowerBoundary = height.size() - 2;
            return;
        }

        //move sector gradually up
        nearestPoint = sector.upperBoundary;
        double dist = distToUpper;
        double newDist = abs(position - height[nearestPoint + 1]);

        while (newDist < dist)
        {
            dist = newDist;
            nearestPoint++;

            if (nearestPoint < (height.size() - 1))
            {
                newDist = abs(position - height[nearestPoint + 1]);
            }
            else
            {
                sector.lowerBoundary = height.size() - 2;
                sector.upperBoundary = height.size() - 1;
                return;
            }
        }
    }

    double distToNearest = position - height[nearestPoint];

    if (distToNearest >= 0)
    {
        sector.lowerBoundary = nearestPoint;
        sector.upperBoundary = nearestPoint + 1;
        return;
    }
    else
    {
        sector.lowerBoundary = nearestPoint - 1;
        sector.upperBoundary = nearestPoint;
        return;
    }
}