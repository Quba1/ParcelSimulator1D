#include "environment.h"
#include "parcel.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

Sector::Sector()
{
    lowerBoundary = 0;
    upperBoundary = 1;
}

Sector::Sector(size_t upperBoundary, size_t lowerBoundary) :
    upperBoundary(upperBoundary),
    lowerBoundary(lowerBoundary)
{}

Location::Location()
{
    position = 0;
    sector = Sector();
}

Location::Location(double position, Sector sector) :
    position(position),
    sector(sector)
{}

Environment::Environment(std::string configurationFileName)
{
    std::ifstream configurationFile(configurationFileName);
    importDataFrom(configurationFile);
    configurationFile.close();
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
