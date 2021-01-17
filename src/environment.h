#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct Sector
{
    size_t upperBoundary, lowerBoundary;

	Sector(size_t upperBoundary, size_t lowerBoundary);
	Sector();
};

struct Location
{
	double position;
	Sector sector;

	Location(double position, Sector sector);
	Location();
};

class Environment
{
private:
	void importDataFrom(std::ifstream& file);
	double getInterpolatedValueofFieldAtLocation(const std::vector<double>& variableField, const Location& location);

public:
	std::vector<double> height, pressure, temperature, dewpoint;

	Environment(std::string configurationFileName);

	double getPressureAtLocation(const Location& location);
};

#endif