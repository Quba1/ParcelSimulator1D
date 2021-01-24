#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct Sector
{
    size_t upperBoundary, lowerBoundary;

	Sector();
};


class Environment
{

public:
	struct Location
	{
		double position;
		Sector sector;

		Location();
		void updateSector();
	};

	static double highestPoint;

	static std::vector<double> height, pressure, temperature, dewpoint;

	Environment(std::string configurationFileName);
	static double getPressureAtLocation(const Location& location);
	static double getTemperatureAtLocation(const Location& location);
	static double getDewpointAtLocation(const Location& location);
	static double getVirtualTemperatureAtLocation(const Location& location);

private:
	void importDataFrom(std::ifstream& file);
	static double getInterpolatedValueofFieldAtLocation(const std::vector<double>& variableField, const Location& location);
};

#endif