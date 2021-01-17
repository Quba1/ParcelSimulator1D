#ifndef PARCEL_H
#define PARCEL_H

#include "environment.h"
#include <map>
#include <memory>
#include <string>

class Parcel
{
private:
	double timeDelta, timeDeltaSquared;
	size_t ascentSteps, currentTimeStep  = 0;

	void calculateConstants();
	void setupVariableFields();
	void setInitialConditionsAndLocation();

	Sector updateSector(Location& location);

public:
	Environment parcelEnvironment;
	Location currentLocation;
    std::unique_ptr<double[]> position, velocity, pressure, temperature, temperatureVirtual, mixingRatio, mixingRatioSaturated;
	const std::map<std::string, std::string> parcelConfiguration;

	Parcel(const std::map<std::string, std::string>& parcelConfiguration, const Environment& parcelEnvironment);
};

#endif
