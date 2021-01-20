#ifndef PARCEL_H
#define PARCEL_H

#include "environment.h"
#include "dynamic_scheme.h"
#include "pseudoadiabatic_scheme.h"
#include <map>
#include <memory>
#include <string>

class Parcel
{
private:
	void calculateConstants();
	void setupVariableFields(); 
	void setInitialConditionsAndLocation();

public:
	std::map<std::string, std::string> parcelConfiguration;
	std::vector<double> position, velocity, pressure, temperature, temperatureVirtual, mixingRatio, mixingRatioSaturated;

	size_t ascentSteps, currentTimeStep;
	double timeDelta, timeDeltaSquared;
	Environment::Location currentLocation;

	Parcel();
	Parcel(const std::map<std::string, std::string>& parcelConfiguration);

	void updateCurrentDynamics();
	void updateCurrentThermodynamicsAdiabatically(double lambda, double gamma);
	void updateCurrentThermodynamicallyPseudoadiabatically();
};

#endif
