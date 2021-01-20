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
	struct Slice
	{
		double position[2] = {};
		double velocity[2] = {};
		double pressure[2] = {};
		double temperature[2] = {};
		double temperatureVirtual[2] = {};
		double mixingRatio[2] = {};
		double mixingRatioSaturated[2] = {};

		Slice() {};
	};

	std::map<std::string, std::string> parcelConfiguration;
	std::string outputFileName;

	std::vector<double> position, velocity, pressure, temperature, temperatureVirtual, mixingRatio, mixingRatioSaturated;

	size_t ascentSteps, currentTimeStep;
	double timeDelta, timeDeltaSquared;
	Environment::Location currentLocation;

	Parcel();
	Parcel(const std::map<std::string, std::string>& parcelConfiguration);

	void updateCurrentDynamicsAndPressure();
	void updateCurrentThermodynamicsAdiabatically(double lambda, double gamma);
	void updateCurrentThermodynamicsPseudoadiabatically();

	Parcel::Slice getCurrentSlice();
};

#endif
