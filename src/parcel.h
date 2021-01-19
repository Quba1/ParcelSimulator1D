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
	size_t ascentSteps, currentTimeStep;

	void calculateParcelSize();
	void setupVariableFields();
	void setInitialConditionsAndLocation();

	void ascentAlongMoistAdiabat();
	void ascentAlongPseudoAdiabat();
	void ascentAlongDryAdiabat();

	double getCurrentBouyancyForce();
	void updateCurrentDynamics(DynamicPair dynamics);
	void updateCurrentThermodynamicAdiabatically(double lambda, double gamma);

public:
	DynamicScheme* dynamicScheme;
	PseudoAdiabaticScheme* pseudadiabaticScheme;

	Environment::Location currentLocation;

    std::unique_ptr<double[]> position, velocity, pressure, temperature, temperatureVirtual, mixingRatio, mixingRatioSaturated;
	const std::map<std::string, std::string> parcelConfiguration;

	Parcel(const std::map<std::string, std::string>& parcelConfiguration, DynamicScheme*& dynamicScheme, PseudoAdiabaticScheme*& pseudadiabaticScheme);
};

#endif
