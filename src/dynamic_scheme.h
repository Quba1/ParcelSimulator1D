#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "environment.h"
#include "parcel.h"
#include "pseudoadiabatic_scheme.h"
#include <memory>

class DynamicScheme
{
public:
	virtual void runSimulationOn(Parcel& passedParcel, size_t pseudoadiabaticSchemeID) = 0;
};

class FiniteDifferenceDynamics : public DynamicScheme
{
private:
	Parcel parcel;

	std::unique_ptr<PseudoAdiabaticScheme> choosePseudoAdiabaticScheme(size_t pseudoadiabaticSchemeID);

	void ascentAlongMoistAdiabat();
	void ascentAlongPseudoAdiabat(size_t pseudoadiabaticSchemeID);
	void ascentAlongDryAdiabat();

	void makeFirstTimeStep();
	void makeTimeStep();

public:
	void runSimulationOn(Parcel& passedParcel, size_t pseudoadiabaticSchemeID);
};

class RungeKuttaDynamics : public DynamicScheme
{
private:
	Parcel parcel;

	std::unique_ptr<PseudoAdiabaticScheme> choosePseudoAdiabaticScheme(size_t pseudoadiabaticSchemeID);

	void ascentAlongMoistAdiabat();
	void ascentAlongPseudoAdiabat(size_t pseudoadiabaticSchemeID);
	void ascentAlongDryAdiabat();

	void makeAdiabaticTimeStep(double lambda, double gamma);
	void makePseudoAdiabaticTimeStep(size_t pseudoadiabaticSchemeID, double wetBulbTemperature);

public:
	void runSimulationOn(Parcel& passedParcel, size_t pseudoadiabaticSchemeID);
};

#endif
