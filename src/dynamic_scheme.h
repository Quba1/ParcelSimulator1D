#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "thermodynamic_calc.h"
#include "environment.h"
#include "parcel.h"
#include "pseudoadiabatic_scheme.h"
#include <memory>

class DynamicScheme
{
public:
	virtual Parcel runSimulationOn(Parcel& passedParcel) = 0;
};

class FiniteDifferenceDynamics : public DynamicScheme
{
private:
	Parcel parcel;

	std::unique_ptr<PseudoAdiabaticScheme> choosePseudoAdiabaticScheme();

	void ascentAlongMoistAdiabat();
	void ascentAlongPseudoAdiabat();
	void ascentAlongDryAdiabat();

	void makeFirstTimeStep();
	void makeTimeStep();

	bool isParcelWithinBounds();

public:
	Parcel runSimulationOn(Parcel& passedParcel);
};

class RungeKuttaDynamics : public DynamicScheme
{
private:
	Parcel parcel;

	std::unique_ptr<PseudoAdiabaticScheme> choosePseudoAdiabaticScheme();

	void ascentAlongMoistAdiabat();
	void ascentAlongPseudoAdiabat();
	void ascentAlongDryAdiabat();

	void makeAdiabaticTimeStep(double lambda, double gamma);
	void makePseudoAdiabaticTimeStep(double wetBulbTemperature);

	bool isParcelWithinBounds();

public:
	Parcel runSimulationOn(Parcel& passedParcel);
};

#endif
