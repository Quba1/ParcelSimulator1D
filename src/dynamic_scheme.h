#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "environment.h"

struct DynamicPair
{
	double velocity;
	Environment::Location location;

	DynamicPair();
};

class DynamicScheme
{
public:
	virtual DynamicPair computeFirstTimeStep(double bouyancyForce, DynamicPair currentDynamicPair) = 0;
	virtual DynamicPair computeTimeStep(double bouyancyForce, DynamicPair currentDynamicPair) = 0;
};

class FiniteDifferenceDynamics : public DynamicScheme
{
private:
	double timeDelta, timeDeltaSquared;

public:
	FiniteDifferenceDynamics(double timeDelta);
	DynamicPair computeFirstTimeStep(double bouyancyForce, DynamicPair currentDynamicPair);
	DynamicPair computeTimeStep(double bouyancyForce, DynamicPair currentDynamicPair);
};

class RungeKuttaDynamics : public DynamicScheme
{
private:
	double timeDelta, timeDeltaSquared;

public:
	RungeKuttaDynamics(double timeDelta);
	DynamicPair computeFirstTimeStep(double bouyancyForce, DynamicPair currentDynamicPair);
	DynamicPair computeTimeStep(double bouyancyForce, DynamicPair currentDynamicPair);
};

#endif
