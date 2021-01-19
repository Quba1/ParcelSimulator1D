#include "dynamic_scheme.h"
#include "environment.h"
#include "parcel.h"

DynamicPair::DynamicPair()
{
	location = Environment::Location();
	velocity = 0;
}

FiniteDifferenceDynamics::FiniteDifferenceDynamics(double timeDelta) :
	timeDelta(timeDelta)
{
	timeDeltaSquared = timeDelta * timeDelta;
}

RungeKuttaDynamics::RungeKuttaDynamics(double timeDelta) :
	timeDelta(timeDelta)
{
	timeDeltaSquared = timeDelta * timeDelta;
}

DynamicPair FiniteDifferenceDynamics::computeFirstTimeStep(double bouyancyForce, DynamicPair currentDynamicPair)
{
	DynamicPair nextDynamicPair;

	nextDynamicPair.location.position = currentDynamicPair.location.position + (currentDynamicPair.velocity * timeDelta);
	nextDynamicPair.velocity = (nextDynamicPair.location.position - currentDynamicPair.location.position) / timeDelta;
	nextDynamicPair.location.updateSector();

	return nextDynamicPair;
}

DynamicPair FiniteDifferenceDynamics::computeTimeStep(double bouyancyForce, DynamicPair currentDynamicPair)
{

}

DynamicPair RungeKuttaDynamics::computeFirstTimeStep(double bouyancyForce, DynamicPair currentDynamicPair)
{
	//redundant function to increase code safety and allow for universal interface of Dynamic Scheme
	return computeTimeStep(bouyancyForce, currentDynamicPair);
}

DynamicPair RungeKuttaDynamics::computeTimeStep(double bouyancyForce, DynamicPair currentDynamicPair)
{

}