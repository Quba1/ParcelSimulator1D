#include "dynamic_scheme.h"
#include "environment.h"
#include "parcel.h"

double FiniteDifferenceDynamics::calculateStep(double bouynacy, double location)
{
	return bouynacy * location;
}

double RungeKuttaDynamics::calculateStep(double bouynacy, double location)
{
	return bouynacy * location;
}