#include "pseudoadiabatic_scheme.h"
#include "thermodynamic_calc.h"

double FiniteDifferencePseudoadiabat::calculateStep(double force)
{
	return force;
}

double RungeKuttaPseudoadiabat::calculateStep(double force)
{
	return force;
}

double NumericalPseudoadiabat::calculateStep(double force)
{
	return force;
}