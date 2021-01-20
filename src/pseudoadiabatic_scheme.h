#ifndef PSEUDOADIABAT_H
#define PSEUDOADIABAT_H

#include "parcel.h"

class PseudoAdiabaticScheme
{
public:
	PseudoAdiabaticScheme() {};
	virtual void calculateCurrentPseudoAdiabaticTemperature(Parcel& parcel) = 0;

private:

};

class FiniteDifferencePseudoadiabat : public PseudoAdiabaticScheme
{
public:
	FiniteDifferencePseudoadiabat() {};
	void calculateCurrentPseudoadiabaticTemperature(Parcel& parcel);

private:

};

class RungeKuttaPseudoadiabat : public PseudoAdiabaticScheme
{
public:
	RungeKuttaPseudoadiabat() {};
	void calculateCurrentPseudoadiabaticTemperature(Parcel& parcel);

private:

};

class NumericalPseudoadiabat : public PseudoAdiabaticScheme
{
public:
	NumericalPseudoadiabat() {};
	void calculateCurrentPseudoadiabaticTemperature(Parcel& parcel);

private:

};
#endif
