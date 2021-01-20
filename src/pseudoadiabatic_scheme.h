#ifndef PSEUDOADIABAT_H
#define PSEUDOADIABAT_H

#include "parcel.h"

class PseudoAdiabaticScheme
{
public:
	PseudoAdiabaticScheme() {};
	virtual double calculateCurrentPseudoadiabaticTemperature(Parcel::Slice parcelSlice, double WetBulbTheta) = 0;

private:

};

class FiniteDifferencePseudoadiabat : public PseudoAdiabaticScheme
{
public:
	FiniteDifferencePseudoadiabat() {};
	double calculateCurrentPseudoadiabaticTemperature(Parcel::Slice parcelSlice, double WetBulbTheta);

private:

};

class RungeKuttaPseudoadiabat : public PseudoAdiabaticScheme
{
public:
	RungeKuttaPseudoadiabat() {};
	double calculateCurrentPseudoadiabaticTemperature(Parcel::Slice parcelSlice, double WetBulbTheta);

private:

};

class NumericalPseudoadiabat : public PseudoAdiabaticScheme
{
public:
	NumericalPseudoadiabat() {};
	double calculateCurrentPseudoadiabaticTemperature(Parcel::Slice parcelSlice, double WetBulbTheta);

private:

};
#endif
