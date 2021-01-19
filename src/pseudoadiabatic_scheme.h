#ifndef PSEUDOADIABAT_H
#define PSEUDOADIABAT_H

class PseudoAdiabaticScheme
{
public:
	PseudoAdiabaticScheme() {};
	virtual double calculateStep(double force) = 0;

private:

};

class FiniteDifferencePseudoadiabat : public PseudoAdiabaticScheme
{
public:
	FiniteDifferencePseudoadiabat() {};
	double calculateStep(double force);

private:

};

class RungeKuttaPseudoadiabat : public PseudoAdiabaticScheme
{
public:
	RungeKuttaPseudoadiabat() {};
	double calculateStep(double force);

private:

};

class NumericalPseudoadiabat : public PseudoAdiabaticScheme
{
public:
	NumericalPseudoadiabat() {};
	double calculateStep(double force);

private:

};
#endif
