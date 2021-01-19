#ifndef DYNAMICS_H
#define DYNAMICS_H

class DynamicScheme
{
private:

public:
	virtual double calculateStep(double bouyancy, double location) = 0;
};

class FiniteDifferenceDynamics : public DynamicScheme
{
private:

public:
	FiniteDifferenceDynamics() {};
	double calculateStep(double bouyancy, double location);
};

class RungeKuttaDynamics : public DynamicScheme
{
private:

public:
	RungeKuttaDynamics() {};
	double calculateStep(double bouyancy, double location);
};

#endif
