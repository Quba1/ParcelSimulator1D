#ifndef DYNAMICS_H
#define DYNAMICS_H

class DynamicScheme
{
private:

public:
	virtual double calculateStep(double bouyancy, double location) = 0;
};

class FiniteDifferenceDynamics : DynamicScheme
{
private:

public:
	FiniteDifferenceDynamics() {};
};

class RungeKuttaDynamics : DynamicScheme
{
private:

public:
	RungeKuttaDynamics() {};
};

#endif
