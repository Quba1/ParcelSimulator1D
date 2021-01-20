#include "dynamic_scheme.h"
#include "environment.h"
#include "parcel.h"
#include "pseudoadiabatic_scheme.h"
#include "thermodynamic_calc.h"

void FiniteDifferenceDynamics::runSimulationOn(Parcel& passedParcel, size_t pseudoadiabaticSchemeID)
{
	parcel = passedParcel;

	ascentAlongMoistAdiabat();
	ascentAlongPseudoAdiabat(pseudoadiabaticSchemeID);
	ascentAlongDryAdiabat();
}

std::unique_ptr<PseudoAdiabaticScheme> FiniteDifferenceDynamics::choosePseudoAdiabaticScheme(size_t pseudoadiabaticSchemeID)
{
	switch (pseudoadiabaticSchemeID)
	{
	case 1:
		return std::make_unique<FiniteDifferencePseudoadiabat>();

	case 2:
		return std::make_unique<RungeKuttaPseudoadiabat>();

	case 3:
		return std::make_unique<NumericalPseudoadiabat>();
	}
}

void FiniteDifferenceDynamics::ascentAlongMoistAdiabat()
{
	//calculate ascent constants
	double gamma = calcGamma(parcel.mixingRatio[parcel.currentTimeStep]);
	double lambda = calcLambda(parcel.temperature[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep], gamma);

	makeFirstTimeStep();

	//update parcel properties
	parcel.currentTimeStep++;
	parcel.updateCurrentDynamicsAndPressure();
	parcel.updateCurrentThermodynamicsAdiabatically(lambda, gamma);

	//loop through next timesteps
	while (parcel.mixingRatioSaturated[parcel.currentTimeStep] > parcel.mixingRatio[parcel.currentTimeStep] && parcel.currentTimeStep < parcel.ascentSteps)
	{
		makeTimeStep();

		//update parcel properties
		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		parcel.updateCurrentThermodynamicsAdiabatically(lambda, gamma);
	}

	//equalise mixing ratio and saturation mixing ratio at the end of adiabatic ascent
	parcel.mixingRatio[parcel.currentTimeStep] = parcel.mixingRatioSaturated[parcel.currentTimeStep];
}

void FiniteDifferenceDynamics::ascentAlongPseudoAdiabat(size_t pseudoadiabaticSchemeID)
{
	//create pseudodynamic scheme
	std::unique_ptr<PseudoAdiabaticScheme> pseudoadiabaticScheme = choosePseudoAdiabaticScheme(pseudoadiabaticSchemeID);

	//calculate wet-bulb potential temperature for pseudoadiabatic ascent
	double wetBulbPotentialTemp = calcWBPotentialTemperature(parcel.temperature[parcel.currentTimeStep], parcel.mixingRatio[parcel.currentTimeStep], parcel.mixingRatioSaturated[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep]);

	//loop through timesteps until point of no moisture
	while (parcel.mixingRatio[parcel.currentTimeStep] > 0.0001 && parcel.currentTimeStep < parcel.ascentSteps)
	{
		makeTimeStep();

		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		pseudoadiabaticScheme->calculateCurrentPseudoAdiabaticTemperature(parcel, wetBulbPotentialTemp);
		parcel.updateCurrentThermodynamicsPseudoadiabatically();
	}
}

void FiniteDifferenceDynamics::ascentAlongDryAdiabat()
{
	double gamma = 1005.7 / 718.0; //simiplified gamma for dry air
	double lambda = calcLambda(parcel.temperature[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep], gamma);

	while (parcel.velocity[parcel.currentTimeStep] > 0 && parcel.currentTimeStep < parcel.ascentSteps)
	{
		makeTimeStep();

		//update parcel properties
		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		parcel.updateCurrentThermodynamicsAdiabatically(lambda, gamma);
	}
}

void FiniteDifferenceDynamics::makeFirstTimeStep()
{
	parcel.position[1] = parcel.position[0] + (parcel.velocity[0] * parcel.timeDelta);
	parcel.velocity[1] = (parcel.position[1] - parcel.position[0]) / parcel.timeDelta;
}

void FiniteDifferenceDynamics::makeTimeStep()
{
	double bouyancyForce = calcBouyancyForce(parcel.temperatureVirtual[parcel.currentTimeStep], Environment::getVirtualTemperatureAtLocation(parcel.currentLocation));

	parcel.position[parcel.currentTimeStep + 1] = (parcel.timeDeltaSquared * bouyancyForce) + (2 * parcel.position[parcel.currentTimeStep]) - parcel.position[parcel.currentTimeStep - 1];
	parcel.velocity[parcel.currentTimeStep + 1] = (parcel.position[parcel.currentTimeStep + 1] - parcel.position[parcel.currentTimeStep]) / parcel.timeDelta;
}

void RungeKuttaDynamics::runSimulationOn(Parcel& passedParcel, size_t pseudoadiabaticSchemeID)
{
	parcel = passedParcel;

	std::unique_ptr<PseudoAdiabaticScheme> pseudoadiabaticScheme;

	switch (pseudoadiabaticSchemeID)
	{
	case 1:
		pseudoadiabaticScheme = std::make_unique<FiniteDifferencePseudoadiabat>();

	case 2:
		pseudoadiabaticScheme = std::make_unique<RungeKuttaPseudoadiabat>();

	case 3:
		pseudoadiabaticScheme = std::make_unique<NumericalPseudoadiabat>();
	}
}

void RungeKuttaDynamics::computeTimeStep()
{
	//algorithm source: https://math.stackexchange.com/a/2023862

	DynamicPair nextDynamicPair;
	
	double halfTemperatureVirtual;
	Environment::Location halfLocation = currentParcelProperties.dynamics.location;

	double C0 = currentParcelProperties.dynamics.velocity;
	double K0 = currentParcelProperties.bouyancyForce;

	double C1 = C0 + (0.5 * timeDelta * K0);
	halfLocation.position = currentParcelProperties.dynamics.location.position + (0.5 * timeDelta * C0);
	halfLocation.updateSector();
	halfTemperatureVirtual = getVirtualTemperatureinStepAt(halfLocation);


	return nextDynamicPair;
}
