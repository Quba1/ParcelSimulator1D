#include "thermodynamic_calc.h"
#include "environment.h"
#include "parcel.h"
#include "dynamic_scheme.h"
#include "pseudoadiabatic_scheme.h"

Parcel FiniteDifferenceDynamics::runSimulationOn(Parcel& passedParcel)
{
	parcel = passedParcel;

	ascentAlongMoistAdiabat();
	ascentAlongPseudoAdiabat();
	ascentAlongDryAdiabat();

	return parcel;
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

void FiniteDifferenceDynamics::ascentAlongPseudoAdiabat()
{
	//create pseudodynamic scheme
	std::unique_ptr<PseudoAdiabaticScheme> pseudoadiabaticScheme = choosePseudoAdiabaticScheme();

	//calculate wet-bulb potential temperature for pseudoadiabatic ascent
	double wetBulbPotentialTemp = calcWBPotentialTemperature(parcel.temperature[parcel.currentTimeStep], parcel.mixingRatio[parcel.currentTimeStep], parcel.mixingRatioSaturated[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep]);

	//loop through timesteps until point of no moisture
	while (parcel.mixingRatio[parcel.currentTimeStep] > parcel.noMoistureTreshold && parcel.currentTimeStep < parcel.ascentSteps)
	{
		makeTimeStep();

		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		double pressureDelta = parcel.pressure[parcel.currentTimeStep] - parcel.pressure[parcel.currentTimeStep - 1];
		parcel.temperature[parcel.currentTimeStep] = pseudoadiabaticScheme->calculateCurrentPseudoadiabaticTemperature(parcel.getSlice(-1), pressureDelta, wetBulbPotentialTemp);
		parcel.updateCurrentThermodynamicsPseudoadiabatically();
	}
}

void FiniteDifferenceDynamics::ascentAlongDryAdiabat()
{
	double gamma = C_P / C_V; //simiplified gamma for dry air
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

	parcel.position[parcel.currentTimeStep + 1] = (parcel.timeDeltaSquared * bouyancyForce) + (2.0 * parcel.position[parcel.currentTimeStep]) - parcel.position[parcel.currentTimeStep - 1];
	parcel.velocity[parcel.currentTimeStep + 1] = (parcel.position[parcel.currentTimeStep + 1] - parcel.position[parcel.currentTimeStep]) / parcel.timeDelta;
}

std::unique_ptr<PseudoAdiabaticScheme> FiniteDifferenceDynamics::choosePseudoAdiabaticScheme()
{
	size_t pseudoadiabaticSchemeID = std::stoi(parcel.parcelConfiguration.at("pseudoadiabatic_scheme"));

	if (pseudoadiabaticSchemeID == 1)
	{
		return std::make_unique<FiniteDifferencePseudoadiabat>();
	}
	else if (pseudoadiabaticSchemeID == 2)
	{
		return std::make_unique<RungeKuttaPseudoadiabat>();
	}
	else if (pseudoadiabaticSchemeID == 3)
	{
		return std::make_unique<NumericalPseudoadiabat>();
	}
	else
	{
		return nullptr;
	}
}
