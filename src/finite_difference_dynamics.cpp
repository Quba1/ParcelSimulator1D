#include "thermodynamic_calc.h"
#include "environment.h"
#include "parcel.h"
#include "dynamic_scheme.h"
#include "pseudoadiabatic_scheme.h"

Parcel FiniteDifferenceDynamics::runSimulationOn(Parcel& passedParcel)
{
	parcel = passedParcel;

	startFromInitialConditions();

	while (isParcelWithinBounds())
	{
		ascentAlongMoistAdiabat();
		ascentAlongPseudoAdiabat();
	}

	return parcel;
}

void FiniteDifferenceDynamics::ascentAlongMoistAdiabat()
{
	//calculate ascent constants
	double gamma = calcGamma(parcel.mixingRatio[parcel.currentTimeStep]);
	double lambda = calcLambda(parcel.temperature[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep], gamma);

	//loop through next timesteps
	do
	{
		if (!isParcelWithinBounds())
		{
			return;
		}

		makeTimeStep();

		//update parcel properties
		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		parcel.updateCurrentThermodynamicsAdiabatically(lambda, gamma);
	} while (parcel.mixingRatioSaturated[parcel.currentTimeStep] > parcel.mixingRatio[parcel.currentTimeStep]);

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
	while (parcel.mixingRatio[parcel.currentTimeStep] > parcel.noMoistureTreshold && parcel.velocity[parcel.currentTimeStep] > 0)
	{
		if (!isParcelWithinBounds())
		{
			return;
		}

		makeTimeStep();

		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		double pressureDelta = parcel.pressure[parcel.currentTimeStep] - parcel.pressure[parcel.currentTimeStep - 1];
		parcel.temperature[parcel.currentTimeStep] = pseudoadiabaticScheme->calculateCurrentPseudoadiabaticTemperature(parcel.getSlice(-1), pressureDelta, wetBulbPotentialTemp);
		parcel.updateCurrentThermodynamicsPseudoadiabatically();
	}
}

void FiniteDifferenceDynamics::startFromInitialConditions()
{
	double gamma = calcGamma(parcel.mixingRatio[parcel.currentTimeStep]);
	double lambda = calcLambda(parcel.temperature[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep], gamma);

	if (!isParcelWithinBounds())
	{
		return;
	}

	makeFirstTimeStep();

	//update parcel properties
	parcel.currentTimeStep++;
	parcel.updateCurrentDynamicsAndPressure();
	parcel.updateCurrentThermodynamicsAdiabatically(lambda, gamma);
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

bool FiniteDifferenceDynamics::isParcelWithinBounds()
{
	if (parcel.position[parcel.currentTimeStep] >= Environment::highestPoint)
	{
		return false;
	}
	else if (parcel.position[parcel.currentTimeStep] <= 0.0)
	{
		return false;
	}
	else if (parcel.currentTimeStep >= parcel.ascentSteps - 1)
	{
		return false;
	}
	else
	{
		return true;
	}

}