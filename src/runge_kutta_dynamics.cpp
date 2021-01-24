#include "thermodynamic_calc.h"
#include "environment.h"
#include "parcel.h"
#include "dynamic_scheme.h"
#include "pseudoadiabatic_scheme.h"

Parcel RungeKuttaDynamics::runSimulationOn(Parcel& passedParcel)
{
	parcel = passedParcel;

	ascentAlongMoistAdiabat();
	ascentAlongPseudoAdiabat();
	ascentAlongDryAdiabat();

	return parcel;
}

void RungeKuttaDynamics::ascentAlongMoistAdiabat()
{
	//calculate ascent constants
	double gamma = calcGamma(parcel.mixingRatio[parcel.currentTimeStep]);
	double lambda = calcLambda(parcel.temperature[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep], gamma);

	//loop through next timesteps
	while (parcel.mixingRatioSaturated[parcel.currentTimeStep] > parcel.mixingRatio[parcel.currentTimeStep] && parcel.currentTimeStep < parcel.ascentSteps)
	{
		if (!isParcelWithinBounds())
		{
			return;
		}

		makeAdiabaticTimeStep(lambda, gamma);

		//update parcel properties
		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		parcel.updateCurrentThermodynamicsAdiabatically(lambda, gamma);
	}

	//equalise mixing ratio and saturation mixing ratio at the end of adiabatic ascent
	parcel.mixingRatio[parcel.currentTimeStep] = parcel.mixingRatioSaturated[parcel.currentTimeStep];
}

void RungeKuttaDynamics::ascentAlongPseudoAdiabat()
{
	//create pseudodynamic scheme
	std::unique_ptr<PseudoAdiabaticScheme> pseudoadiabaticScheme = choosePseudoAdiabaticScheme();

	//calculate wet-bulb potential temperature for pseudoadiabatic ascent
	double wetBulbPotentialTemp = calcWBPotentialTemperature(parcel.temperature[parcel.currentTimeStep], parcel.mixingRatio[parcel.currentTimeStep], parcel.mixingRatioSaturated[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep]);

	//loop through timesteps until point of no moisture
	while (parcel.mixingRatio[parcel.currentTimeStep] > parcel.noMoistureTreshold && parcel.currentTimeStep < parcel.ascentSteps)
	{
		if (!isParcelWithinBounds())
		{
			return;
		}

		makePseudoAdiabaticTimeStep(wetBulbPotentialTemp);

		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		double pressureDelta = parcel.pressure[parcel.currentTimeStep] - parcel.pressure[parcel.currentTimeStep - 1];
		parcel.temperature[parcel.currentTimeStep] = pseudoadiabaticScheme->calculateCurrentPseudoadiabaticTemperature(parcel.getSlice(-1), pressureDelta, wetBulbPotentialTemp);
		parcel.updateCurrentThermodynamicsPseudoadiabatically();
	}
}

void RungeKuttaDynamics::ascentAlongDryAdiabat()
{
	//calculate ascent constants
	double gamma = C_P / C_V; //simiplified gamma for dry air
	double lambda = calcLambda(parcel.temperature[parcel.currentTimeStep], parcel.pressure[parcel.currentTimeStep], gamma);

	//loop through next timesteps
	while (parcel.velocity[parcel.currentTimeStep] > 0 && parcel.currentTimeStep < parcel.ascentSteps)
	{
		if (!isParcelWithinBounds())
		{
			return;
		}

		makeAdiabaticTimeStep(lambda, gamma);

		//update parcel properties
		parcel.currentTimeStep++;
		parcel.updateCurrentDynamicsAndPressure();
		parcel.updateCurrentThermodynamicsAdiabatically(lambda, gamma);
	}
}

void RungeKuttaDynamics::makeAdiabaticTimeStep(double lambda, double gamma)
{
	//algorithm source: https://math.stackexchange.com/a/2023862

	double stepTemperature, stepTemperatureVirtual, stepPressure;
	Environment::Location stepLocation = parcel.currentLocation;

	double C0 = parcel.velocity[parcel.currentTimeStep];
	double K0 = calcBouyancyForce(parcel.temperatureVirtual[parcel.currentTimeStep], Environment::getVirtualTemperatureAtLocation(parcel.currentLocation));

	double C1 = C0 + (0.5 * parcel.timeDelta * K0);
	stepLocation.position = parcel.currentLocation.position + (0.5 * parcel.timeDelta * C0);
	stepLocation.updateSector();
	stepPressure = Environment::getPressureAtLocation(stepLocation);
	stepTemperature = calcTemperatureInAdiabat(stepPressure, gamma, lambda);
	stepTemperatureVirtual = calcVirtualTemperature(stepTemperature, parcel.mixingRatio[parcel.currentTimeStep]);
	double K1 = calcBouyancyForce(stepTemperatureVirtual, Environment::getVirtualTemperatureAtLocation(stepLocation));

	double C2 = C0 + (0.5 * parcel.timeDelta * K1);
	stepLocation.position = parcel.currentLocation.position + (0.5 * parcel.timeDelta * C1);
	stepLocation.updateSector();
	stepPressure = Environment::getPressureAtLocation(stepLocation);
	stepTemperature = calcTemperatureInAdiabat(stepPressure, gamma, lambda);
	stepTemperatureVirtual = calcVirtualTemperature(stepTemperature, parcel.mixingRatio[parcel.currentTimeStep]);
	double K2 = calcBouyancyForce(stepTemperatureVirtual, Environment::getVirtualTemperatureAtLocation(stepLocation));

	double C3 = C0 + (0.5 * parcel.timeDelta * K2);
	stepLocation.position = parcel.currentLocation.position + (0.5 * parcel.timeDelta * C2);
	stepLocation.updateSector();
	stepPressure = Environment::getPressureAtLocation(stepLocation);
	stepTemperature = calcTemperatureInAdiabat(stepPressure, gamma, lambda);
	stepTemperatureVirtual = calcVirtualTemperature(stepTemperature, parcel.mixingRatio[parcel.currentTimeStep]);
	double K3 = calcBouyancyForce(stepTemperatureVirtual, Environment::getVirtualTemperatureAtLocation(stepLocation));

	parcel.position[parcel.currentTimeStep + 1] = parcel.position[parcel.currentTimeStep] + ((parcel.timeDelta / 6.0) * (C0 + 2 * C1 + 2 * C2 + C3));
	parcel.velocity[parcel.currentTimeStep + 1] = parcel.velocity[parcel.currentTimeStep] + ((parcel.timeDelta / 6.0) * (K0 + 2 * K1 + 2 * K2 + K3));

}

void RungeKuttaDynamics::makePseudoAdiabaticTimeStep(double wetBulbTemperature)
{
	std::unique_ptr<PseudoAdiabaticScheme> pseudoadiabaticScheme = choosePseudoAdiabaticScheme();

	double stepTemperature, stepTemperatureVirtual, stepPressure, deltaPressure, stepMixingRatio;
	Environment::Location stepLocation = parcel.currentLocation;
	Parcel::Slice stepSlice = parcel.getSlice(0);

	double C0 = parcel.velocity[parcel.currentTimeStep];
	double K0 = calcBouyancyForce(parcel.temperatureVirtual[parcel.currentTimeStep], Environment::getVirtualTemperatureAtLocation(parcel.currentLocation));

	double C1 = C0 + (0.5 * parcel.timeDelta * K0);
	stepLocation.position = parcel.currentLocation.position + (0.5 * parcel.timeDelta * C0);
	stepLocation.updateSector();
	stepPressure = Environment::getPressureAtLocation(stepLocation);
	deltaPressure = stepPressure - stepSlice.pressure;
	stepTemperature = pseudoadiabaticScheme->calculateCurrentPseudoadiabaticTemperature(stepSlice, deltaPressure, wetBulbTemperature);
	stepMixingRatio = calcMixingRatio(stepTemperature, stepPressure);
	stepTemperatureVirtual = calcVirtualTemperature(stepTemperature, stepMixingRatio);
	double K1 = calcBouyancyForce(stepTemperatureVirtual, Environment::getVirtualTemperatureAtLocation(stepLocation));

	double C2 = C0 + (0.5 * parcel.timeDelta * K1);
	stepLocation.position = parcel.currentLocation.position + (0.5 * parcel.timeDelta * C1);
	stepLocation.updateSector();
	stepPressure = Environment::getPressureAtLocation(stepLocation);
	deltaPressure = stepPressure - stepSlice.pressure;
	stepTemperature = pseudoadiabaticScheme->calculateCurrentPseudoadiabaticTemperature(stepSlice, deltaPressure, wetBulbTemperature);
	stepMixingRatio = calcMixingRatio(stepTemperature, stepPressure);
	stepTemperatureVirtual = calcVirtualTemperature(stepTemperature, stepMixingRatio);
	double K2 = calcBouyancyForce(stepTemperatureVirtual, Environment::getVirtualTemperatureAtLocation(stepLocation));

	double C3 = C0 + (0.5 * parcel.timeDelta * K2);
	stepLocation.position = parcel.currentLocation.position + (0.5 * parcel.timeDelta * C2);
	stepLocation.updateSector();
	stepPressure = Environment::getPressureAtLocation(stepLocation);
	deltaPressure = stepPressure - stepSlice.pressure;
	stepTemperature = pseudoadiabaticScheme->calculateCurrentPseudoadiabaticTemperature(stepSlice, deltaPressure, wetBulbTemperature);
	stepMixingRatio = calcMixingRatio(stepTemperature, stepPressure);
	stepTemperatureVirtual = calcVirtualTemperature(stepTemperature, stepMixingRatio);
	double K3 = calcBouyancyForce(stepTemperatureVirtual, Environment::getVirtualTemperatureAtLocation(stepLocation));

	parcel.position[parcel.currentTimeStep + 1] = parcel.position[parcel.currentTimeStep] + ((parcel.timeDelta / 6.0) * (C0 + 2.0 * C1 + 2.0 * C2 + C3));
	parcel.velocity[parcel.currentTimeStep + 1] = parcel.velocity[parcel.currentTimeStep] + ((parcel.timeDelta / 6.0) * (K0 + 2.0 * K1 + 2.0 * K2 + K3));
}

std::unique_ptr<PseudoAdiabaticScheme> RungeKuttaDynamics::choosePseudoAdiabaticScheme()
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

bool RungeKuttaDynamics::isParcelWithinBounds()
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