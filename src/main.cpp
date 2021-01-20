#include "dynamic_scheme.h"
#include "environment.h"
#include "parcel.h"
#include "pseudoadiabatic_scheme.h"
#include "thermodynamic_calc.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

//define functions

void updateRKStepAdiabatic(double& halfPos, double& halfPres, double& halfTemp, double& halfTempV, Sector& halfSector, const Environment& environment, double mixr, const double gamma, const double lambda)
{
    updateSector(halfPos, halfSector, environment.height);
    halfPres = getPressureAtLocation(halfPos, halfSector, environment.height, environment.pressure);
    halfTemp = calcTemperatureInAdiabat(halfPres, gamma, lambda);
    halfTempV = calcVirtualTemperature(halfTemp, mixr);
}

void updateRKStepPseudo(double& halfPos, double& halfPres, double& halfTemp, double& halfTempV, double& halfMixr, Sector& halfSector, const Environment& environment, const double thetaW)
{
    updateSector(halfPos, halfSector, environment.height);
    halfPres = getPressureAtLocation(halfPos, halfSector, environment.height, environment.pressure);
    halfTemp = calcTemperatureInPseudoC(halfPres, thetaW);
    halfMixr = calcMixingRatio(halfTemp, halfPres);
    halfTempV = calcVirtualTemperature(halfTemp, halfMixr);
}

std::map<std::string, std::string> readConfiguration()
{
    std::map<std::string, std::string> config;
    std::ifstream configFile("config/model.conf");
    std::string line;

    //read file into a map
    while (getline(configFile, line))
    {
        if (line[0] != '#' && line[0] != '\n' && line[0] != '\0') //read only lines containing variables
        {
            std::stringstream lineStream(line);
            std::string key, value;

            getline(lineStream, key, '=');
            getline(lineStream, value, '=');

            config.insert({ key, value });
        }
    }

    config["profile_filename"] = "input/" + config["profile_file"];
    config["output_filename"] = "output/" + config["output_filename"];

    configFile.close();

    return config;
}

void outputData(const Parcel& parcel, std::string filename, size_t steps)
{
    std::ofstream output(filename);
    output << std::fixed << std::setprecision(5);

    output << "position; velocity; pressure; temperature; virtual_temperature; mixing_ratio; saturation_mixing_ratio;" << "\n";

    for (size_t i = 0; i < steps; i++)
    {
        output << parcel.position[i] << "; "
            << parcel.velocity[i] << "; "
            << parcel.pressure[i] << "; "
            << parcel.temperature[i] << "; "
            << parcel.temperatureVirtual[i] << "; "
            << parcel.mixingRatio[i] << "; "
            << parcel.mixingRatioSaturated[i] << ";" << "\n";
    }

    output.close();
}

int main()
{
    //read model configuration
    const std::map<std::string, std::string> configuration = readConfiguration();  

    //create environment from given profile file
    Environment environment(configuration.at("profile_filename"));

    //create parcel
    Parcel parcel(configuration);

    //create instances of schemes
    size_t dynamicSchemeID = stoi(configuration.at("dynamic_scheme"));
    size_t pseudoadiabaticSchemeID = stoi(configuration.at("pseudoadiabatic_scheme"));

    std::unique_ptr<DynamicScheme> dynamicScheme;
    std::unique_ptr<PseudoAdiabaticScheme> pseudoadiabaticScheme;

    switch (dynamicSchemeID)
    {
    case 1:
        dynamicScheme = std::make_unique<FiniteDifferenceDynamics>(parcel);

    case 2:
        dynamicScheme = std::make_unique<RungeKuttaDynamics>(parcel);

    default:
        printf("Incorect value of dynamic_scheme in model.conf\n");
        return -1;
    }

    //compute moist adiabatic ascent
 
    if (scheme == 2) //https://math.stackexchange.com/a/2023862
    {
        //loop through timesteps
        //running two parallel Runge-Kutta schemes
        while (parcel.mixingRatioSaturated[itr] > parcel.mixingRatio[itr])
        {
            Sector halfSector = sector;
            double halfPres, halfTemp, halfTempV, halfPos;

            double C0 = parcel.velocity[itr];
            double K0 = calcBouyancyForce(parcel.temperatureVirtual[itr], getEnvVirtualTemperature(parcel.position[itr], sector, environment));

            double C1 = parcel.velocity[itr] + (0.5 * dt * K0);
            halfPos = parcel.position[itr] + (0.5 * dt * C0);
            updateRKStepAdiabatic(halfPos, halfPres, halfTemp, halfTempV, halfSector, environment, parcel.mixingRatio[itr], gamma, lambda);
            double K1 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            double C2 = parcel.velocity[itr] + (0.5 * dt * K1);
            halfPos = parcel.position[itr] + (0.5 * dt * C1);
            updateRKStepAdiabatic(halfPos, halfPres, halfTemp, halfTempV, halfSector, environment, parcel.mixingRatio[itr], gamma, lambda);
            double K2 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            double C3 = parcel.velocity[itr] + (dt * K2);
            halfPos = parcel.position[itr] + (dt * C2);
            updateRKStepAdiabatic(halfPos, halfPres, halfTemp, halfTempV, halfSector, environment, parcel.mixingRatio[itr], gamma, lambda);
            double K3 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            parcel.position[itr + 1] = parcel.position[itr] + ((dt / 6.0) * (C0 + 2 * C1 + 2 * C2 + C3));
            parcel.velocity[itr + 1] = parcel.velocity[itr] + ((dt / 6.0) * (K0 + 2 * K1 + 2 * K2 + K3));

            updateParcelThermodynamicsAdiabatic((itr + 1), parcel, sector, environment, gamma, lambda);

        }
    }

    //------------------------ compute pseudoadiabatic ascent ---------------------//
    if (scheme == 2)
    {
        //loop through timesteps until point of no moisture
        while (parcel.mixingRatio[itr] > 0.0001)
        {
            Sector halfSector = sector;
            double halfPres, halfTemp, halfTempV, halfMixr, halfPos;

            double C0 = parcel.velocity[itr];
            double K0 = calcBouyancyForce(parcel.temperatureVirtual[itr], getEnvVirtualTemperature(parcel.position[itr], sector, environment));

            double C1 = parcel.velocity[itr] + (0.5 * dt * K0);
            halfPos = parcel.position[itr] + (0.5 * dt * C0);
            updateRKStepPseudo(halfPos, halfPres, halfTemp, halfTempV, halfMixr, halfSector, environment, thetaW);
            double K1 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            double C2 = parcel.velocity[itr] + (0.5 * dt * K1);
            halfPos = parcel.position[itr] + (0.5 * dt * C1);
            updateRKStepPseudo(halfPos, halfPres, halfTemp, halfTempV, halfMixr, halfSector, environment, thetaW);
            double K2 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            double C3 = parcel.velocity[itr] + (dt * K2);
            halfPos = parcel.position[itr] + (dt * C2);
            updateRKStepPseudo(halfPos, halfPres, halfTemp, halfTempV, halfMixr, halfSector, environment, thetaW);
            double K3 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            parcel.position[itr + 1] = parcel.position[itr] + ((dt / 6.0) * (C0 + 2 * C1 + 2 * C2 + C3));
            parcel.velocity[itr + 1] = parcel.velocity[itr] + ((dt / 6.0) * (K0 + 2 * K1 + 2 * K2 + K3));

            updateParcelThermodynamicsPseudo((itr + 1), parcel, sector, environment, thetaW);

            itr++;

            if (itr == steps)
            {
                outputData(parcel, outputFilename, itr);
                return 0;
            }
        }
    }

    //------------------------ compute dry adiabatic ascent ---------------------//

    outputData(parcel, outputFilename, itr);
    return 0;
}
