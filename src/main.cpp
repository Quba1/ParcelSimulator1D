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



double getTemperatureAtLocation(double location, const Sector& sctr, const std::vector<double>& heightField, const std::vector<double>& temperatureField)
{
    //input in m; output in K

    //do linear interpolation of pressure within the sector
    double b = (temperatureField[sctr.upper] - temperatureField[sctr.lower]) / (heightField[sctr.upper] - heightField[sctr.lower]);
    double a = temperatureField[sctr.lower] - (heightField[sctr.lower] * b);

    //caluclate and return pressure
    return (a + (b * location)) + 273.15;
}

double getEnvVirtualTemperature(double location, const Sector& sctr, const Environment& env)
{
    double temperature = getTemperatureAtLocation(location, sctr, env.height, env.temperature);
    double dewpoint = getTemperatureAtLocation(location, sctr, env.height, env.dewpoint);
    double pressure = getPressureAtLocation(location, sctr, env.height, env.pressure);
    double mixingRatio = calcMixingRatio(dewpoint, pressure);
    double virtualTemperature = calcVirtualTemperature(temperature, mixingRatio);

    return virtualTemperature;
}

double calcTemperatureInPseudoC(double pressure, double wbTemp)
{
    //input in Pa & K; output in K
        //Bakhshaii & Stull (2013)
        //necesary unit conversions
    pressure = pressure / 1000.0;
    wbTemp -= 273.15;

    if (wbTemp > -30.0 && wbTemp <= 4.0)
    {
        double g1 = -20.3313 - (0.0253 * pressure);
        double g2 = sin(pow(wbTemp + pressure, 0.5)) + (wbTemp / pressure) + pressure - 2.8565;
        double g3 = cos(19.6836 + pow(1.0 + exp(-wbTemp), -(1.0 / 3.0)) + (pressure / 15.0252));
        double g4 = (4.4653 * sin(pow(pressure, 0.5))) - 71.9358;
        double g5 = pow(exp(wbTemp - (2.71828 * cos(pressure / 18.5219))), (1.0 / 6.0));
        double g6 = wbTemp - sin(pow(pressure + wbTemp + atan(wbTemp) + 6.6165, 0.5));

        return g1 + g2 + g3 + g4 + g5 + g6 + 273.15;

    }
    else if (wbTemp > 4.0 && wbTemp <= 21.0)
    {
        double g1 = -9.6285 + cos(log(atan(atan(exp((-9.2121 * wbTemp) / pressure)))));
        double g2 = wbTemp - ((19.9563 / pressure) * atan(wbTemp)) + (pow(wbTemp, 2.0) / (5.47162 * pressure));
        double g3 = sin(log(8.0 * pow(pressure, 3.0))) * log(2.0 * pow(pressure, 1.5));
        double g4 = wbTemp + (((pressure * wbTemp) - pressure + wbTemp) / (pressure - 190.2578));
        double g5 = pressure - ((pressure - 383.0292) / ((15.4014 * pressure) - pow(pressure, 2.0)));
        double g6 = ((1.0 / 3.0) * log(339.0316 - pressure)) + atan(wbTemp - pressure + 95.9839);
        double g7 = -(log(pressure) * ((298.2909 + (16.5109 * pressure)) / (pressure - 2.2183)));

        return g1 + g2 + g3 + g4 + g5 + g6 + g7 + 273.15;
    }
    else if (wbTemp > 21.0 && wbTemp < 45.0)
    {
        double g1 = 0.3919 * pow(wbTemp, 7 / 3) * pow(pressure * (pressure + 15.8148), -1.0);
        double g2 = (19.9724 + (797.7921 / pressure)) * sin(-19.9724 / wbTemp);
        double g3 = pow(log(-3.927765 + wbTemp + pressure) * cos(log(wbTemp + pressure)), 3.0);
        double g4 = pow(exp(pow(wbTemp + pow(1.0 + exp(-pressure), -1.0), 0.5) - 1.5603), 0.5);
        double g5 = pow(pressure + wbTemp, 0.5) * exp(atan((pressure + wbTemp) / 7.9081));
        double g6 = ((pressure / pow(wbTemp, 2)) * std::min(9.6112, pressure - wbTemp)) - 13.73;
        double g7 = sin(pow(sin(std::min(pressure, 17.3170)), 3.0) - pow(pressure, 0.5) + (25.5113 / wbTemp));

        return g1 + g2 + g3 + g4 + g5 + g6 + g7 + 273.15;
    }

    return -999.0;
}

void updateParcelThermodynamicsAdiabatic(size_t i, Parcel& parcel, Sector& sector, const Environment& environment, const double gamma, const double lambda)
{
    updateSector(parcel.position[i], sector, environment.height);

    //calculate thermodynamic properties for adiabatic ascent
    parcel.pressure[i] = getPressureAtLocation(parcel.position[i], sector, environment.height, environment.pressure); //pressure of parcel always equalises with atmosphere
    parcel.mixingRatio[i] = parcel.mixingRatio[i - 1]; //mixing ratio is conservative during adiabatic ascent
    parcel.temperature[i] = calcTemperatureInAdiabat(parcel.pressure[i], gamma, lambda);
    parcel.mixingRatioSaturated[i] = calcMixingRatio(parcel.temperature[i], parcel.pressure[i]);
    parcel.temperatureVirtual[i] = calcVirtualTemperature(parcel.temperature[i], parcel.mixingRatio[i]);
}

void updateParcelThermodynamicsPseudo(size_t i, Parcel& parcel, Sector& sector, const Environment& environment, const double thetaW)
{
    updateSector(parcel.position[i], sector, environment.height);

    parcel.pressure[i] = getPressureAtLocation(parcel.position[i], sector, environment.height, environment.pressure);
    parcel.temperature[i] = calcTemperatureInPseudoC(parcel.pressure[i], thetaW);
    parcel.mixingRatioSaturated[i] = calcMixingRatio(parcel.temperature[i], parcel.pressure[i]);
    parcel.mixingRatio[i] = parcel.mixingRatioSaturated[i];
    parcel.temperatureVirtual[i] = calcVirtualTemperature(parcel.temperature[i], parcel.mixingRatio[i]);
}

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

    //create instances of schemes
    size_t dynamicSchemeID = stoi(configuration.at("dynamic_scheme"));
    size_t pseudoadiabaticSchemeID = stoi(configuration.at("pseudoadiabatic_scheme"));

    DynamicScheme* dynamicScheme = nullptr;
    PseudoAdiabaticScheme* pseudoadiabaticScheme = nullptr;

    switch (dynamicSchemeID)
    {
    case 1:
        dynamicScheme = new FiniteDifferenceDynamics();

    case 2:
        dynamicScheme = new RungeKuttaDynamics();

    default:
        printf("Incorect value of dynamic_scheme in model.conf\n");
        return -1;
    }

    switch (pseudoadiabaticSchemeID)
    {
    case 1:
        pseudoadiabaticScheme = new FiniteDifferencePseudoadiabat();

    case 2:
        pseudoadiabaticScheme = new RungeKuttaPseudoadiabat();

    case 3:
        pseudoadiabaticScheme = new NumericalPseudoadiabat();

    default:
        printf("Incorect value of pseudoadiabatic_scheme in model.conf\n");
        return -1;
    }


    //create parcel
    Parcel parcel(configuration, dynamicScheme, pseudoadiabaticScheme);

    //------------------------ compute moist adiabatic ascent ---------------------//

    //constants for adaibatic ascent
    

    if (scheme == 1)
    {
        //compute first timestep
        //using initial velocity and forward-time finite difference (1st order)
        itr++;

        parcel.position[itr] = parcel.position[0] + (parcel.velocity[itr - 1] * dt);
        updateParcelThermodynamicsAdiabatic(itr, parcel, sector, environment, gamma, lambda);

        //loop through adiabatic ascent
        //using centered-time finite difference (2nd order) for second derivative
        while (parcel.mixingRatioSaturated[itr] > parcel.mixingRatio[itr])
        {
            itr++;

            double envVirtualTemperature = getEnvVirtualTemperature(parcel.position[itr - 1], sector, environment);
            double bouyancyForce = calcBouyancyForce(parcel.temperatureVirtual[itr - 1], envVirtualTemperature);

            parcel.position[itr] = (dt2 * bouyancyForce) + (2 * parcel.position[itr - 1]) - parcel.position[itr - 2];
            parcel.velocity[itr - 1] = calcVelocity(parcel.position[itr - 1], parcel.position[itr], dt);
            updateParcelThermodynamicsAdiabatic(itr, parcel, sector, environment, gamma, lambda);

            if (itr == steps)
            {
                outputData(parcel, outputFilename, itr);
                return 0;
            }
        }
    }
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

            double C3 = parcel.velocity[itr] +  (dt * K2);
            halfPos = parcel.position[itr] + (dt * C2);
            updateRKStepAdiabatic(halfPos, halfPres, halfTemp, halfTempV, halfSector, environment, parcel.mixingRatio[itr], gamma, lambda);
            double K3 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            parcel.position[itr + 1] = parcel.position[itr] + ((dt / 6.0) * (C0 + 2 * C1 + 2 * C2 + C3));
            parcel.velocity[itr + 1] = parcel.velocity[itr] + ((dt / 6.0) * (K0 + 2 * K1 + 2 * K2 + K3));

            updateParcelThermodynamicsAdiabatic((itr + 1), parcel, sector, environment, gamma, lambda);

            itr++;

            if (itr == steps)
            {
                outputData(parcel, outputFilename, itr);
                return 0;
            }
        }
    }

    //equalise mixing ratio and saturation mixing ratio
    parcel.mixingRatio[itr] = parcel.mixingRatioSaturated[itr];

    //------------------------ compute pseudoadiabatic ascent ---------------------//

    //calculate wet-bulb potential temperature for pseudoadiabatic ascent
    const double thetaW = calcWBPotentialTemperature(parcel.temperature[itr], parcel.mixingRatio[itr], parcel.mixingRatioSaturated[itr], parcel.pressure[itr]);

    if (scheme == 1)
    {
        //loop through timesteps until point of no moisture
        while (parcel.mixingRatio[itr] > 0.0001)
        {
            itr++;

            double envVirtualTemperature = getEnvVirtualTemperature(parcel.position[itr - 1], sector, environment);
            double bouyancyForce = calcBouyancyForce(parcel.temperatureVirtual[itr - 1], envVirtualTemperature);

            parcel.position[itr] = (dt2 * bouyancyForce) + (2 * parcel.position[itr - 1]) - parcel.position[itr - 2];
            parcel.velocity[itr - 1] = calcVelocity(parcel.position[itr - 1], parcel.position[itr], dt);
            updateParcelThermodynamicsPseudo(itr, parcel, sector, environment, thetaW);

            if (itr == steps)
            {
                outputData(parcel, outputFilename, itr);
                return 0;
            }
        }
        
    }
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

    //constants for adaibatic ascent
    gamma = 1005.7 / 718.0;
    lambda = pow(parcel.pressure[itr], 1 - gamma) * pow(parcel.temperature[itr], gamma);

    if (scheme == 1)
    {
        //loop through timesteps
        while (parcel.velocity[itr - 1] > 0)
        {
            itr++;

            double envTemperature = getTemperatureAtLocation(parcel.position[itr - 1], sector, environment.height, environment.temperature);
            double bouyancyForce = calcBouyancyForce(parcel.temperatureVirtual[itr - 1], envTemperature);

            parcel.position[itr] = (dt2 * bouyancyForce) + (2 * parcel.position[itr - 1]) - parcel.position[itr - 2];
            parcel.velocity[itr - 1] = calcVelocity(parcel.position[itr - 1], parcel.position[itr], dt);
            updateParcelThermodynamicsAdiabatic(itr, parcel, sector, environment, gamma, lambda);

            if (itr == steps)
            {
                outputData(parcel, outputFilename, itr);
                return 0;
            }
        } 
    }

    outputData(parcel, outputFilename, itr);
    return 0;
}
