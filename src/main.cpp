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

//define structures

struct Environment
{
    std::vector<double> height, pressure, temperature, dewpoint;
};

struct Sector
{
    size_t upper = 1;
    size_t lower = 0;
};

struct Parcel
{
    std::unique_ptr<double[]> posZ;
    std::unique_ptr<double[]> velZ;
    std::unique_ptr<double[]> pres;
    std::unique_ptr<double[]> temp;
    std::unique_ptr<double[]> tmpV;
    std::unique_ptr<double[]> mxrt;
    std::unique_ptr<double[]> smxr;

    Parcel(size_t size)
    {
        posZ = std::make_unique<double[]>(size + 1);
        velZ = std::make_unique<double[]>(size + 1);
        pres = std::make_unique<double[]>(size + 1);
        temp = std::make_unique<double[]>(size + 1);
        tmpV = std::make_unique<double[]>(size + 1);
        mxrt = std::make_unique<double[]>(size + 1);
        smxr = std::make_unique<double[]>(size + 1);
    }
};

//define functions

void updateSector(double location, Sector& sctr, const std::vector<double>& heightField)
{
    //assuming sorted array of ascending values in heightField and location within bounds of heightField

    size_t nearestPoint;

    //check which sector boundary is closest
    double distToUpper = abs(location - heightField[sctr.upper]);
    double distToLower = abs(location - heightField[sctr.lower]);

    if (distToLower == distToUpper)
    {
        //no need for sector correction
        return;
    }
    else if (distToLower < distToUpper)
    {
        if (sctr.lower == 0)
        {
            //location is in the lowest sector
            sctr.upper = 1;
            return;
        }

        //move sector gradually down
        nearestPoint = sctr.lower;
        double dist = distToLower;
        double newDist = abs(location - heightField[nearestPoint - 1]);

        while (newDist < dist)
        {
            dist = newDist;
            nearestPoint--;

            if (nearestPoint > 0)
            {
                newDist = abs(location - heightField[nearestPoint - 1]);
            }
        }
    }
    else if (distToLower > distToUpper)
    {
        if (sctr.upper == (heightField.size() - 1))
        {
            //location is in the highest sector
            sctr.lower = heightField.size() - 2;
            return;
        }

        //move sector gradually up
        nearestPoint = sctr.upper;
        double dist = distToUpper;
        double newDist = abs(location - heightField[nearestPoint + 1]);

        while (newDist < dist)
        {
            dist = newDist;
            nearestPoint++;

            if (nearestPoint < (heightField.size() - 1))
            {
                newDist = abs(location - heightField[nearestPoint + 1]);
            }
        }
    }

    double distToNearest = location - heightField[nearestPoint];

    if (distToNearest >= 0)
    {
        sctr.lower = nearestPoint;
        sctr.upper = nearestPoint + 1;
        return;
    }
    else
    {
        sctr.lower = nearestPoint - 1;
        sctr.upper = nearestPoint;
        return;
    }

}

double getPressureAtLocation(double location, const Sector& sctr, const std::vector<double>& heightField, const std::vector<double>& pressureField)
{
    //input in m; output in Pa

    //do linear interpolation of pressure within the sector
    double b = (pressureField[sctr.upper] - pressureField[sctr.lower]) / (heightField[sctr.upper] - heightField[sctr.lower]);
    double a = pressureField[sctr.lower] - (heightField[sctr.lower] * b);

    //caluclate and return pressure
    return (a + (b * location)) * 100;
}

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
    updateSector(parcel.posZ[i], sector, environment.height);

    //calculate thermodynamic properties for adiabatic ascent
    parcel.pres[i] = getPressureAtLocation(parcel.posZ[i], sector, environment.height, environment.pressure); //pressure of parcel always equalises with atmosphere
    parcel.mxrt[i] = parcel.mxrt[i - 1]; //mixing ratio is conservative during adiabatic ascent
    parcel.temp[i] = calcTemperatureInAdiabat(parcel.pres[i], gamma, lambda);
    parcel.smxr[i] = calcMixingRatio(parcel.temp[i], parcel.pres[i]);
    parcel.tmpV[i] = calcVirtualTemperature(parcel.temp[i], parcel.mxrt[i]);
}

void updateParcelThermodynamicsPseudo(size_t i, Parcel& parcel, Sector& sector, const Environment& environment, const double thetaW)
{
    updateSector(parcel.posZ[i], sector, environment.height);

    parcel.pres[i] = getPressureAtLocation(parcel.posZ[i], sector, environment.height, environment.pressure);
    parcel.temp[i] = calcTemperatureInPseudoC(parcel.pres[i], thetaW);
    parcel.smxr[i] = calcMixingRatio(parcel.temp[i], parcel.pres[i]);
    parcel.mxrt[i] = parcel.smxr[i];
    parcel.tmpV[i] = calcVirtualTemperature(parcel.temp[i], parcel.mxrt[i]);
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

void readConfiguration(std::map<std::string, std::string>& config)
{
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

    configFile.close();
}

void importEnvData(std::string filename, Environment& env)
{
    std::ifstream txtProfile(filename);
    std::string line;

    //dummy getlines for skipping the header lines
    getline(txtProfile, line);
    getline(txtProfile, line);


    //read variables line by line and push into separate vectors
    while (getline(txtProfile, line))
    {
        std::stringstream lineStream(line);
        std::string var;

        getline(lineStream, var, ';');
        env.height.push_back(stod(var));
        getline(lineStream, var, ';');
        env.pressure.push_back(stod(var));
        getline(lineStream, var, ';');
        env.temperature.push_back(stod(var));
        getline(lineStream, var, ';');
        env.dewpoint.push_back(stod(var));
    }

    txtProfile.close();
}

void outputData(const Parcel& parcel, std::string filename, size_t steps)
{
    std::ofstream output(filename);
    output << std::fixed << std::setprecision(5);

    output << "position; velocity; pressure; temperature; virtual_temperature; mixing_ratio; saturation_mixing_ratio;" << "\n";

    for (size_t i = 0; i < steps; i++)
    {
        output << parcel.posZ[i] << "; "
            << parcel.velZ[i] << "; "
            << parcel.pres[i] << "; "
            << parcel.temp[i] << "; "
            << parcel.tmpV[i] << "; "
            << parcel.mxrt[i] << "; "
            << parcel.smxr[i] << ";" << "\n";
    }

    output.close();
}

int main()
{
    //------------------------ read model configuration ---------------------//
    std::map<std::string, std::string> configuration;
    readConfiguration(configuration);

    unsigned char scheme = stoi(configuration["dynamic_scheme"]);
    unsigned char pseudo_scheme = stoi(configuration["pseudoadiabatic_scheme"]);
    std::string outputFilename = "output/" + configuration["output_file"];

    //------------------------ import environmental variables ---------------------//
    
    Environment environment;
    importEnvData("input/" + configuration["profile_file"], environment);

    //------------------------ set parcel and initial conditions ---------------------//

    //read values from file into variables and make calculations
    double const dt = stod(configuration["timestep"]);
    double const period = stod(configuration["period"]);
    size_t itr = 0;

    unsigned int const steps = static_cast<int>(floor((period * 3600) / dt));
    double const dt2 = dt * dt;

    //create parcel
    Parcel parcel(steps);

    //write initial conditions into parcel and convert to SI units
    parcel.posZ[itr] = stod(configuration["init_height"]);
    parcel.velZ[itr] = stod(configuration["init_velocity"]);
    parcel.temp[itr] = stod(configuration["init_temp"]) + 273.15;

    //determine initial sector parcel is in
    Sector sector;
    updateSector(parcel.posZ[0], sector, environment.height);

    //convert initial variables to other intermediate variables (in SI units)
    parcel.pres[itr] = getPressureAtLocation(parcel.posZ[itr], sector, environment.height, environment.pressure);
    parcel.mxrt[itr] = calcMixingRatio((stod(configuration["init_dewpoint"]) + 273.15), parcel.pres[itr]);
    parcel.tmpV[itr] = calcVirtualTemperature(parcel.temp[itr], parcel.mxrt[itr]);
    parcel.smxr[itr] = calcMixingRatio(parcel.temp[itr], parcel.pres[itr]);

    //------------------------ compute moist adiabatic ascent ---------------------//

    //constants for adaibatic ascent
    double gamma = (1005.7 * ((1 + (parcel.mxrt[itr] * (1870.0 / 1005.7))) / (1 + parcel.mxrt[itr]))) / (718.0 * ((1 + (parcel.mxrt[itr] * (1410.0 / 718.0))) / (1 + parcel.mxrt[itr]))); //Bailyn (1994)
    double lambda = pow(parcel.pres[itr], 1 - gamma) * pow(parcel.temp[itr], gamma);

    if (scheme == 1)
    {
        //compute first timestep
        //using initial velocity and forward-time finite difference (1st order)
        itr++;

        parcel.posZ[itr] = parcel.posZ[0] + (parcel.velZ[itr - 1] * dt);
        updateParcelThermodynamicsAdiabatic(itr, parcel, sector, environment, gamma, lambda);

        //loop through adiabatic ascent
        //using centered-time finite difference (2nd order) for second derivative
        while (parcel.smxr[itr] > parcel.mxrt[itr])
        {
            itr++;

            double envVirtualTemperature = getEnvVirtualTemperature(parcel.posZ[itr - 1], sector, environment);
            double bouyancyForce = calcBouyancyForce(parcel.tmpV[itr - 1], envVirtualTemperature);

            parcel.posZ[itr] = (dt2 * bouyancyForce) + (2 * parcel.posZ[itr - 1]) - parcel.posZ[itr - 2];
            parcel.velZ[itr - 1] = calcVelocity(parcel.posZ[itr - 1], parcel.posZ[itr], dt);
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
        while (parcel.smxr[itr] > parcel.mxrt[itr])
        {
            Sector halfSector = sector;
            double halfPres, halfTemp, halfTempV, halfPos;

            double C0 = parcel.velZ[itr];
            double K0 = calcBouyancyForce(parcel.tmpV[itr], getEnvVirtualTemperature(parcel.posZ[itr], sector, environment));

            double C1 = parcel.velZ[itr] + (0.5 * dt * K0);           
            halfPos = parcel.posZ[itr] + (0.5 * dt * C0);
            updateRKStepAdiabatic(halfPos, halfPres, halfTemp, halfTempV, halfSector, environment, parcel.mxrt[itr], gamma, lambda);
            double K1 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            double C2 = parcel.velZ[itr] + (0.5 * dt * K1);
            halfPos = parcel.posZ[itr] + (0.5 * dt * C1);
            updateRKStepAdiabatic(halfPos, halfPres, halfTemp, halfTempV, halfSector, environment, parcel.mxrt[itr], gamma, lambda);
            double K2 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            double C3 = parcel.velZ[itr] +  (dt * K2);
            halfPos = parcel.posZ[itr] + (dt * C2);
            updateRKStepAdiabatic(halfPos, halfPres, halfTemp, halfTempV, halfSector, environment, parcel.mxrt[itr], gamma, lambda);
            double K3 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            parcel.posZ[itr + 1] = parcel.posZ[itr] + ((dt / 6.0) * (C0 + 2 * C1 + 2 * C2 + C3));
            parcel.velZ[itr + 1] = parcel.velZ[itr] + ((dt / 6.0) * (K0 + 2 * K1 + 2 * K2 + K3));

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
    parcel.mxrt[itr] = parcel.smxr[itr];

    //------------------------ compute pseudoadiabatic ascent ---------------------//

    //calculate wet-bulb potential temperature for pseudoadiabatic ascent
    const double thetaW = calcWBPotentialTemperature(parcel.temp[itr], parcel.mxrt[itr], parcel.smxr[itr], parcel.pres[itr]);

    if (scheme == 1)
    {
        //loop through timesteps until point of no moisture
        while (parcel.mxrt[itr] > 0.0001)
        {
            itr++;

            double envVirtualTemperature = getEnvVirtualTemperature(parcel.posZ[itr - 1], sector, environment);
            double bouyancyForce = calcBouyancyForce(parcel.tmpV[itr - 1], envVirtualTemperature);

            parcel.posZ[itr] = (dt2 * bouyancyForce) + (2 * parcel.posZ[itr - 1]) - parcel.posZ[itr - 2];
            parcel.velZ[itr - 1] = calcVelocity(parcel.posZ[itr - 1], parcel.posZ[itr], dt);
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
        while (parcel.mxrt[itr] > 0.0001)
        {
            Sector halfSector = sector;
            double halfPres, halfTemp, halfTempV, halfMixr, halfPos;

            double C0 = parcel.velZ[itr];
            double K0 = calcBouyancyForce(parcel.tmpV[itr], getEnvVirtualTemperature(parcel.posZ[itr], sector, environment));

            double C1 = parcel.velZ[itr] + (0.5 * dt * K0);
            halfPos = parcel.posZ[itr] + (0.5 * dt * C0);
            updateRKStepPseudo(halfPos, halfPres, halfTemp, halfTempV, halfMixr, halfSector, environment, thetaW);
            double K1 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            double C2 = parcel.velZ[itr] + (0.5 * dt * K1);
            halfPos = parcel.posZ[itr] + (0.5 * dt * C1);
            updateRKStepPseudo(halfPos, halfPres, halfTemp, halfTempV, halfMixr, halfSector, environment, thetaW);
            double K2 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            double C3 = parcel.velZ[itr] + (dt * K2);
            halfPos = parcel.posZ[itr] + (dt * C2);
            updateRKStepPseudo(halfPos, halfPres, halfTemp, halfTempV, halfMixr, halfSector, environment, thetaW);
            double K3 = calcBouyancyForce(halfTempV, getEnvVirtualTemperature(halfPos, halfSector, environment));

            parcel.posZ[itr + 1] = parcel.posZ[itr] + ((dt / 6.0) * (C0 + 2 * C1 + 2 * C2 + C3));
            parcel.velZ[itr + 1] = parcel.velZ[itr] + ((dt / 6.0) * (K0 + 2 * K1 + 2 * K2 + K3));

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
    lambda = pow(parcel.pres[itr], 1 - gamma) * pow(parcel.temp[itr], gamma);

    if (scheme == 1)
    {
        //loop through timesteps
        while (parcel.velZ[itr - 1] > 0)
        {
            itr++;

            double envTemperature = getTemperatureAtLocation(parcel.posZ[itr - 1], sector, environment.height, environment.temperature);
            double bouyancyForce = calcBouyancyForce(parcel.tmpV[itr - 1], envTemperature);

            parcel.posZ[itr] = (dt2 * bouyancyForce) + (2 * parcel.posZ[itr - 1]) - parcel.posZ[itr - 2];
            parcel.velZ[itr - 1] = calcVelocity(parcel.posZ[itr - 1], parcel.posZ[itr], dt);
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
