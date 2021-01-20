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

std::map<std::string, std::string> readConfigurationFromFile(std::string filename)
{
    std::map<std::string, std::string> config;
    std::ifstream configFile("config/" + filename);
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

    if(config.find("profile_filename") != config.end())
    {
        config["profile_filename"] = "input/" + config["profile_file"];
    }
        

    if (config.find("output_filename") != config.end())
    {
        config["output_filename"] = "output/" + config["output_filename"];
    }
       
    configFile.close();

    return config;
}

void outputDataFrom(const Parcel& parcel)
{
    std::ofstream output(parcel.outputFileName);
    output << std::fixed << std::setprecision(5);

    output << "position; velocity; pressure; temperature; virtual_temperature; mixing_ratio; saturation_mixing_ratio;" << "\n";

    for (size_t i = 0; i < parcel.ascentSteps; i++)
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
    //read model and parcel configuration
    const std::map<std::string, std::string> modelConfiguration = readConfigurationFromFile("model.conf");
    const std::map<std::string, std::string> parcelConfiguration = readConfigurationFromFile("parcel.conf");

    //create environment from given profile file
    Environment environment(modelConfiguration.at("profile_filename"));

    //create parcel
    Parcel parcel(parcelConfiguration);

    //create instances of schemes
    std::unique_ptr<DynamicScheme> dynamicScheme;

    size_t dynamicSchemeID = stoi(modelConfiguration.at("dynamic_scheme"));
    size_t pseudoAdiabaticSchemeID = stoi(parcelConfiguration.at("pseudoadiabatic_scheme"));

    switch (dynamicSchemeID)
    {
    case 1:
        dynamicScheme = std::make_unique<FiniteDifferenceDynamics>();

    case 2:
        dynamicScheme = std::make_unique<RungeKuttaDynamics>();

    default:
        printf("Incorect value of dynamic_scheme in model.conf\n");
        return -1;
    }

    dynamicScheme->runSimulationOn(parcel, pseudoAdiabaticSchemeID);

    outputDataFrom(parcel);
    return 0;
}
