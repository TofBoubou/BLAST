#include <yaml-cpp/yaml.h>
#include <iostream>
#include <string>

int main() {
    try {
        std::cout << "Test de lecture du fichier YAML default.yaml...\n";
        
        // Charger le fichier YAML
        YAML::Node config = YAML::LoadFile("config/default.yaml");
        std::cout << "Fichier chargé avec succès!\n";
        
        // Tester la lecture de quelques valeurs
        std::cout << "\nValeurs extraites:\n";
        
        if (config["simulation"]) {
            auto sim = config["simulation"];
            std::cout << "- Body type: " << sim["body_type"].as<std::string>() << "\n";
            std::cout << "- Only stagnation point: " << sim["only_stagnation_point"].as<bool>() << "\n";
            std::cout << "- Diffusion type: " << sim["diffusion_type"].as<std::string>() << "\n";
        }
        
        if (config["numerical"]) {
            auto num = config["numerical"];
            std::cout << "- n_eta: " << num["n_eta"].as<int>() << "\n";
            std::cout << "- eta_max: " << num["eta_max"].as<double>() << "\n";
            std::cout << "- Max iterations: " << num["max_iterations"].as<int>() << "\n";
        }
        
        if (config["mixture"]) {
            auto mix = config["mixture"];
            std::cout << "- Mixture name: " << mix["name"].as<std::string>() << "\n";
            std::cout << "- Database: " << mix["thermodynamic_database"].as<std::string>() << "\n";
        }
        
        if (config["output"]) {
            auto out = config["output"];
            std::cout << "- Output directory: " << out["output_directory"].as<std::string>() << "\n";
            
            if (out["x_stations"]) {
                std::cout << "- X stations: ";
                for (const auto& x : out["x_stations"]) {
                    std::cout << x.as<double>() << " ";
                }
                std::cout << "\n";
            }
        }
        
        std::cout << "\nTest réussi! La lecture YAML fonctionne correctement.\n";
        
    } catch (const YAML::Exception& e) {
        std::cerr << "Erreur YAML: " << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Erreur: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}