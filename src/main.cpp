#include "blast/io/config_manager.hpp"
#include "blast/boundary_layer/grid/grid.hpp"
#include <iostream>

int main() {
    blast::io::ConfigurationManager manager;
    auto config_result = manager.load("config/default.yaml");
    if (!config_result) {
        std::cerr << "Failed to load configuration: " << config_result.error().message() << std::endl;
        return 1;
    }
    std::cout << "Configuration loaded successfully" << std::endl;

    auto config = config_result.value();
    auto stagnation_grid_result = blast::boundary_layer::grid::BoundaryLayerGrid::create_stagnation_grid(config.numerical, config.outer_edge);

    if (!stagnation_grid_result) {
        std::cerr << "Failed to create a grid: " << stagnation_grid_result.error().message() << std::endl;
        return 1;
    }
    std::cout << "Grid loaded successfully" << std::endl;

    auto grid = stagnation_grid_result.value();
    auto grid_eta = grid.eta_coordinates();
    std::cout << grid_eta[50] << std::endl;




    return 0;
}










/* 
$cppFiles = Get-ChildItem -Recurse -Filter *.cpp -Path src | ForEach-Object { $_.FullName }

g++ -std=c++23 `
    -DYAML_CPP_STATIC_DEFINE `
    -Iinclude `
    -Ilibs/yaml-cpp/include `
    -Ilibs/eigen `
    -Llibs/yaml-cpp/build `
    $cppFiles `
    -lyaml-cpp `
    -static-libgcc `
    -static-libstdc++ `
    -o blast_main.exe
 */
