#include "blast/io/config_manager.hpp"
#include <iostream>

int main() {
    blast::io::ConfigurationManager manager;
    auto config_result = manager.load("config/default.yaml");
    if (!config_result) {
        std::cerr << "Failed to load configuration: " << config_result.error().message() << std::endl;
        return 1;
    }
    std::cout << "Configuration loaded successfully" << std::endl;
    return 0;
}
