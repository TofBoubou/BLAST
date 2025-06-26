#include "blast/io/config_manager.hpp"
#include <iostream>

int main() {
    std::cout << "ABCDEFG" << std::endl;
    blast::io::ConfigurationManager manager;
    auto config_result = manager.load("config/default.yaml");
    if (!config_result) {
        std::cerr << "Failed to load configuration: " << config_result.error().message() << std::endl;
        return 1;
    }
    std::cout << "Configuration loaded successfully" << std::endl;
    return 0;
}

/* 
$cppFiles = Get-ChildItem -Recurse -Filter *.cpp -Path src | ForEach-Object { $_.FullName }
>>
>> g++ -std=c++23 `
>>     -DYAML_CPP_STATIC_DEFINE `
>>     -Iinclude `
>>     -Ilibs/yaml-cpp/include `
>>     -Ilibs/eigen `
>>     -Llibs/yaml-cpp/build `
>>     $cppFiles `
>>     -lyaml-cpp `
>>     -o blast_main.exe
 */
