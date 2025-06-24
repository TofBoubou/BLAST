# Ajouter les sources IO
target_sources(blast_lib PRIVATE
    io/config_manager.cpp
    io/yaml_parser.cpp
    io/validation.cpp
    io/config_types.cpp
)

# Ajouter yaml-cpp dependency
find_package(yaml-cpp REQUIRED)
target_link_libraries(blast_lib PRIVATE yaml-cpp)