# Ajouter tests configuration
add_executable(test_config
    unit/test_config_manager.cpp
    unit/test_yaml_parser.cpp
    unit/test_validation.cpp
)

target_link_libraries(test_config 
    PRIVATE 
    blast_lib
    GTest::gtest_main
)

# Copier les fixtures de test
configure_file(fixtures/valid_config.yaml 
    ${CMAKE_BINARY_DIR}/test_fixtures/valid_config.yaml COPYONLY)