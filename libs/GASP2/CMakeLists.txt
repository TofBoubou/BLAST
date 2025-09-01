cmake_minimum_required(VERSION 3.20)
project(GASP2 LANGUAGES CXX)

# --- Standard ---
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(GASP2_ENABLE_WARNINGS "Enable extra warnings" ON)

find_package(Eigen3 3.3 REQUIRED)

# --- Sorgenti libreria ---
set(GASP2_SOURCES
  src/catalysis/gamma/gamma_dispatcher.cpp
  src/catalysis/gamma/gamma_build.cpp
  src/catalysis/gamma/gamma_core.cpp
  src/catalysis/gamma/bose_flux.cpp
  src/catalysis/gamma/rini_flux.cpp
    src/catalysis/gamma/gamma_model.cpp
    src/catalysis/finite_rate/finite_rate.cpp
    src/catalysis/finite_rate/finite_rate_input_parser.cpp
    src/catalysis/finite_rate/finite_rate_reaction_validation.cpp
    src/catalysis/finite_rate/finite_rate_reaction_properties.cpp
    src/catalysis/finite_rate/finite_rate_flux_computation.cpp
    src/catalysis/finite_rate/steady_surface_coverage.cpp
    src/catalysis/finite_rate/stoichiometry.cpp
    src/catalysis/gamma/gamma_input_parser.cpp
    src/catalysis/gamma/gamma_reaction_validation.cpp
    src/utils/input_models.cpp
  src/gasp2.cpp           # rimuovi questa riga se non lo usi
)

add_library(gasp2 ${GASP2_SOURCES})

target_include_directories(gasp2
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(gasp2 PRIVATE Eigen3::Eigen)

if(GASP2_ENABLE_WARNINGS)
  if(MSVC)
    target_compile_options(gasp2 PRIVATE /W4 /permissive- /Zc:preprocessor /utf-8)
  else()
    target_compile_options(gasp2 PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()

# --- Main executable ---
add_executable(gasp2_run
  tests/main.cpp
)
target_link_libraries(gasp2_run PRIVATE gasp2)

# --- Test ---
enable_testing()
include(CTest)
add_executable(gasp2_test
  tests/gamma_model/test_gasp2.cpp
)
target_link_libraries(gasp2_test PRIVATE gasp2)

file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/inputs" INPUTS_DIR_CMAKE)
target_compile_definitions(gasp2
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
target_compile_definitions(gasp2_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
target_compile_definitions(gasp2_run
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)

add_test(NAME gasp2_test COMMAND gasp2_test)

add_executable(gamma_validation_test
  tests/gamma_model/test_invalid_gamma.cpp
)
target_link_libraries(gamma_validation_test PRIVATE gasp2)
target_compile_definitions(gamma_validation_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME gamma_validation_test COMMAND gamma_validation_test)


add_executable(mixed_gamma_test
  tests/gamma_model/test_mixed_gamma.cpp
)
target_link_libraries(mixed_gamma_test PRIVATE gasp2)
target_compile_definitions(mixed_gamma_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME mixed_gamma_test COMMAND mixed_gamma_test)

add_executable(mass_balance_test
  tests/gamma_model/test_mass_balance.cpp
)
target_link_libraries(mass_balance_test PRIVATE gasp2)
target_compile_definitions(mass_balance_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME mass_balance_test COMMAND mass_balance_test)

add_executable(finite_rate_test
  tests/finite_rate/test_finite_rate.cpp
)
target_link_libraries(finite_rate_test PRIVATE gasp2)
target_compile_definitions(finite_rate_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_test COMMAND finite_rate_test)

add_executable(finite_rate_input_validation_test
  tests/finite_rate/test_input_validation.cpp
)
target_link_libraries(finite_rate_input_validation_test PRIVATE gasp2)
target_compile_definitions(finite_rate_input_validation_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_input_validation_test COMMAND finite_rate_input_validation_test)

add_executable(finite_rate_ion_test
  tests/finite_rate/test_ion_rejection.cpp
)
target_link_libraries(finite_rate_ion_test PRIVATE gasp2)
target_compile_definitions(finite_rate_ion_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_ion_test COMMAND finite_rate_ion_test)

add_executable(finite_rate_missing_adsorption_test
  tests/finite_rate/test_missing_adsorption.cpp
)
target_link_libraries(finite_rate_missing_adsorption_test PRIVATE gasp2)
target_compile_definitions(finite_rate_missing_adsorption_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_missing_adsorption_test COMMAND finite_rate_missing_adsorption_test)

add_executable(finite_rate_fractional_coeff_test
  tests/finite_rate/test_fractional_coeff.cpp
)
target_link_libraries(finite_rate_fractional_coeff_test PRIVATE gasp2)
target_compile_definitions(finite_rate_fractional_coeff_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_fractional_coeff_test COMMAND finite_rate_fractional_coeff_test)

add_executable(finite_rate_missing_parameter_test
  tests/finite_rate/test_missing_parameter.cpp
)
target_link_libraries(finite_rate_missing_parameter_test PRIVATE gasp2)
target_compile_definitions(finite_rate_missing_parameter_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_missing_parameter_test COMMAND finite_rate_missing_parameter_test)

add_executable(finite_rate_useless_parameter_test
  tests/finite_rate/test_useless_parameter.cpp
)
target_link_libraries(finite_rate_useless_parameter_test PRIVATE gasp2)
target_compile_definitions(finite_rate_useless_parameter_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_useless_parameter_test COMMAND finite_rate_useless_parameter_test)

add_executable(finite_rate_valid_adsorption_test
  tests/finite_rate/test_valid_adsorption.cpp
)
target_link_libraries(finite_rate_valid_adsorption_test PRIVATE gasp2)
target_compile_definitions(finite_rate_valid_adsorption_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_valid_adsorption_test COMMAND finite_rate_valid_adsorption_test)

add_executable(finite_rate_er_bare_site_test
  tests/finite_rate/test_er_bare_site.cpp
)
target_link_libraries(finite_rate_er_bare_site_test PRIVATE gasp2)
target_compile_definitions(finite_rate_er_bare_site_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME finite_rate_er_bare_site_test COMMAND finite_rate_er_bare_site_test)

add_executable(invalid_surface_model_test
  tests/gamma_model/test_invalid_surface_model.cpp
)
target_link_libraries(invalid_surface_model_test PRIVATE gasp2)
target_compile_definitions(invalid_surface_model_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME invalid_surface_model_test COMMAND invalid_surface_model_test)

add_executable(non_catalytic_test
  tests/gamma_model/test_non_catalytic.cpp
)
target_link_libraries(non_catalytic_test PRIVATE gasp2)
target_compile_definitions(non_catalytic_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME non_catalytic_test COMMAND non_catalytic_test)

add_executable(carbon_gamma_test
  tests/gamma_model/test_carbon_gamma.cpp
)
target_link_libraries(carbon_gamma_test PRIVATE gasp2)
target_compile_definitions(carbon_gamma_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME carbon_gamma_test COMMAND carbon_gamma_test)

add_executable(gamma_consistent_test
  tests/gamma_model/test_gamma_consistent.cpp
)
target_link_libraries(gamma_consistent_test PRIVATE gasp2)
target_compile_definitions(gamma_consistent_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME gamma_consistent_test COMMAND gamma_consistent_test)

add_executable(duplicate_reaction_test
  tests/gamma_model/test_duplicate_reaction.cpp
)
target_link_libraries(duplicate_reaction_test PRIVATE gasp2)
target_compile_definitions(duplicate_reaction_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME duplicate_reaction_test COMMAND duplicate_reaction_test)

add_executable(trivial_reaction_test
  tests/gamma_model/test_trivial_reaction.cpp
)
target_link_libraries(trivial_reaction_test PRIVATE gasp2)
target_compile_definitions(trivial_reaction_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME trivial_reaction_test COMMAND trivial_reaction_test)

add_executable(bose_gamma_test
  tests/gamma_model/test_bose_gamma.cpp
)
target_link_libraries(bose_gamma_test PRIVATE gasp2)
target_compile_definitions(bose_gamma_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME bose_gamma_test COMMAND bose_gamma_test)

add_executable(supercatalytic_test
  tests/gamma_model/test_supercatalytic.cpp
)
target_link_libraries(supercatalytic_test PRIVATE gasp2)
target_compile_definitions(supercatalytic_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME supercatalytic_test COMMAND supercatalytic_test)

add_executable(rini_gamma_test
  tests/gamma_model/test_rini_gamma.cpp
)
target_link_libraries(rini_gamma_test PRIVATE gasp2)
target_compile_definitions(rini_gamma_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME rini_gamma_test COMMAND rini_gamma_test)

add_executable(rini_gamma_three_test
  tests/gamma_model/test_rini_gamma_three.cpp
)
target_link_libraries(rini_gamma_three_test PRIVATE gasp2)
target_compile_definitions(rini_gamma_three_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME rini_gamma_three_test COMMAND rini_gamma_three_test)

add_executable(ion_reaction_test
  tests/gamma_model/test_ion_reaction.cpp
)
target_link_libraries(ion_reaction_test PRIVATE gasp2)
target_compile_definitions(ion_reaction_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME ion_reaction_test COMMAND ion_reaction_test)

add_executable(charge_balance_test
  tests/gamma_model/test_charge_balance.cpp
)
target_link_libraries(charge_balance_test PRIVATE gasp2)
target_compile_definitions(charge_balance_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME charge_balance_test COMMAND charge_balance_test)

add_executable(invalid_arrow_test
  tests/gamma_model/test_invalid_arrow.cpp
)
target_link_libraries(invalid_arrow_test PRIVATE gasp2)
target_compile_definitions(invalid_arrow_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME invalid_arrow_test COMMAND invalid_arrow_test)
add_executable(electron_reaction_test
  tests/gamma_model/test_electron_reaction.cpp
)
target_link_libraries(electron_reaction_test PRIVATE gasp2)
target_compile_definitions(electron_reaction_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME electron_reaction_test COMMAND electron_reaction_test)

add_executable(input_validation_test
  tests/gamma_model/test_input_validation.cpp
)
target_link_libraries(input_validation_test PRIVATE gasp2)
target_compile_definitions(input_validation_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME input_validation_test COMMAND input_validation_test)

add_executable(limiting_flux_test
  tests/gamma_model/test_limiting_flux.cpp
)
target_link_libraries(limiting_flux_test PRIVATE gasp2)
target_compile_definitions(limiting_flux_test
  PRIVATE INPUTS_DIR="${INPUTS_DIR_CMAKE}"
)
add_test(NAME limiting_flux_test COMMAND limiting_flux_test)
