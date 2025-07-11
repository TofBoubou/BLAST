# BLAST Modern Makefile with Output System
# Updated to include HDF5 and comprehensive output capabilities

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++23 -O3 -DNDEBUG -march=native -w
DEBUG_FLAGS = -std=c++23 -Wall -Wextra -g -O0 -DDEBUG
INCLUDE_FLAGS = -Iinclude -I/usr/include/hdf5/serial


# Library paths and flags
EIGEN_PATH = libs/eigen
MUTATIONPP_PATH = libs/mutationpp
YAML_CPP_PATH = libs/yaml-cpp

# External library flags
HDF5_FLAGS = $(shell pkg-config --cflags --libs hdf5)
ZLIB_FLAGS = -lz

# Include paths
INCLUDES = $(INCLUDE_FLAGS) \
           -I$(EIGEN_PATH) \
           -I$(MUTATIONPP_PATH)/install/include \
           -I$(YAML_CPP_PATH)/include

# Library linking
LIBS = -L$(MUTATIONPP_PATH)/install/lib -lmutation++ \
       -L$(YAML_CPP_PATH)/build -lyaml-cpp \
       $(HDF5_FLAGS) \
       $(ZLIB_FLAGS) \
       -lpthread -ldl

# Directories
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = .

# Source files organized by module
BOUNDARY_LAYER_SOURCES = \
    $(SRC_DIR)/boundary_layer/coefficients/coefficient_calculator.cpp \
    $(SRC_DIR)/boundary_layer/coefficients/diffusion.cpp \
    $(SRC_DIR)/boundary_layer/coefficients/xi_derivatives.cpp \
    $(SRC_DIR)/boundary_layer/conditions/boundary_interpolator.cpp \
    $(SRC_DIR)/boundary_layer/equations/continuity.cpp \
    $(SRC_DIR)/boundary_layer/equations/energy.cpp \
    $(SRC_DIR)/boundary_layer/equations/momentum.cpp \
    $(SRC_DIR)/boundary_layer/equations/species.cpp \
    $(SRC_DIR)/boundary_layer/grid/coordinate_transform.cpp \
    $(SRC_DIR)/boundary_layer/grid/grid.cpp \
    $(SRC_DIR)/boundary_layer/solver/boundary_layer_solver.cpp \
    $(SRC_DIR)/boundary_layer/solvers/tridiagonal_solvers.cpp \
    $(SRC_DIR)/boundary_layer/thermodynamics/enthalpy_temperature_solver.cpp

IO_SOURCES = \
    $(SRC_DIR)/io/config_manager.cpp \
    $(SRC_DIR)/io/yaml_parser.cpp

# Output system sources
OUTPUT_SOURCES = \
    $(SRC_DIR)/io/output/output_writer.cpp \
    $(SRC_DIR)/io/output/hdf5_writer.cpp

THERMOPHYSICS_SOURCES = \
    $(SRC_DIR)/thermophysics/mutation_mixture.cpp

MAIN_SOURCE = $(SRC_DIR)/main.cpp

# All sources
ALL_SOURCES = $(BOUNDARY_LAYER_SOURCES) $(IO_SOURCES) $(OUTPUT_SOURCES) $(THERMOPHYSICS_SOURCES) $(MAIN_SOURCE)

# Object files
BOUNDARY_LAYER_OBJECTS = $(BOUNDARY_LAYER_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
IO_OBJECTS = $(IO_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
OUTPUT_OBJECTS = $(OUTPUT_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
THERMOPHYSICS_OBJECTS = $(THERMOPHYSICS_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
MAIN_OBJECT = $(MAIN_SOURCE:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

ALL_OBJECTS = $(BOUNDARY_LAYER_OBJECTS) $(IO_OBJECTS) $(OUTPUT_OBJECTS) $(THERMOPHYSICS_OBJECTS) $(MAIN_OBJECT)

# Targets
TARGET = $(BIN_DIR)/blast
DEBUG_TARGET = $(BIN_DIR)/blast_debug

# Default target
.PHONY: all
all: $(TARGET)

# Debug target
.PHONY: debug
debug: CXXFLAGS = $(DEBUG_FLAGS)
debug: $(DEBUG_TARGET)

# Production target
$(TARGET): $(ALL_OBJECTS) | check_dependencies
	@echo "Linking production executable..."
	$(CXX) $(ALL_OBJECTS) -o $@ $(LIBS)
	@echo "✓ Production build complete: $@"

# Debug target
$(DEBUG_TARGET): $(ALL_OBJECTS) | check_dependencies
	@echo "Linking debug executable..."
	$(CXX) $(ALL_OBJECTS) -o $@ $(LIBS)
	@echo "✓ Debug build complete: $@"

# Object file compilation rules
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create build directories
$(BUILD_DIR):
	@echo "Creating build directories..."
	@mkdir -p $(BUILD_DIR)/boundary_layer/coefficients
	@mkdir -p $(BUILD_DIR)/boundary_layer/conditions
	@mkdir -p $(BUILD_DIR)/boundary_layer/equations
	@mkdir -p $(BUILD_DIR)/boundary_layer/grid
	@mkdir -p $(BUILD_DIR)/boundary_layer/solver
	@mkdir -p $(BUILD_DIR)/boundary_layer/solvers
	@mkdir -p $(BUILD_DIR)/boundary_layer/thermodynamics
	@mkdir -p $(BUILD_DIR)/io
	@mkdir -p $(BUILD_DIR)/io/output
	@mkdir -p $(BUILD_DIR)/thermophysics

# Dependency checking
.PHONY: check_dependencies
check_dependencies:
	@echo "Checking dependencies..."
	@command -v pkg-config >/dev/null 2>&1 || { echo "Error: pkg-config not found"; exit 1; }
	@pkg-config --exists hdf5 || { echo "Error: HDF5 development libraries not found. Install libhdf5-dev"; exit 1; }
	@test -d $(MUTATIONPP_PATH) || { echo "Error: Mutation++ not found at $(MUTATIONPP_PATH)"; exit 1; }
	@test -d $(YAML_CPP_PATH) || { echo "Error: yaml-cpp not found at $(YAML_CPP_PATH)"; exit 1; }
	@test -d $(EIGEN_PATH) || { echo "Error: Eigen not found at $(EIGEN_PATH)"; exit 1; }
	@echo "✓ All dependencies found"

# Library building
.PHONY: build_libs
build_libs:
	@echo "Building external libraries..."
	@if [ ! -f $(MUTATIONPP_PATH)/install/lib/libmutation++.a ]; then \
		echo "Building Mutation++..."; \
		cd $(MUTATIONPP_PATH) && mkdir -p build && cd build && \
		cmake .. -DCMAKE_INSTALL_PREFIX=../install && \
		make -j$(shell nproc) && make install; \
	fi
	@if [ ! -f $(YAML_CPP_PATH)/build/libyaml-cpp.a ]; then \
		echo "Building yaml-cpp..."; \
		cd $(YAML_CPP_PATH) && mkdir -p build && cd build && \
		cmake .. -DYAML_CPP_BUILD_TESTS=OFF && \
		make -j$(shell nproc); \
	fi
	@echo "✓ External libraries built"

# System dependency installation (Ubuntu/Debian)
.PHONY: install_deps_ubuntu
install_deps_ubuntu:
	@echo "Installing system dependencies (Ubuntu/Debian)..."
	sudo apt-get update
	sudo apt-get install -y \
		build-essential \
		cmake \
		pkg-config \
		libhdf5-dev \
		libboost-all-dev \
		zlib1g-dev \
		git
	@echo "✓ System dependencies installed"

# System dependency installation (macOS)
.PHONY: install_deps_macos
install_deps_macos:
	@echo "Installing system dependencies (macOS)..."
	brew install \
		cmake \
		pkg-config \
		hdf5 \
		boost \
		zlib
	@echo "✓ System dependencies installed"

# Testing targets
.PHONY: test
test: $(TARGET)
	@echo "Running tests..."
	@if [ -f config/default.yaml ]; then \
		./$(TARGET) config/default.yaml test_output; \
	else \
		echo "No test configuration found. Create config/default.yaml"; \
	fi

# Post-processing setup
.PHONY: setup_postprocess
setup_postprocess:
	@echo "Setting up HDF5 post-processing environment..."
	@command -v python3 >/dev/null 2>&1 || { echo "Error: Python 3 not found"; exit 1; }
	python3 -m pip install numpy matplotlib pandas h5py seaborn pathlib
	@echo "✓ HDF5 post-processing setup complete"

# Example run
.PHONY: example
example: $(TARGET)
	@echo "Running example simulation..."
	@mkdir -p example_output
	./$(TARGET) config/default.yaml example_simulation
	@echo "✓ Example complete. HDF5 output available for analysis"

# Performance profiling
.PHONY: profile
profile: debug
	@echo "Running performance profile..."
	@command -v valgrind >/dev/null 2>&1 || { echo "Install valgrind for profiling"; exit 1; }
	valgrind --tool=callgrind --callgrind-out-file=profile.out ./$(DEBUG_TARGET) config/default.yaml profile_test
	@echo "Profile saved to profile.out"

# Memory check
.PHONY: memcheck
memcheck: debug
	@echo "Running memory check..."
	@command -v valgrind >/dev/null 2>&1 || { echo "Install valgrind for memory checking"; exit 1; }
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./$(DEBUG_TARGET) config/default.yaml memcheck_test

# Documentation generation
.PHONY: docs
docs:
	@echo "Generating documentation..."
	@command -v doxygen >/dev/null 2>&1 || { echo "Install doxygen for documentation"; exit 1; }
	doxygen Doxyfile
	@echo "✓ Documentation generated in docs/"

# Static analysis
.PHONY: analyze
analyze:
	@echo "Running static analysis..."
	@command -v cppcheck >/dev/null 2>&1 || { echo "Install cppcheck for static analysis"; exit 1; }
	cppcheck --enable=all --std=c++23 --inconclusive --xml --xml-version=2 \
		--suppress=missingIncludeSystem \
		$(SRC_DIR) 2> analysis.xml
	@echo "Static analysis complete. Results in analysis.xml"

# Format code
.PHONY: format
format:
	@echo "Formatting code..."
	@command -v clang-format >/dev/null 2>&1 || { echo "Install clang-format for code formatting"; exit 1; }
	find $(SRC_DIR) include -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i
	@echo "✓ Code formatted"

# Clean targets
.PHONY: clean
clean:
	@echo "Cleaning build files..."
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET) $(DEBUG_TARGET)
	@echo "✓ Clean complete"

.PHONY: clean_all
clean_all: clean
	@echo "Cleaning all generated files..."
	rm -rf example_output/
	rm -f *.out *.xml
	rm -rf docs/
	@echo "✓ Deep clean complete"

# Install target
.PHONY: install
install: $(TARGET)
	@echo "Installing BLAST..."
	@mkdir -p $(HOME)/bin
	cp $(TARGET) $(HOME)/bin/
	@echo "✓ BLAST installed to $(HOME)/bin"

# Help target
.PHONY: help
help:
	@echo "BLAST Modern Build System"
	@echo "========================="
	@echo ""
	@echo "Main targets:"
	@echo "  all              Build production executable"
	@echo "  debug            Build debug executable"
	@echo "  test             Run test simulation"
	@echo "  example          Run example simulation"
	@echo "  clean            Remove build files"
	@echo "  clean_all        Remove all generated files"
	@echo ""
	@echo "Setup targets:"
	@echo "  install_deps_ubuntu   Install Ubuntu/Debian dependencies"
	@echo "  install_deps_macos    Install macOS dependencies"
	@echo "  build_libs            Build external libraries"
	@echo "  setup_postprocess     Setup HDF5 post-processing"
	@echo ""
	@echo "Development targets:"
	@echo "  format           Format source code"
	@echo "  analyze          Run static analysis"
	@echo "  profile          Profile performance"
	@echo "  memcheck         Check for memory leaks"
	@echo "  docs             Generate documentation"
	@echo ""
	@echo "Example workflow:"
	@echo "  make install_deps_ubuntu  # Install system dependencies"
	@echo "  make build_libs           # Build external libraries"
	@echo "  make setup_postprocess    # Setup HDF5 post-processing"
	@echo "  make all                  # Build BLAST"
	@echo "  make example              # Run example"

# Print configuration
.PHONY: config
config:
	@echo "Build Configuration:"
	@echo "==================="
	@echo "CXX: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "INCLUDES: $(INCLUDES)"
	@echo "LIBS: $(LIBS)"
	@echo "TARGET: $(TARGET)"
	@echo ""
	@echo "Library paths:"
	@echo "  Eigen: $(EIGEN_PATH)"
	@echo "  Mutation++: $(MUTATIONPP_PATH)"
	@echo "  yaml-cpp: $(YAML_CPP_PATH)"
	@echo ""
	@echo "Source files: $(words $(ALL_SOURCES)) files"
	@echo "Object files: $(words $(ALL_OBJECTS)) objects"

# Automatic dependency generation
-include $(ALL_OBJECTS:.o=.d)

$(BUILD_DIR)/%.d: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MM -MT $(@:.d=.o) $< > $@

.PHONY: deps
deps: $(ALL_OBJECTS:.o=.d)

# Default target
.DEFAULT_GOAL := help