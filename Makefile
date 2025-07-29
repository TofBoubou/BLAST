# BLAST Modern Makefile with Output System
# Updated to include HDF5 and comprehensive output capabilities

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++23 -O3 -DNDEBUG -march=native -w
SHELL := /bin/bash
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

# Automatic source detection
ALL_SOURCES := $(shell find $(SRC_DIR) -name '*.cpp' -type f)
TOTAL_FILES := $(words $(ALL_SOURCES))

# Object files
ALL_OBJECTS = $(ALL_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Target
TARGET = $(BIN_DIR)/blast

# Default target
.PHONY: all
all:
	@echo "Building BLAST..."
	@rm -f $(BUILD_DIR)/.progress 2>/dev/null || true
	@$(MAKE) -s $(TARGET)

# Build target
$(TARGET): $(ALL_OBJECTS) | check_dependencies
	@printf "▸ Linking         "
	@$(CXX) $(ALL_OBJECTS) -o $@ $(LIBS) 2>/dev/null && printf "✓\n" || (printf "✗\n" && exit 1)
	@printf "\n✓ Complete → %s\n" "$@"

# Object file compilation rules with progress
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	@to_build=0; \
	for src in $(ALL_SOURCES); do \
		obj=$$(echo "$$src" | sed 's|$(SRC_DIR)/|$(BUILD_DIR)/|; s|\.cpp$$|.o|'); \
		if [ "$$src" -nt "$$obj" ] 2>/dev/null || [ ! -f "$$obj" ]; then \
			to_build=$$((to_build + 1)); \
		fi; \
	done; \
	echo "$$$$" >> $(BUILD_DIR)/.progress 2>/dev/null || echo "$$$$" > $(BUILD_DIR)/.progress; \
	current=$$(wc -l < $(BUILD_DIR)/.progress 2>/dev/null); \
	percent=$$((current * 100 / to_build)); \
	filled=$$((percent / 4)); \
	bar=$$(printf "%*s" $$filled | tr ' ' '='); \
	printf "\r[%-25s] %3d%% (%d/%d) %s\n" \
		"$$bar" $$percent $$current $$to_build "$(notdir $<)"; \
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ 2>/dev/null

# Create build directories automatically
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Dependency checking
.PHONY: check_dependencies
check_dependencies:
	@printf "▸ Dependencies    "
	@command -v pkg-config >/dev/null 2>&1 || { printf "✗\nError: pkg-config not found\n"; exit 1; }
	@pkg-config --exists hdf5 || { printf "✗\nError: HDF5 development libraries not found. Install libhdf5-dev\n"; exit 1; }
	@test -d $(MUTATIONPP_PATH) || { printf "✗\nError: Mutation++ not found at $(MUTATIONPP_PATH)\n"; exit 1; }
	@test -d $(YAML_CPP_PATH) || { printf "✗\nError: yaml-cpp not found at $(YAML_CPP_PATH)\n"; exit 1; }
	@test -d $(EIGEN_PATH) || { printf "✗\nError: Eigen not found at $(EIGEN_PATH)\n"; exit 1; }
	@printf "✓\n"

# Library building
.PHONY: build_libs
build_libs:
	@echo "Building external libraries..."
	@if [ ! -f $(MUTATIONPP_PATH)/install/lib/libmutation++.a ]; then \
		echo "Building Mutation++..."; \
		cd $(MUTATIONPP_PATH) && mkdir -p build && cd build && \
		cmake .. -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_BUILD_TYPE=Release && \
		make -j$(shell nproc) && make install; \
	else \
		echo "✓ Mutation++ already built"; \
	fi
	@if [ ! -f $(YAML_CPP_PATH)/build/libyaml-cpp.a ]; then \
		echo "Building yaml-cpp..."; \
		cd $(YAML_CPP_PATH) && mkdir -p build && cd build && \
		cmake .. -DYAML_CPP_BUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release && \
		make -j$(shell nproc); \
	else \
		echo "✓ yaml-cpp already built"; \
	fi
	@echo "✓ External libraries ready"

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
	rm -f $(TARGET)
	@echo "✓ Clean complete"


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
	@echo "  clean            Remove build files"
	@echo ""
	@echo "Setup targets:"
	@echo "  install_deps_ubuntu   Install Ubuntu/Debian dependencies"
	@echo "  install_deps_macos    Install macOS dependencies"
	@echo "  build_libs            Build external libraries"
	@echo ""
	@echo "Development targets:"
	@echo "  format           Format source code"
	@echo ""
	@echo "Example workflow:"
	@echo "  make install_deps_ubuntu  # Install system dependencies"
	@echo "  make build_libs           # Build external libraries"
	@echo "  make all                  # Build BLAST"

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

# Dependency generation (disabled for cleaner build)
# $(BUILD_DIR)/%.d: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
#	@mkdir -p $(dir $@)
#	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MM -MT $(@:.d=.o) $< > $@
# -include $(ALL_OBJECTS:.o=.d)

# Default target
.DEFAULT_GOAL := all