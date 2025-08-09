# BLAST Modern

A high-performance C++ boundary layer solver for hypersonic flows with chemical non-equilibrium effects.

## Overview

BLAST is a computational fluid dynamics solver specifically designed for analyzing boundary layer flows. The solver incorporates advanced thermophysical modeling including chemical non-equilibrium effects, making it particularly suitable for atmospheric entry vehicle analysis and hypersonic flow studies.

**Key Features:**
- Boundary layer equations solving for axisymmetric geometries
- Chemical non-equilibrium modeling with finite-rate chemistry
- Stefan-Maxwell multicomponent diffusion
- Comprehensive thermodynamic property calculations via Mutation++
- HDF5-based output system with complete flow field data
- YAML-driven configuration system
- Modern C++23 implementation with optimized performance

## Features

### Flow Physics
- **Boundary Layer Equations**: Full Navier-Stokes boundary layer equations
- **Chemical Non-Equilibrium**: Finite-rate chemistry with species transport
- **Diffusion Models**: Stefan-Maxwell multicomponent diffusion
- **Thermal Effects**: Energy equation with chemical energy coupling

### Numerical Methods
- **Finite Difference Discretization**: High-order accurate schemes
- **Stagnation point method**: Main numerical methode to ensure convergence
- **Adaptive Relaxation**: Automatic convergence acceleration
- **Tridiagonal Solvers**: Efficient linear system solution
- **Grid Transformation**: Coordinate transformation

### Thermophysics
- **Mutation++ Integration**: State-of-the-art thermophysical property library
- **Multiple Databases**: RRHO, NASA polynomials, and other thermodynamic models
- **Transport Properties**: Chapman-Enskog theory
- **Mixture Models**: Support for arbitrary gas mixtures

## Requirements

### System Requirements
- **Compiler**: GCC 11+ or Clang 13+ (C++23 support required)
- **Build System**: Make
- **Platform**: Linux, macOS, Windows (WSL recommended)

### Dependencies
- **Eigen 3.4+**: Linear algebra library (must be cloned)
- **Mutation++ 1.0+**: Thermophysical properties library (must be cloned)
- **yaml-cpp 0.7+**: YAML parsing library (must be cloned)
- **HDF5 1.10+**: High-performance data format library
- **pkg-config**: Dependency management
- **Zlib**: Compression library

## Installation

### Quick Start

1. **Clone the repository:**
```bash
git clone https://github.com/TofBoubou/BLAST.git
cd BLAST
```

2. **Clone external libraries:**
```bash
# Clone required libraries into libs/ directory
git clone https://gitlab.com/libeigen/eigen.git libs/eigen
git clone https://github.com/mutationpp/Mutationpp.git libs/mutationpp
git clone https://github.com/jbeder/yaml-cpp.git libs/yaml-cpp
```

3. **Install system dependencies:**

**Ubuntu/Debian:**
```bash
make install_deps_ubuntu
```

**macOS:**
```bash
make install_deps_macos
```

4. **Build external libraries:**
```bash
make build_libs
```

5. **Compile BLAST:**
```bash
make all
```

### Detailed Installation

#### System Dependencies

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    pkg-config \
    libhdf5-dev \
    libboost-all-dev \
    zlib1g-dev \
    git
```

**macOS:**
```bash
brew install cmake pkg-config hdf5 boost zlib
```

#### Building External Libraries

The project includes all necessary libraries as git submodules. Build them with:

```bash
# Build Mutation++ thermophysics library
cd libs/mutationpp
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_BUILD_TYPE=Release
make -j$(nproc) && make install
cd ../../..

# Build yaml-cpp configuration parser
cd libs/yaml-cpp
mkdir build && cd build
cmake .. -DYAML_CPP_BUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
cd ../../..
```

Or use the automated build:
```bash
make build_libs
```

#### Compilation Options

**Production build (default):**
```bash
make all
```

**Debug build with profiling:**
```bash
make profile
```

**Check build configuration:**
```bash
make config
```

## Usage

### Basic Command Line

```bash
./blast <config_file.yaml>
```

**Examples:**
```bash
# Run with standard CO2 configuration
./blast config/CO2_5.yaml

# Run with custom output name
./blast config/config_new_standard.yaml
```

### Configuration Files

BLAST uses YAML configuration files to specify simulation parameters. Two example configurations are provided:

- `config/CO2_5.yaml`: CO2 mixture simulation
- `config/config_new_standard.yaml`: Standard configuration template

### Output

The solver generates HDF5 output files containing:
- Complete flow field data (velocity, temperature, species concentrations)
- Boundary layer profiles at specified stations
- Thermodynamic properties
- Transport coefficients

Output files are saved as `simulation_YYYYMMDD_HHMMSS.h5` in the `test_outputs/` directory.

## Configuration

### YAML Configuration Structure

```yaml
simulation:
  body_type: axisymmetric              # Body geometry type
  only_stagnation_point: false        # Full body vs stagnation point only
  diffusion_type: stefan_maxwell      # Diffusion model
  consider_thermal_diffusion: false   # Thermal diffusion effects
  chemical_non_equilibrium: true      # Enable finite-rate chemistry

numerical:
  n_eta: 50                           # Grid points in normal direction
  eta_max: 6.0                        # Boundary layer edge location
  convergence_tolerance: 1.0e-06      # Convergence criterion
  max_iterations: 100000              # Maximum iterations

mixture:
  name: CO2_5                         # Mixture name (from Mutation++ database)
  thermodynamic_database: RRHO        # Thermodynamic model
  viscosity_algorithm: chapman_enskog_ldlt  # Transport property method

output:
  x_stations: [0.0, 0.002, 0.004, ...]  # Output stations along body
```

### Key Parameters

#### Simulation Parameters
- `body_type`: Geometry type (`axisymmetric`)
- `only_stagnation_point`: Boolean for stagnation point analysis only
- `diffusion_type`: Diffusion model (`stefan_maxwell`, `fick`)
- `chemical_non_equilibrium`: Enable chemical reactions

#### Numerical Parameters
- `n_eta`: Number of grid points across boundary layer
- `eta_max`: Boundary layer edge location (similarity coordinate)
- `convergence_tolerance`: Residual tolerance for convergence
- `max_iterations`: Maximum solver iterations

#### Mixture Properties
- `name`: Gas mixture name (must exist in Mutation++ database)
- `thermodynamic_database`: Property calculation method
- `viscosity_algorithm`: Transport property calculation method

## Examples

### Example 1: CO2 Stagnation Point Flow

```bash
# Run CO2 simulation
./blast config/CO2_5.yaml

# Results will be in test_outputs/simulation_*.h5
```

### Post-Processing

BLAST includes a comprehensive Python post-processing script that generates plots from HDF5 simulation results.

**Basic usage:**
```bash
cd scripts
python3 postprocess_blast.py --input ../test_outputs/simulation_YYYYMMDD_HHMMSS.h5 --plots all
```

**Available plot types:**
```bash
# Generate all plots (profiles for up to 5 stations + F/g map + summary JSON)
python3 postprocess_blast.py --input simulation.h5 --plots all

# Boundary layer profiles at a specific station
python3 postprocess_blast.py --input simulation.h5 --plots profiles --station 0

# Generate only the 2D interpolated F(η, ξ) and g(η, ξ) maps
python3 postprocess_blast.py --input simulation.h5 --plots f_g_map

# Generate all plots and export summary.json (same as 'all', explicit)
python3 postprocess_blast.py --input simulation.h5 --plots summary
```

**Generated outputs:**
- Boundary layer velocity, temperature, and species profiles
- Heat flux and wall properties
- Flow field contour plots
- Convergence history
- Summary statistics in JSON format

## Development

### Project Structure

```
BLAST_MODERN/
├── include/blast/           # Header files
│   ├── boundary_layer/      # Core solver components
│   ├── core/               # Utility classes
│   ├── io/                 # Input/output systems
│   └── thermophysics/      # Thermophysical interfaces
├── src/                    # Source files
├── libs/                   # External libraries (submodules)
│   ├── eigen/              # Linear algebra
│   ├── mutationpp/         # Thermophysics
│   └── yaml-cpp/           # Configuration parsing
├── config/                 # Example configurations
├── scripts/                # Post-processing tools
└── test_outputs/           # Simulation results
```

### Building and Testing

**Format code:**
```bash
make format
```

**Performance profiling:**
```bash
make profile CONFIG=config/CO2_5.yaml
# Results in profile_report.txt
```

**Clean builds:**
```bash
make clean          # Clean all build files
make clean-profile  # Clean profile files only
```

### Dependencies

This project builds upon several excellent open-source libraries:

- **[Mutation++](https://github.com/mutationpp/Mutationpp)**: Thermophysical properties library (LGPL v3)
- **[Eigen](https://eigen.tuxfamily.org/)**: Linear algebra library (MPL2)
- **[yaml-cpp](https://github.com/jbeder/yaml-cpp)**: YAML parsing library (MIT)
- **[HDF5](https://www.hdfgroup.org/solutions/hdf5/)**: High-performance data format (BSD-style)

**BLAST** - High-Performance Boundary Layer Analysis for Hypersonic Flows
