# BLAST Modern

A high-performance C++ boundary layer solver for hypersonic flows with chemical non-equilibrium effects.

## Overview

BLAST is a computational fluid dynamics solver specifically designed for analyzing boundary layer flows. The solver incorporates advanced thermophysical modeling including chemical non-equilibrium effects, making it particularly suitable for atmospheric entry vehicle analysis and hypersonic flow studies.

**Key Features:**
- Boundary layer equations solving for axisymmetric geometries
- Chemical non-equilibrium modeling with finite-rate chemistry
- Stefan-Maxwell multicomponent diffusion with thermal effects
- Catalytic wall modeling for surface reactions
- Numerical continuation methods for challenging convergence cases
- Automated abaque generation for thermal protection system design
- Comprehensive thermodynamic property calculations via Mutation++
- HDF5-based output system with complete flow field data
- YAML-driven configuration system
- Modern C++23 implementation with optimized performance

## Features

### Flow Physics
- **Boundary Layer Equations**: Full Navier-Stokes boundary layer equations
- **Chemical Non-Equilibrium**: Finite-rate chemistry with species transport
- **Advanced Diffusion Models**: Stefan-Maxwell multicomponent diffusion
- **Thermal Diffusion Effects**: Soret effect (thermal diffusion) and Dufour effect
- **Catalytic Walls**: Surface chemical reactions with controllable catalyticity
- **Multiple Chemical Modes**: Equilibrium, frozen, or non-equilibrium chemistry

### Numerical Methods
- **Finite Difference Discretization**: High-order accurate schemes
- **Stagnation Point Method**: Main numerical method to ensure convergence
- **Continuation Method**: Advanced numerical stabilization for difficult cases
- **Adaptive Relaxation**: Automatic convergence acceleration
- **Tridiagonal Solvers**: Efficient linear system solution

### Engineering Tools
- **Abaque Generation**: Automated heat flux maps for thermal protection design
- **Multi-Station Analysis**: Complete boundary layer development along vehicle surface
- **Comprehensive Post-Processing**: Specialized Python tools for data analysis and visualization

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

## Usage

### Basic Command Line

```bash
./blast <config_file.yaml>
```

**Examples:**
```bash
# Run with standard CO2 configuration
./blast config/CO2_5.yaml

# Run with custom configuration
./blast config/my_simulation.yaml
```

### Advanced Physical Modeling

#### Thermal Diffusion Effects
BLAST supports both thermal diffusion (Soret effect) and Dufour effect for accurate modeling of species transport in high-temperature flows:

```yaml
simulation:
  consider_thermal_diffusion: true   # Species migration due to temperature gradients
  consider_dufour_effect: true       # Energy transport due to concentration gradients
```

These effects are particularly important for hypersonic flows where strong temperature gradients exist.

#### Catalytic Wall Modeling
For atmospheric entry applications, surface catalysis significantly affects heat transfer:

```yaml
simulation:
  catalytic_wall: true               # Enable surface chemistry
```

This enables finite-rate surface reactions where dissociated species can recombine at the wall, releasing chemical energy and increasing heat flux.

#### Chemical Modes
Choose the appropriate chemistry model for your application:

```yaml
simulation:
  chemical_mode: "non_equilibrium"   # Options: "equilibrium", "frozen", "non_equilibrium"
```

- **Equilibrium**: Instantaneous chemical equilibrium (fast chemistry limit)
- **Frozen**: No chemical reactions (slow chemistry limit)  
- **Non-equilibrium**: Finite-rate chemistry with species transport

### Continuation Method

For challenging cases with extreme conditions (high temperatures, pressures, or strong catalytic effects), BLAST includes a continuation method that helps achieve convergence:

```yaml
continuation:
  wall_temperature_stable: 3000.0   # Starting wall temperature [K]
  edge_temperature_stable: 5905.5   # Starting edge temperature [K]
  pressure_stable: 7000.0           # Starting pressure [Pa]
```

**How it works**: The solver first converges a "stable" case with moderate conditions, then gradually transitions to your target conditions. This prevents numerical instabilities that can occur when starting directly with extreme parameters.

**When to use**: 
- Very high wall temperatures (> 2500K)
- Strong catalytic walls
- High-pressure flows
- Cases where direct solving fails to converge

### Abaque Generation

BLAST can automatically generate abaques (heat flux maps) for thermal protection system design:

```yaml
abaque:
  enabled: true
  catalyticity_values: [0.0, 0.01, 0.1, 1.0]     # Range of catalytic efficiencies
  temperature_range: [300, 3000]                  # Wall temperature range [K]
  temperature_points: 50                          # Resolution
```

**Purpose**: Abaques provide heat flux as a function of wall temperature and catalyticity, essential for:
- Thermal protection system sizing
- Material selection
- Mission design optimization
- Parametric studies

**Output**: 2D maps showing how heat flux varies with wall conditions, allowing engineers to quickly assess thermal loads across different scenarios.

## Configuration

### Complete YAML Configuration Structure

```yaml
simulation:
  body_type: axisymmetric
  only_stagnation_point: true
  diffusion_type: stefan_maxwell
  consider_thermal_diffusion: false
  consider_dufour_effect: false
  chemical_mode: "non_equilibrium"     # "frozen" or "equilibrium" or "non_equilibrium"
  catalytic_wall: false

numerical:
  n_eta: 200                       
  eta_max: 6.0
  convergence_tolerance: 1.0e-9    
  max_iterations: 100000

mixture:
  name: CO2_5
  thermodynamic_database: RRHO
  viscosity_algorithm: chapman_enskog_ldlt
  state_model: ChemNonEq1T

output:
  x_stations: [0.0]
  output_directory: test_outputs

outer_edge:
  edge_points:
    - x: 0.0
      radius: 1
      velocity: 0
      temperature: 5905.5
      pressure: 7000
      boundary_override: false  # or true to use specific mass fractions
      mass_fraction_condition: [0.1, 0.2, 0.3, 0.2, 0.2]
  
  velocity_gradient_stagnation: 4027.1213622326
  freestream_density: 0.0019382429
  freestream_velocity: 500000

wall_parameters:
  temperatures: [300]

abaque:
  enabled: false
  catalyticity_values: [0.0, 0.01]
  temperature_range: [3500, 5500]
  temperature_points: 20

continuation:
  wall_temperature_stable: 3000.0
  edge_temperature_stable: 5905.5
  pressure_stable: 7000.0
```

### Key Configuration Sections

#### Simulation Parameters
- `body_type`: Geometry type (`axisymmetric`)
- `only_stagnation_point`: Boolean for stagnation point analysis only
- `chemical_mode`: Chemistry treatment (`equilibrium`, `frozen`, `non_equilibrium`)
- `consider_thermal_diffusion`: Enable Soret effect
- `consider_dufour_effect`: Enable Dufour effect
- `catalytic_wall`: Enable surface chemical reactions

#### Physical Effects Control
- `diffusion_type`: Multicomponent diffusion model (`stefan_maxwell`)
- `thermal_diffusion`/`dufour_effect`: Coupled heat and mass transfer
- `catalytic_wall`: Surface chemistry for atmospheric entry applications

#### Numerical Parameters
- `n_eta`: Number of grid points across boundary layer
- `eta_max`: Boundary layer edge location (similarity coordinate)
- `convergence_tolerance`: Residual tolerance for convergence
- `max_iterations`: Maximum solver iterations

## Post-Processing

BLAST includes comprehensive Python post-processing tools for analyzing simulation results:

### Core Visualization Script

```bash
cd scripts
python3 postprocess_blast.py --input ../test_outputs/simulation_YYYYMMDD_HHMMSS.h5 --plots all
```

**Features:**
- Boundary layer velocity, temperature, and species profiles
- Heat flux analysis (conductive, diffusive, total)
- Wall properties and gradients
- Flow field contour plots
- Convergence history
- Summary statistics export (JSON)

**Usage options:**
```bash
# Generate all plots and summary
python3 postprocess_blast.py --input simulation.h5 --plots all

# Specific station analysis
python3 postprocess_blast.py --input simulation.h5 --plots profiles --station 0

# Heat flux maps only
python3 postprocess_blast.py --input simulation.h5 --plots f_g_map
```

### Temperature Analysis Script

```bash
python3 postprocess_blast_temperatures.py --input simulation.h5 --plots all
```

**Specialized for:**
- Detailed temperature profile analysis
- Thermal boundary layer characteristics
- Temperature gradient visualization
- Multi-station temperature comparison
- Wall temperature effects

### Abaque Processing Script

```bash
python3 postprocess_abaque.py --input simulation_abaque.h5
```

**Purpose:**
- Process abaque simulation results
- Generate heat flux contour maps
- Export data for thermal protection system design
- Create catalyticity effect visualizations
- Produce engineering-ready design charts

**Outputs:**
- 2D contour plots of heat flux vs. wall temperature and catalyticity
- Data tables for interpolation
- Design curves for specific operating points
- Comparison plots between different gas compositions

## Examples

### Example 1: CO2 Stagnation Point Flow with Catalytic Wall

```yaml
simulation:
  body_type: axisymmetric
  only_stagnation_point: true
  chemical_mode: "non_equilibrium"
  catalytic_wall: true
  consider_thermal_diffusion: true

wall_parameters:
  temperatures: [300, 500, 1000, 1500, 2000]

abaque:
  enabled: true
  catalyticity_values: [0.0, 0.01, 0.1, 1.0]
  temperature_range: [300, 2500]
  temperature_points: 50
```

```bash
# Run simulation
./blast config/CO2_catalytic.yaml

# Process results
python3 scripts/postprocess_abaque.py --input test_outputs/simulation_abaque.h5
```

### Example 2: High-Temperature Case with Continuation

```yaml
simulation:
  chemical_mode: "non_equilibrium"
  catalytic_wall: true

outer_edge:
  edge_points:
    - temperature: 8000.0    # Very high temperature
      pressure: 50000.0     # High pressure

continuation:
  wall_temperature_stable: 1000.0
  edge_temperature_stable: 4000.0
  pressure_stable: 5000.0
```

### Example 3: Multi-Station Analysis

```yaml
simulation:
  only_stagnation_point: false

output:
  x_stations: [0.0, 0.002, 0.004, 0.006, 0.008, 0.01]
```

```bash
# Run multi-station simulation
./blast config/multi_station.yaml

# Analyze temperature profiles
python3 scripts/postprocess_blast_temperatures.py --input test_outputs/simulation_*.h5 --plots all
```

## Output

The solver generates HDF5 output files containing:
- Complete flow field data (velocity, temperature, species concentrations)
- Boundary layer profiles at specified stations
- Heat flux components (conductive, diffusive, catalytic)
- Thermodynamic properties and transport coefficients
- Convergence metrics and simulation metadata

Output files are saved as `simulation_YYYYMMDD_HHMMSS.h5` in the configured output directory.

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