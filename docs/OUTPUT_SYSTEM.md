# BLAST Output System Documentation

## Overview

The BLAST boundary layer solver features a comprehensive, professional-grade output system designed for modern CFD post-processing workflows. The system supports multiple file formats, metadata preservation, and seamless integration with popular analysis tools.

## Supported Formats

### HDF5 (Recommended)
- **Extension**: `.h5`, `.hdf5`
- **Features**: Hierarchical data organization, compression, metadata, cross-platform
- **Use Cases**: Primary data storage, large simulations, programmatic analysis
- **Tools**: Python (h5py, pandas), MATLAB, HDFView, VisIt

### VTK XML
- **Extension**: `.vts`, `.vtu`
- **Features**: ParaView native format, binary compression, time series
- **Use Cases**: 3D visualization, streamlines, vector fields
- **Tools**: ParaView, VisIt, Mayavi

### CSV
- **Extension**: `.csv`
- **Features**: Simple text format, Excel compatibility, human-readable
- **Use Cases**: Quick analysis, spreadsheet import, plotting
- **Tools**: Excel, Python (pandas), MATLAB, Origin

## Quick Start

### Basic Usage
```cpp
#include "blast/io/output/output_writer.hpp"

// Configure output
blast::io::output::OutputConfig config;
config.primary_format = blast::io::output::OutputFormat::HDF5;
config.additional_formats = {blast::io::output::OutputFormat::VTK_XML};
config.base_directory = "simulation_results";

// Create writer
blast::io::output::OutputWriter writer(config);

// Write solution
auto result = writer.write_solution(solution, sim_config, mixture, "hypersonic_flow");
```

### Command Line
```bash
# Run simulation with automatic output
./blast config.yaml simulation_name

# Post-process results
python3 scripts/postprocess_blast.py --input simulation_name.h5 --plots all
```

## Configuration Options

### OutputConfig Structure
```cpp
struct OutputConfig {
    std::filesystem::path base_directory = "BLAST_outputs";
    std::string case_name = "simulation";
    OutputFormat primary_format = OutputFormat::HDF5;
    std::vector<OutputFormat> additional_formats;
    
    bool save_metadata = true;
    bool save_derivatives = true;
    bool save_coefficients = false;
    bool compress_data = true;
    bool include_timestamp = true;
    
    struct VariableSelection {
        bool flow_variables = true;
        bool species_concentrations = true;
        bool transport_properties = false;
        bool chemical_rates = false;
        bool diffusion_fluxes = false;
    } variables;
};
```

### Variable Selection Guide
- **flow_variables**: F, g, V, Temperature, Pressure, Density
- **species_concentrations**: Mass fractions for all species
- **transport_properties**: Viscosity, thermal conductivity, diffusion coefficients
- **chemical_rates**: Production rates, reaction rates
- **diffusion_fluxes**: Stefan-Maxwell fluxes

## File Organization

### HDF5 Structure
```
simulation.h5
├── metadata/
│   ├── simulation_config
│   ├── grid_info
│   ├── mixture_properties
│   └── convergence_info
├── stations/
│   ├── station_000/
│   │   ├── coordinates (xi, x, eta, y)
│   │   ├── flow_variables (F, g, V, T, P, rho)
│   │   ├── species_concentrations
│   │   ├── transport_properties (optional)
│   │   └── derivatives/ (optional)
│   ├── station_001/
│   └── ...
├── wall/
│   ├── x_positions
│   ├── temperatures
│   ├── heat_flux
│   └── shear_stress
└── integrated/
    ├── displacement_thickness
    ├── momentum_thickness
    ├── shape_factor
    └── total_enthalpy_thickness
```

### CSV Structure
```
simulation_results/
├── metadata.txt
├── station_000.csv
├── station_001.csv
├── ...
├── wall_properties.csv
└── integrated_quantities.csv
```

## Advanced Usage

### Streaming Output for Large Simulations
```cpp
// Create stream writer for memory-efficient output
auto stream_writer = writer.create_stream_writer(config, mixture, "large_simulation");

// Add stations progressively
for (int station = 0; station < n_stations; ++station) {
    // Solve station...
    stream_writer->add_station(solution_state, xi, x_physical, temperature);
}

// Finalize output
auto files = stream_writer->finalize();
```

### Custom Format Configuration
```cpp
// HDF5 with maximum compression
blast::io::output::HDF5Config hdf5_config;
hdf5_config.compression_level = 9;
hdf5_config.use_shuffle_filter = true;
hdf5_config.use_fletcher32 = true;

auto hdf5_writer = std::make_unique<blast::io::output::HDF5Writer>(hdf5_config);

// VTK with ASCII output for debugging
blast::io::output::VTKConfig vtk_config;
vtk_config.data_mode = blast::io::output::VTKConfig::DataMode::ASCII;
vtk_config.precision = 12;

auto vtk_writer = std::make_unique<blast::io::output::VTKXMLWriter>(vtk_config);
```

### Progress Monitoring
```cpp
auto progress_callback = [](double progress, const std::string& stage) {
    std::cout << "\r" << stage << " [" << std::setw(6) << std::fixed 
              << std::setprecision(1) << (progress * 100.0) << "%]" << std::flush;
};

writer.write_solution(solution, config, mixture, "simulation", progress_callback);
```

## Post-Processing

### Python Analysis
```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Load HDF5 data
with h5py.File('simulation.h5', 'r') as f:
    # Access metadata
    n_eta = f['metadata/grid/n_eta'][()]
    species_names = f['metadata/mixture/species_names'][:]
    
    # Access station data
    station_0 = f['stations/station_000']
    eta = station_0['eta'][:]
    F = station_0['F'][:]
    temperature = station_0['temperature'][:]
    
    # Plot boundary layer profile
    plt.figure()
    plt.plot(F, eta, label='F')
    plt.xlabel('F (Stream Function)')
    plt.ylabel('η (Similarity Coordinate)')
    plt.grid(True)
    plt.show()
```

### Automated Post-Processing
```bash
# Generate all plots
python3 scripts/postprocess_blast.py --input simulation.h5 --plots all --output results/

# Generate specific plots
python3 scripts/postprocess_blast.py --input simulation.h5 --plots profiles --station 0
python3 scripts/postprocess_blast.py --input simulation.h5 --plots wall
python3 scripts/postprocess_blast.py --input simulation.h5 --plots integrated

# Create summary report
python3 scripts/postprocess_blast.py --input simulation.h5 --plots summary --output report/
```

## Performance Considerations

### File Size Optimization
```cpp
// Minimal output for quick runs
config.variables.transport_properties = false;
config.variables.chemical_rates = false;
config.variables.diffusion_fluxes = false;
config.save_derivatives = false;

// Maximum compression
config.compress_data = true;
config.compression_level = 9.0;
```

### Memory Management
- Use streaming output for simulations with >100 stations
- Consider CSV format for simple profiles (<50 MB total)
- HDF5 chunking optimizes random access patterns

### Parallel I/O
The output system is designed for future parallel I/O extension:
```cpp
// Future: MPI-parallel output
writer.set_parallel_mode(true);
writer.set_mpi_communicator(MPI_COMM_WORLD);
```

## Integration Examples

### Research Workflow
```cpp
// High-fidelity research simulation
OutputConfig research_config;
research_config.primary_format = OutputFormat::HDF5;
research_config.additional_formats = {OutputFormat::VTK_XML, OutputFormat::CSV};
research_config.save_derivatives = true;
research_config.save_coefficients = true;
research_config.variables.transport_properties = true;
research_config.variables.chemical_rates = true;
research_config.variables.diffusion_fluxes = true;
```

### Engineering Analysis
```cpp
// Fast engineering analysis
OutputConfig engineering_config;
engineering_config.primary_format = OutputFormat::CSV;
engineering_config.additional_formats = {OutputFormat::VTK_XML};
engineering_config.save_derivatives = false;
engineering_config.variables.transport_properties = false;
engineering_config.variables.chemical_rates = false;
```

### Visualization Focus
```cpp
// Visualization-optimized output
OutputConfig viz_config;
viz_config.primary_format = OutputFormat::VTK_XML;
viz_config.save_metadata = false;
viz_config.variables.flow_variables = true;
viz_config.variables.species_concentrations = true;
```

## Error Handling

### Common Issues and Solutions

**Error**: `HDF5 library not found`
```bash
# Ubuntu/Debian
sudo apt-get install libhdf5-dev

# macOS
brew install hdf5

# Check installation
pkg-config --exists hdf5
```

**Error**: `Cannot create output directory`
```cpp
// Ensure directory exists and is writable
std::filesystem::create_directories(config.base_directory);
auto perms = std::filesystem::status(config.base_directory).permissions();
```

**Error**: `File size too large`
```cpp
// Reduce output size
config.save_derivatives = false;
config.variables.diffusion_fluxes = false;
config.compression_level = 9.0;
```

### Validation
```cpp
// Validate output configuration
auto validation = writer.validate_config();
if (!validation) {
    std::cerr << "Config error: " << validation.error().message() << std::endl;
}

// Check output files
for (const auto& file_path : output_files) {
    if (file_path.extension() == ".h5") {
        auto validate_result = blast::io::output::hdf5::validate_file(file_path);
        if (!validate_result) {
            std::cerr << "HDF5 validation failed: " << validate_result.error().message() << std::endl;
        }
    }
}
```

## Best Practices

### File Naming
- Use descriptive case names: `mach6_stagnation`, `reentry_trajectory`
- Include timestamps for parametric studies
- Avoid spaces and special characters

### Data Organization
- Group related simulations in subdirectories
- Use consistent variable selection within studies
- Document simulation parameters in metadata

### Version Control
- Store configuration files in version control
- Exclude large output files (`.h5`, `.vts`) from git
- Use `.gitignore` patterns:
```gitignore
# BLAST outputs
*.h5
*.hdf5
*.vts
*.vtk
*_outputs/
BLAST_outputs/
```

### Archival Storage
- Compress HDF5 files for long-term storage
- Include README files with simulation descriptions
- Consider converting large datasets to read-only formats

## Troubleshooting

### Performance Issues
1. **Slow writing**: Enable compression, reduce variable count
2. **Large files**: Use streaming output, disable derivatives
3. **Memory usage**: Process stations individually

### Compatibility Issues
1. **HDF5 version**: Ensure compatible versions across tools
2. **VTK format**: Use XML format for best compatibility
3. **Endianness**: HDF5 handles automatically, VTK uses little-endian

### Data Integrity
1. **Validation**: Always validate output files
2. **Backup**: Keep multiple format copies for important runs
3. **Checksums**: Use HDF5 Fletcher32 filter for error detection

## API Reference

### Key Classes
- `OutputWriter`: Main interface for writing simulation data
- `OutputConfig`: Configuration structure for output options
- `OutputDataset`: Complete dataset structure
- `HDF5Writer`: HDF5-specific implementation
- `VTKXMLWriter`: VTK XML implementation
- `CSVWriter`: CSV implementation

### Key Functions
- `write_solution()`: Write complete solution
- `create_stream_writer()`: Create streaming writer
- `validate_config()`: Validate configuration
- `get_output_info()`: Get planned output information

For detailed API documentation, see the header files in `include/blast/io/output/`.

## Support

For issues, feature requests, or questions:
1. Check this documentation
2. Review example configurations
3. Run validation checks
4. Contact development team

The output system is designed to be extensible - new formats can be added by implementing the `FormatWriter` interface.