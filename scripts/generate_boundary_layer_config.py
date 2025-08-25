#!/usr/bin/env python3
"""
Generate YAML configuration files for boundary layer simulations with multiple edge points.
Based on CO2_5_reconstruction.yaml structure.
"""

import yaml
import numpy as np
from pathlib import Path

def generate_config(
    num_points=5,
    x_range=(0.0, 1),  # Start from 0.001 to avoid x=0 issue
    temperature_range=(5905.5, 3000),
    pressure_range=(7000, 7000),
    velocity_range=(0, 100),  # Non-zero velocities
    radius_range=(1, 1),
    wall_temperature_range=(3000, 3000),
    output_file="boundary_layer_config.yaml"
):
    """
    Generate a YAML configuration file for boundary layer simulation.
    
    Parameters:
    -----------
    num_points : int
        Number of edge points to generate
    x_range : tuple
        (min, max) for x coordinates
    temperature_range : tuple
        (min, max) for temperature values
    pressure_range : tuple
        (min, max) for pressure values
    velocity_range : tuple
        (min, max) for velocity values
    radius_range : tuple
        (min, max) for radius values
    wall_temperature_range : tuple
        (min, max) for wall temperature values
    output_file : str
        Name of the output YAML file
    """
    
    # Generate linearly spaced values
    x_values = np.linspace(x_range[0], x_range[1], num_points)
    temperatures = np.linspace(temperature_range[0], temperature_range[1], num_points)
    pressures = np.linspace(pressure_range[0], pressure_range[1], num_points)
    velocities = np.linspace(velocity_range[0], velocity_range[1], num_points)
    radii = np.linspace(radius_range[0], radius_range[1], num_points)
    wall_temperatures = np.linspace(wall_temperature_range[0], wall_temperature_range[1], num_points)
    
    # Create edge points without boundary_override
    edge_points = []
    for i in range(num_points):
        point = {
            'x': float(x_values[i]),
            'radius': float(radii[i]),
            'velocity': float(velocities[i]),
            'temperature': float(temperatures[i]),
            'pressure': float(pressures[i])
        }
        edge_points.append(point)
    
    # Build the complete configuration following CO2_5_reconstruction structure
    config = {}
    
    # Add header comments as plain strings (will be handled separately)
    config['_comment1'] = "# ============================================================================="
    config['_comment2'] = "# CONTINUATION METHOD PARAMETERS FOR CONVERGENCE"
    config['_comment3'] = "# ============================================================================="
    
    # Continuation section - keep fixed values
    config['continuation'] = {
        'wall_temperature_stable': 3000,
        'edge_temperature_stable': 5905.5,
        'pressure_stable': 7000
    }
    
    config['_comment4'] = "# ============================================================================="
    config['_comment5'] = "# BASE MODE"
    config['_comment6'] = "# ============================================================================="
    
    # Base section
    config['base'] = {
        'enabled': True,
        
        'simulation': {
            'body_type': 'axisymmetric',
            'only_stagnation_point': False,  # Changed to False for full boundary layer
            'finite_thickness': False,
            'diffusion_type': 'stefan_maxwell',
            'consider_thermal_diffusion': False,
            'consider_dufour_effect': False,
            'chemical_mode': 'non_equilibrium',
            'catalytic_wall': False,
            'wall_mode': 'imposed_temperature'
        },
        
        'numerical': {
            'n_eta': 20,
            'eta_max': 6.0,
            'convergence_tolerance': 1.0e-9,
            'max_iterations': 100000
        },
        
        'mixture': {
            'name': 'CO2_5',
            'thermodynamic_database': 'RRHO',
            'viscosity_algorithm': 'chapman_enskog_ldlt',
            'thermal_conductivity_algorithm': 'chapman_enskog_ldlt',
            'state_model': 'ChemNonEq1T'
        },
        
        'output': {
            'x_stations': x_values.tolist(),
            'output_directory': 'test_outputs'
        },
        
        'outer_edge': {
            'edge_points': edge_points,
            
            'velocity_gradient_stagnation': 4027.1213622326,
            'freestream_density': 0.0019382429,
            'freestream_velocity': 500000,
            
            'finite_thickness_params': {
                'v_edge': -111.31796000,
                'd2_ue_dxdy': -290979.9698604651,
                'delta_bl': 0.0105250000
            }
        },
        
        'wall_parameters': {
            'temperatures': wall_temperatures.tolist(),
            'emissivity': 0,
            'environment_temperature': 300
        }
    }
    
    # Write to YAML file with proper formatting
    output_path = Path(output_file)
    
    # Custom YAML representation to handle comments
    yaml_content = []
    yaml_content.append("# =============================================================================")
    yaml_content.append("# CONTINUATION METHOD PARAMETERS FOR CONVERGENCE")
    yaml_content.append("# =============================================================================")
    yaml_content.append("")
    
    # Continuation section
    yaml_content.append("continuation:")
    yaml_content.append(f"  wall_temperature_stable: {config['continuation']['wall_temperature_stable']}")
    yaml_content.append(f"  edge_temperature_stable: {config['continuation']['edge_temperature_stable']}")
    yaml_content.append(f"  pressure_stable: {config['continuation']['pressure_stable']}")
    yaml_content.append("")
    
    yaml_content.append("# =============================================================================")
    yaml_content.append("# BASE MODE")
    yaml_content.append("# =============================================================================")
    yaml_content.append("")
    
    # Convert base section to YAML
    base_yaml = yaml.dump({'base': config['base']}, default_flow_style=False, sort_keys=False)
    yaml_content.append(base_yaml)
    
    # Write complete file
    with open(output_path, 'w') as f:
        f.write('\n'.join(yaml_content))
    
    print(f"Configuration file generated: {output_path}")
    print(f"Number of edge points: {num_points}")
    print(f"X range: {x_range}")
    print(f"Temperature range: {temperature_range} K")
    print(f"Pressure range: {pressure_range} Pa")
    print(f"Velocity range: {velocity_range} m/s")
    

if __name__ == "__main__":
    # Use default parameters - modify if needed
    generate_config(
        output_file="config/boundary_layer_auto.yaml"
    )