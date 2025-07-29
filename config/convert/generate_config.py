import numpy as np
import yaml

def generate_config():
    N = 51
    x_min, x_max = 0.0, 1.0
    v_min, v_max = 0.0, 100.0
    temp_min, temp_max = 3000.0, 5000.0
    pressure_base = 7000.0
    output_dir = "test_outputs"
    output_file = "config_new_standard.yaml"

    x_vals = np.linspace(x_min, x_max, N).round(4).tolist()
    v_vals = np.linspace(v_min, v_max, N).round(4).tolist()
    temp_vals = np.linspace(temp_min, temp_max, N).round(1).tolist()
    config = {
        "simulation": {
            "body_type": "axisymmetric",
            "only_stagnation_point": False,
            "diffusion_type": "stefan_maxwell",
            "consider_thermal_diffusion": False,
            "chemical_non_equilibrium": False,
        },
        "numerical": {
            "n_eta": 50,
            "eta_max": 6.0,
            "convergence_tolerance": 1.0e-6,
            "max_iterations": 100000,
            "under_relaxation": 1,
            "step_control": {"lower_bound": 5, "upper_bound": 20},
            "solvers": {
                "h2t_tolerance": 1.0e-12,
                "h2t_max_iterations": 50000,
                "stefan_tolerance": 1.0e-12,
                "stefan_max_iterations": 50000,
            },
        },
        "mixture": {
            "name": "CO2_5",
            "thermodynamic_database": "RRHO",
            "viscosity_algorithm": "chapman_enskog_ldlt",
            "reference_temperature": 0.0,
        },
        "output": {
            "x_stations": x_vals,
            "output_directory": output_dir,
            "generate_lookup_table": False,
        },
        "outer_edge": {
            "edge_points": [
                {
                    "x": x,
                    "radius": 1,
                    "velocity": v,
                    "temperature": temp,
                    "pressure": pressure_base,
                }
                for x, v, temp in zip(x_vals, v_vals, temp_vals)
            ],
            "velocity_gradient_stagnation": 4027.1213622326,
            "freestream_density": 0.0019382429,
            "freestream_velocity": 500000,
        },
        "wall_parameters": {
            "temperatures": [3100]
        },
    }

    with open(output_file, 'w') as f:
        yaml.safe_dump(config, f, sort_keys=False, default_flow_style=False)

    print(f"Configuration generated: {output_file}")
    print(f"Number of points generated: {N}")
    print(f"X range: [{x_min}, {x_max}]")
    print(f"Velocity range: [{v_min}, {v_max}]")
    print(f"Temperature range: [{temp_min}, {temp_max}]")

if __name__ == "__main__":
    generate_config()