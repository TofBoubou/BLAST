import numpy as np
import yaml

def generate_config():
    # Paramètres embarqués
    N = 501
    x_min, x_max = 0.0, 1.0
    v_min, v_max = 0.0, 100.0
    t_min, t_max = 500.0, 500.0
    output_dir = "test_outputs"
    output_file = "config_50points.yaml"

    # Génération des listes
    x_vals = np.linspace(x_min, x_max, N).round(4).tolist()
    v_vals = np.linspace(v_min, v_max, N).round(4).tolist()
    temps = np.linspace(t_min, t_max, N).round(1).tolist()

    # Construction de la configuration YAML
    config = {
        "simulation": {
            "body_type": "axisymmetric",
            "only_stagnation_point": False,
            "diffusion_type": "stefan_maxwell",
            "consider_thermal_diffusion": True,
            "chemical_non_equilibrium": False,
        },
        "numerical": {
            "n_eta": 20,
            "eta_max": 4.0,
            "convergence_tolerance": 1.0e-6,
            "max_iterations": 1000000,
            "under_relaxation": 0.01,
            "step_control": {"lower_bound": 5, "upper_bound": 20},
            "solvers": {
                "h2t_tolerance": 1.0e-12,
                "h2t_max_iterations": 50000,
                "stefan_tolerance": 1.0e-12,
                "stefan_max_iterations": 50000,
            },
        },
        "mixture": {
            "name": "air_5",
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
                    "radius": 1.0,
                    "velocity": v,
                    "enthalpy": 800000.0,
                    "pressure": 1000.0,
                    "density": 0.01,
                    "viscosity": 2.0e-5,
                    "species": [0.0, 0.0, 0.0, 0.767, 0.233],
                }
                for x, v in zip(x_vals, v_vals)
            ],
            "velocity_gradient_stagnation": 100000.0,
            "freestream_density": 0.01,
            "freestream_velocity": 3000.0,
        },
        "wall_parameters": {
            "temperatures": temps
        },
    }

    # Sauvegarde du fichier YAML
    with open(output_file, 'w') as f:
        yaml.safe_dump(config, f, sort_keys=False)

    print(f"Config généré : {output_file}")

if __name__ == "__main__":
    generate_config()