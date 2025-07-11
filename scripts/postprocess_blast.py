#!/usr/bin/env python3
"""
BLAST Boundary Layer Post-Processing Script

This script provides comprehensive post-processing capabilities for BLAST simulation results.
It can read HDF5, VTK, and CSV output formats and generate publication-quality plots.

Usage:
    python postprocess_blast.py --input simulation.h5 --plots all
    python postprocess_blast.py --input simulation.h5 --plots profiles --station 0
    python postprocess_blast.py --input csv_dir/ --format csv --plots wall
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import vtk
from vtk.util import numpy_support
import seaborn as sns
from pathlib import Path
import json
from typing import Dict, List, Optional, Union, Tuple
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality plot style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")


class BLASTReader:
    """Unified reader for BLAST output formats"""
    
    def __init__(self, input_path: Union[str, Path]):
        self.input_path = Path(input_path)
        self.format = self._detect_format()
        self.data = None
        self.metadata = None
        
    def _detect_format(self) -> str:
        """Auto-detect input format"""
        if self.input_path.suffix.lower() in ['.h5', '.hdf5']:
            return 'hdf5'
        elif self.input_path.suffix.lower() in ['.vts', '.vtk']:
            return 'vtk'
        elif self.input_path.is_dir():
            # Check for CSV files
            csv_files = list(self.input_path.glob('*.csv'))
            if csv_files:
                return 'csv'
        raise ValueError(f"Cannot detect format for: {self.input_path}")
    
    def load_data(self) -> Dict:
        """Load data based on detected format"""
        if self.format == 'hdf5':
            return self._load_hdf5()
        elif self.format == 'csv':
            return self._load_csv()
        elif self.format == 'vtk':
            return self._load_vtk()
        else:
            raise NotImplementedError(f"Format {self.format} not implemented")
    
    def _load_hdf5(self) -> Dict:
        """Load HDF5 data with comprehensive error handling"""
        try:
            with h5py.File(self.input_path, 'r') as f:
                data = {}
                
                # Load metadata
                if 'metadata' in f:
                    data['metadata'] = self._read_hdf5_group(f['metadata'])
                
                # Load stations
                if 'stations' in f:
                    data['stations'] = {}
                    for station_name in f['stations'].keys():
                        station_data = self._read_hdf5_group(f['stations'][station_name])
                        data['stations'][station_name] = station_data
                
                # Load wall data
                if 'wall' in f:
                    data['wall'] = self._read_hdf5_group(f['wall'])
                
                # Load integrated quantities
                if 'integrated' in f:
                    data['integrated'] = self._read_hdf5_group(f['integrated'])
                
                self.data = data
                self.metadata = data.get('metadata', {})
                return data
                
        except Exception as e:
            raise RuntimeError(f"Failed to load HDF5 file: {e}")
    
    def _read_hdf5_group(self, group) -> Dict:
        """Recursively read HDF5 group"""
        result = {}
        for key, item in group.items():
            if isinstance(item, h5py.Group):
                result[key] = self._read_hdf5_group(item)
            elif isinstance(item, h5py.Dataset):
                result[key] = np.array(item)
                # Read attributes
                if item.attrs:
                    result[f'{key}_attrs'] = dict(item.attrs)
        return result
    
    def _load_csv(self) -> Dict:
        """Load CSV directory structure"""
        data = {'stations': {}}
        
        # Load metadata if available
        metadata_file = self.input_path / 'metadata.txt'
        if metadata_file.exists():
            data['metadata'] = self._parse_metadata_file(metadata_file)
        
        # Load station files
        station_files = sorted(self.input_path.glob('station_*.csv'))
        for i, station_file in enumerate(station_files):
            station_data = pd.read_csv(station_file, comment='#')
            data['stations'][f'station_{i:03d}'] = station_data.to_dict('series')
        
        # Load wall properties
        wall_file = self.input_path / 'wall_properties.csv'
        if wall_file.exists():
            data['wall'] = pd.read_csv(wall_file, comment='#').to_dict('series')
        
        # Load integrated quantities
        integrated_file = self.input_path / 'integrated_quantities.csv'
        if integrated_file.exists():
            data['integrated'] = pd.read_csv(integrated_file, comment='#').to_dict('series')
        
        self.data = data
        return data
    
    def _parse_metadata_file(self, file_path: Path) -> Dict:
        """Parse BLAST metadata file"""
        metadata = {}
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Try to convert to appropriate type
                    try:
                        if value.lower() in ['true', 'false']:
                            metadata[key] = value.lower() == 'true'
                        elif ',' in value:
                            metadata[key] = [v.strip() for v in value.split(',')]
                        else:
                            metadata[key] = float(value) if '.' in value else int(value)
                    except ValueError:
                        metadata[key] = value
        return metadata
    
    def _load_vtk(self) -> Dict:
        """Load VTK data (basic implementation)"""
        # This would require more complex VTK parsing
        raise NotImplementedError("VTK loading not fully implemented")
    
    def get_station_data(self, station_index: int) -> Dict:
        """Get data for specific station"""
        if self.data is None:
            self.load_data()
        
        station_key = f'station_{station_index:03d}'
        if station_key in self.data['stations']:
            return self.data['stations'][station_key]
        else:
            available = list(self.data['stations'].keys())
            raise KeyError(f"Station {station_index} not found. Available: {available}")
    
    def get_species_names(self) -> List[str]:
        """Extract species names from metadata or data"""
        if self.metadata and 'mixture' in self.metadata:
            if 'species_names' in self.metadata['mixture']:
                return list(self.metadata['mixture']['species_names'])
        
        # Fallback: extract from column names
        if self.data and 'stations' in self.data:
            station_data = next(iter(self.data['stations'].values()))
            species_cols = [col for col in station_data.keys() if col.startswith('c_')]
            return [col[2:] for col in species_cols]  # Remove 'c_' prefix
        
        return []


class BLASTPlotter:
    """Publication-quality plotting for BLAST results"""
    
    def __init__(self, reader: BLASTReader):
        self.reader = reader
        self.fig_size = (12, 8)
        self.dpi = 300
        
    def plot_profiles(self, station_index: int = 0, save_path: Optional[Path] = None) -> None:
        """Plot boundary layer profiles for a given station"""
        
        station_data = self.reader.get_station_data(station_index)
        species_names = self.reader.get_species_names()
        
        # Create subplots
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'Boundary Layer Profiles - Station {station_index}', fontsize=16, fontweight='bold')
        
        # Get eta coordinates
        eta = np.array(station_data['eta'])
        
        # Plot 1: Velocity profile (F)
        axes[0, 0].plot(station_data['F'], eta, 'b-', linewidth=2, label='F')
        axes[0, 0].set_xlabel('F (Dimensionless Stream Function)')
        axes[0, 0].set_ylabel('η (Similarity Coordinate)')
        axes[0, 0].grid(True, alpha=0.3)
        axes[0, 0].set_title('Velocity Profile')
        
        # Plot 2: Temperature profile (g)
        axes[0, 1].plot(station_data['g'], eta, 'r-', linewidth=2, label='g')
        axes[0, 1].set_xlabel('g (Dimensionless Enthalpy)')
        axes[0, 1].set_ylabel('η')
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_title('Enthalpy Profile')
        
        # Plot 3: Physical temperature
        if 'Temperature' in station_data:
            axes[0, 2].plot(station_data['Temperature'], eta, 'orange', linewidth=2)
            axes[0, 2].set_xlabel('Temperature (K)')
            axes[0, 2].set_ylabel('η')
            axes[0, 2].grid(True, alpha=0.3)
            axes[0, 2].set_title('Temperature Profile')
        
        # Plot 4: Major species concentrations
        colors = plt.cm.tab10(np.linspace(0, 1, len(species_names)))
        for i, (species, color) in enumerate(zip(species_names[:6], colors)):
            col_name = f'c_{species}' if f'c_{species}' in station_data else species
            if col_name in station_data:
                axes[1, 0].plot(station_data[col_name], eta, color=color, 
                               linewidth=2, label=species)
        
        axes[1, 0].set_xlabel('Mass Fraction')
        axes[1, 0].set_ylabel('η')
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        axes[1, 0].set_title('Species Concentrations')
        
        # Plot 5: Density and pressure
        if 'Density' in station_data:
            ax_rho = axes[1, 1]
            ax_rho.plot(station_data['Density'], eta, 'g-', linewidth=2, label='Density')
            ax_rho.set_xlabel('Density (kg/m³)')
            ax_rho.set_ylabel('η')
            ax_rho.grid(True, alpha=0.3)
            ax_rho.set_title('Density Profile')
            
            if 'Pressure' in station_data:
                ax_p = ax_rho.twinx()
                ax_p.plot(station_data['Pressure'], eta, 'm--', linewidth=2, label='Pressure')
                ax_p.set_xlabel('Pressure (Pa)')
                ax_p.legend(loc='upper right')
        
        # Plot 6: V profile
        if 'V' in station_data:
            axes[1, 2].plot(station_data['V'], eta, 'purple', linewidth=2)
            axes[1, 2].set_xlabel('V (Velocity Component)')
            axes[1, 2].set_ylabel('η')
            axes[1, 2].grid(True, alpha=0.3)
            axes[1, 2].set_title('V Profile')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
            print(f"Profiles saved to: {save_path}")
        else:
            plt.show()
    
    def plot_wall_properties(self, save_path: Optional[Path] = None) -> None:
        """Plot wall properties along the body"""
        
        if 'wall' not in self.reader.data:
            print("No wall data available")
            return
        
        wall_data = self.reader.data['wall']
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Wall Properties', fontsize=16, fontweight='bold')
        
        x_pos = np.array(wall_data['x_positions'])
        
        # Wall temperature
        if 'temperatures' in wall_data:
            axes[0, 0].plot(x_pos, wall_data['temperatures'], 'r-', linewidth=2)
            axes[0, 0].set_xlabel('x (m)')
            axes[0, 0].set_ylabel('Temperature (K)')
            axes[0, 0].set_title('Wall Temperature')
            axes[0, 0].grid(True, alpha=0.3)
        
        # Heat flux
        if 'heat_flux' in wall_data:
            axes[0, 1].plot(x_pos, wall_data['heat_flux'], 'b-', linewidth=2)
            axes[0, 1].set_xlabel('x (m)')
            axes[0, 1].set_ylabel('Heat Flux (W/m²)')
            axes[0, 1].set_title('Wall Heat Flux')
            axes[0, 1].grid(True, alpha=0.3)
        
        # Shear stress
        if 'shear_stress' in wall_data:
            axes[1, 0].plot(x_pos, wall_data['shear_stress'], 'g-', linewidth=2)
            axes[1, 0].set_xlabel('x (m)')
            axes[1, 0].set_ylabel('Shear Stress (Pa)')
            axes[1, 0].set_title('Wall Shear Stress')
            axes[1, 0].grid(True, alpha=0.3)
        
        # Combined plot
        if all(key in wall_data for key in ['temperatures', 'heat_flux']):
            ax1 = axes[1, 1]
            color = 'tab:red'
            ax1.set_xlabel('x (m)')
            ax1.set_ylabel('Temperature (K)', color=color)
            ax1.plot(x_pos, wall_data['temperatures'], color=color, linewidth=2)
            ax1.tick_params(axis='y', labelcolor=color)
            
            ax2 = ax1.twinx()
            color = 'tab:blue'
            ax2.set_ylabel('Heat Flux (W/m²)', color=color)
            ax2.plot(x_pos, wall_data['heat_flux'], color=color, linewidth=2)
            ax2.tick_params(axis='y', labelcolor=color)
            
            axes[1, 1].set_title('Temperature & Heat Flux')
            axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
            print(f"Wall properties saved to: {save_path}")
        else:
            plt.show()
    
    def plot_integrated_quantities(self, save_path: Optional[Path] = None) -> None:
        """Plot boundary layer integral parameters"""
        
        if 'integrated' not in self.reader.data:
            print("No integrated quantities data available")
            return
        
        integrated_data = self.reader.data['integrated']
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Boundary Layer Integral Parameters', fontsize=16, fontweight='bold')
        
        # Get xi coordinates
        if 'metadata' in self.reader.data and 'grid' in self.reader.data['metadata']:
            xi = np.array(self.reader.data['metadata']['grid']['xi_coordinates'])
        else:
            # Create dummy xi coordinates
            n_points = len(integrated_data['displacement_thickness'])
            xi = np.linspace(0, 1, n_points)
        
        # Displacement thickness
        if 'displacement_thickness' in integrated_data:
            axes[0, 0].plot(xi, integrated_data['displacement_thickness'], 'b-', linewidth=2)
            axes[0, 0].set_xlabel('ξ')
            axes[0, 0].set_ylabel('δ* (m)')
            axes[0, 0].set_title('Displacement Thickness')
            axes[0, 0].grid(True, alpha=0.3)
        
        # Momentum thickness
        if 'momentum_thickness' in integrated_data:
            axes[0, 1].plot(xi, integrated_data['momentum_thickness'], 'r-', linewidth=2)
            axes[0, 1].set_xlabel('ξ')
            axes[0, 1].set_ylabel('θ (m)')
            axes[0, 1].set_title('Momentum Thickness')
            axes[0, 1].grid(True, alpha=0.3)
        
        # Shape factor
        if 'shape_factor' in integrated_data:
            axes[1, 0].plot(xi, integrated_data['shape_factor'], 'g-', linewidth=2)
            axes[1, 0].set_xlabel('ξ')
            axes[1, 0].set_ylabel('H = δ*/θ')
            axes[1, 0].set_title('Shape Factor')
            axes[1, 0].grid(True, alpha=0.3)
        
        # Enthalpy thickness
        if 'total_enthalpy_thickness' in integrated_data:
            axes[1, 1].plot(xi, integrated_data['total_enthalpy_thickness'], 'purple', linewidth=2)
            axes[1, 1].set_xlabel('ξ')
            axes[1, 1].set_ylabel('δ_h (m)')
            axes[1, 1].set_title('Enthalpy Thickness')
            axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
            print(f"Integrated quantities saved to: {save_path}")
        else:
            plt.show()
    
    def plot_convergence_history(self, save_path: Optional[Path] = None) -> None:
        """Plot convergence history if available"""
        
        # This would require convergence data to be saved during simulation
        print("Convergence history plotting not implemented (requires convergence data)")
    
    def create_summary_report(self, output_dir: Path) -> None:
        """Create a comprehensive summary report"""
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        print(f"Creating summary report in: {output_dir}")
        
        # Plot profiles for first few stations
        n_stations = len(self.reader.data['stations'])
        stations_to_plot = min(3, n_stations)
        
        for i in range(stations_to_plot):
            self.plot_profiles(i, save_path=output_dir / f'profiles_station_{i:03d}.png')
        
        # Plot wall properties
        self.plot_wall_properties(save_path=output_dir / 'wall_properties.png')
        
        # Plot integrated quantities
        self.plot_integrated_quantities(save_path=output_dir / 'integrated_quantities.png')
        
        # Create summary data file
        self._create_summary_json(output_dir / 'summary.json')
        
        print("Summary report completed!")
    
    def _create_summary_json(self, file_path: Path) -> None:
        """Create JSON summary of key results"""
        
        summary = {
            'metadata': self.reader.metadata,
            'n_stations': len(self.reader.data['stations']),
            'species_names': self.reader.get_species_names()
        }
        
        # Add key results
        if 'wall' in self.reader.data:
            wall_data = self.reader.data['wall']
            if 'temperatures' in wall_data:
                summary['max_wall_temperature'] = float(np.max(wall_data['temperatures']))
                summary['min_wall_temperature'] = float(np.min(wall_data['temperatures']))
        
        # Convert numpy arrays to lists for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, dict):
                return {key: convert_numpy(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(item) for item in obj]
            return obj
        
        summary = convert_numpy(summary)
        
        with open(file_path, 'w') as f:
            json.dump(summary, f, indent=2)


def main():
    """Main function with command-line interface"""
    
    parser = argparse.ArgumentParser(description='BLAST Post-Processing Tool')
    parser.add_argument('--input', '-i', type=str, required=True,
                       help='Input file or directory (HDF5, VTK, or CSV)')
    parser.add_argument('--plots', '-p', type=str, default='all',
                       choices=['all', 'profiles', 'wall', 'integrated', 'summary'],
                       help='Type of plots to generate')
    parser.add_argument('--station', '-s', type=int, default=0,
                       help='Station index for profile plots')
    parser.add_argument('--output', '-o', type=str, default='blast_plots',
                       help='Output directory for plots')
    parser.add_argument('--format', '-f', type=str, choices=['hdf5', 'csv', 'vtk'],
                       help='Force input format (auto-detect if not specified)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for saved figures')
    parser.add_argument('--show', action='store_true',
                       help='Show plots interactively instead of saving')
    
    args = parser.parse_args()
    
    # Initialize reader
    try:
        reader = BLASTReader(args.input)
        if args.format:
            reader.format = args.format
        
        print(f"Loading {reader.format.upper()} data from: {reader.input_path}")
        reader.load_data()
        print(f"Loaded data with {len(reader.data['stations'])} stations")
        
    except Exception as e:
        print(f"Error loading data: {e}")
        return 1
    
    # Initialize plotter
    plotter = BLASTPlotter(reader)
    plotter.dpi = args.dpi
    
    # Create output directory
    output_path = Path(args.output) if not args.show else None
    if output_path:
        output_path.mkdir(exist_ok=True)
    
    # Generate plots
    try:
        if args.plots == 'all':
            if args.show:
                plotter.plot_profiles(args.station)
                plotter.plot_wall_properties()
                plotter.plot_integrated_quantities()
            else:
                plotter.create_summary_report(output_path)
                
        elif args.plots == 'profiles':
            save_path = output_path / f'profiles_station_{args.station:03d}.png' if output_path else None
            plotter.plot_profiles(args.station, save_path)
            
        elif args.plots == 'wall':
            save_path = output_path / 'wall_properties.png' if output_path else None
            plotter.plot_wall_properties(save_path)
            
        elif args.plots == 'integrated':
            save_path = output_path / 'integrated_quantities.png' if output_path else None
            plotter.plot_integrated_quantities(save_path)
            
        elif args.plots == 'summary':
            if not output_path:
                output_path = Path('blast_summary')
            plotter.create_summary_report(output_path)
        
        print("Post-processing completed successfully!")
        
    except Exception as e:
        print(f"Error during plotting: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())