#!/usr/bin/env python3
"""
BLAST Boundary Layer Post-Processing Script

This script provides comprehensive post-processing capabilities for BLAST simulation results.
It can read HDF5 output format and generate publication-quality plots with robust error handling.
Now includes dedicated heat flux analysis.

Usage:
    python postprocess_blast.py --input simulation.h5 --plots all
    python postprocess_blast.py --input simulation.h5 --plots profiles --station 0
    python postprocess_blast.py --input simulation.h5 --plots heat_flux --station 0
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import h5py
import seaborn as sns
from pathlib import Path
import json
from typing import Dict, List, Optional, Union, Tuple, Set
import warnings
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
warnings.filterwarnings('ignore')

# Set publication-quality plot style
plt.style.use('classic')

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'lines.linewidth': 1.5,
    'lines.markersize': 5,
    'axes.linewidth': 0.8,
    'grid.linewidth': 0.5,
    'grid.alpha': 0.3,
    'figure.dpi': 300,
    'savefig.dpi': 600,
})

class BLASTFileError(Exception):
    """Custom exception for BLAST file-related errors"""
    pass

class BLASTReader:
    """Unified reader for BLAST output formats with robust error handling"""
    
    # Required fields for profile plotting
    REQUIRED_PROFILE_FIELDS = {'F', 'g', 'V', 'temperature', 'eta'}
    
    # Required fields for heat flux plotting  
    REQUIRED_HEAT_FLUX_FIELDS = {'eta', 'heat_flux'}
    
    def __init__(self, input_path: Union[str, Path]):
        print("\n" + "="*60)
        print("BLAST POST-PROCESSING TOOL")
        print("="*60)
        print(f"\nInitializing HDF5 reader for: {input_path}")
        self.input_path = Path(input_path)
        self.format = self._detect_format()
        self.data = None
        self.metadata = None
        self.valid_stations = set()
        self.valid_heat_flux_stations = set()
        
    def _detect_format(self) -> str:
        """Auto-detect input format"""
        if not self.input_path.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_path}")
            
        if self.input_path.suffix.lower() in ['.h5', '.hdf5']:
            return 'hdf5'
        raise ValueError(f"Only HDF5 format is supported. Got: {self.input_path}")
    
    def _verify_file_integrity(self) -> None:
        """Test file integrity and required structure"""
        try:
            with h5py.File(self.input_path, 'r') as f:
                # Check for required top-level groups
                if 'stations' not in f:
                    raise BLASTFileError("Missing 'stations' group in HDF5 file")
                
                if len(f['stations'].keys()) == 0:
                    raise BLASTFileError("No station data found in 'stations' group")
                
                # Optional metadata check
                if 'metadata' not in f:
                    print("Warning: No 'metadata' group found - some features may be limited")
                
                print(f"File integrity check passed: {len(f['stations'])} stations found")
                
        except Exception as e:
            if isinstance(e, BLASTFileError):
                raise
            raise BLASTFileError(f"File integrity check failed: {e}")
    
    def load_data(self) -> Dict:
        """Load data with integrity verification"""
        # First verify file integrity
        self._verify_file_integrity()
        
        if self.format == 'hdf5':
            return self._load_hdf5()
        else:
            raise NotImplementedError(f"Format {self.format} not implemented")
    
    def _load_hdf5(self) -> Dict:
        """Load HDF5 data with comprehensive error handling"""
        try:
            with h5py.File(self.input_path, 'r') as f:
                data = {}
                
                # Load metadata (optional)
                if 'metadata' in f:
                    data['metadata'] = self._read_hdf5_group(f['metadata'])
                else:
                    data['metadata'] = {}
                
                # Load stations
                data['stations'] = {}
                for station_name in f['stations'].keys():
                    try:
                        station_data = self._read_hdf5_group(f['stations'][station_name])
                        data['stations'][station_name] = station_data
                        
                        station_idx = int(station_name.split('_')[-1])
                        
                        # Check if this station has required fields for profile plotting
                        if self._validate_station_data(station_data):
                            self.valid_stations.add(station_idx)
                        
                        # Check if this station has required fields for heat flux plotting
                        if self._validate_heat_flux_data(station_data):
                            self.valid_heat_flux_stations.add(station_idx)
                        
                    except Exception as e:
                        print(f"Warning: Failed to load station {station_name}: {e}")
                
                if len(data['stations']) == 0:
                    raise BLASTFileError("No valid station data could be loaded")
                
                print(f"\nLoaded {len(data['stations'])} stations:")
                print(f"   - {len(self.valid_stations)} valid for profile plotting")
                print(f"   - {len(self.valid_heat_flux_stations)} valid for heat flux plotting")
                
                self.data = data
                self.metadata = data.get('metadata', {})
                return data
                
        except Exception as e:
            if isinstance(e, BLASTFileError):
                raise
            raise RuntimeError(f"Failed to load HDF5 file: {e}")
    
    def _validate_station_data(self, station_data: Dict) -> bool:
        """Check if station has all required fields for profile plotting"""
        available_fields = set(station_data.keys())
        missing_fields = self.REQUIRED_PROFILE_FIELDS - available_fields
        
        if missing_fields:
            return False
        
        # Check if arrays have data and compatible shapes
        try:
            eta_len = len(station_data['eta'])
            for field in self.REQUIRED_PROFILE_FIELDS:
                if field != 'eta' and len(station_data[field]) != eta_len:
                    return False
            return eta_len > 0
        except:
            return False
    
    def _validate_heat_flux_data(self, station_data: Dict) -> bool:
        """Check if station has all required fields for heat flux plotting"""
        try:
            # Check for eta coordinates
            if 'eta' not in station_data or len(station_data['eta']) == 0:
                return False
            
            # Check for heat_flux group
            if 'heat_flux' not in station_data:
                return False
            
            heat_flux = station_data['heat_flux']
            
            # Check for required subgroups
            required_groups = ['dimensional', 'nondimensional']
            for group in required_groups:
                if group not in heat_flux:
                    return False
                
                # Check for required fields in each group
                required_fields = ['q_conductive', 'q_diffusive', 'q_total']
                for field in required_fields:
                    if field not in heat_flux[group]:
                        return False
                    
                    # Check array length compatibility
                    if len(heat_flux[group][field]) != len(station_data['eta']):
                        return False
            
            return True
            
        except Exception:
            return False
    
    def _read_hdf5_group(self, group) -> Dict:
        """Recursively read HDF5 group"""
        result = {}
        for key, item in group.items():
            if isinstance(item, h5py.Group):
                result[key] = self._read_hdf5_group(item)
            elif isinstance(item, h5py.Dataset):
                # Handle string datasets specially
                if h5py.check_string_dtype(item.dtype):
                    data = item[()]
                    if isinstance(data, np.ndarray):
                        result[key] = [s.decode('utf-8') if isinstance(s, bytes) else s for s in data]
                    else:
                        result[key] = data.decode('utf-8') if isinstance(data, bytes) else data
                else:
                    result[key] = np.array(item)
                
                # Read attributes
                if item.attrs:
                    result[f'{key}_attrs'] = dict(item.attrs)
        return result
    
    def get_station_data(self, station_index: int) -> Dict:
        """Get data for specific station with validation"""
        if self.data is None:
            self.load_data()
        
        station_key = f'station_{station_index:03d}'
        if station_key not in self.data['stations']:
            available = sorted([int(k.split('_')[-1]) for k in self.data['stations'].keys()])
            raise KeyError(f"Station {station_index} not found. Available stations: {available}")
        
        station_data = self.data['stations'][station_key].copy()
        
        # Validate required fields for profiles
        if not self._validate_station_data(station_data):
            missing = self.REQUIRED_PROFILE_FIELDS - set(station_data.keys())
            raise ValueError(f"Station {station_index} missing required fields: {missing}")
        
        # Process species concentrations if available
        if 'species_concentrations' in station_data:
            species_names = self.get_species_names()
            conc_matrix = station_data['species_concentrations']
            
            for i, species in enumerate(species_names):
                if i < conc_matrix.shape[0]:
                    station_data[f'c_{species}'] = conc_matrix[i, :]
        
        return station_data
    
    def get_station_heat_flux_data(self, station_index: int) -> Dict:
        """Get heat flux data for specific station with validation"""
        if self.data is None:
            self.load_data()
        
        station_key = f'station_{station_index:03d}'
        if station_key not in self.data['stations']:
            available = sorted([int(k.split('_')[-1]) for k in self.data['stations'].keys()])
            raise KeyError(f"Station {station_index} not found. Available stations: {available}")
        
        station_data = self.data['stations'][station_key].copy()
        
        # Validate required fields for heat flux
        if not self._validate_heat_flux_data(station_data):
            raise ValueError(f"Station {station_index} missing required heat flux fields")
        
        return station_data
    
    def get_valid_stations(self) -> Set[int]:
        """Get set of station indices that have valid data for profile plotting"""
        if self.data is None:
            self.load_data()
        return self.valid_stations.copy()
    
    def get_valid_heat_flux_stations(self) -> Set[int]:
        """Get set of station indices that have valid data for heat flux plotting"""
        if self.data is None:
            self.load_data()
        return self.valid_heat_flux_stations.copy()
    
    def get_species_names(self) -> List[str]:
        """Extract species names from metadata or data"""
        if self.metadata and 'mixture' in self.metadata:
            if 'species_names' in self.metadata['mixture']:
                names = self.metadata['mixture']['species_names']
                if isinstance(names, (list, np.ndarray)):
                    return [str(name) for name in names]
                else:
                    return [str(names)]
        
        # Fallback: extract from column names
        if self.data and 'stations' in self.data:
            station_data = next(iter(self.data['stations'].values()))
            species_cols = [col for col in station_data.keys() if col.startswith('c_')]
            return [col[2:] for col in species_cols]
        
        return []


class BLASTPlotter:
    """Publication-quality plotting for BLAST results with robust error handling"""
    
    def __init__(self, reader: BLASTReader):
        self.reader = reader
        self.fig_size = (12, 8)
        self.dpi = 300
        
    def plot_profiles(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot boundary layer profiles for a given station
        
        Returns:
            bool: True if plot was created successfully, False otherwise
        """
        
        try:
            # Check if station is valid for plotting
            valid_stations = self.reader.get_valid_stations()
            if station_index not in valid_stations:
                available = sorted(valid_stations)
                print(f"Error: Station {station_index} does not have required fields for plotting.")
                print(f"Required fields: {', '.join(self.reader.REQUIRED_PROFILE_FIELDS)}")
                print(f"Available stations with valid data: {available}")
                return False
            
            station_data = self.reader.get_station_data(station_index)
            species_names = self.reader.get_species_names()
            
        except Exception as e:
            print(f"Error loading station {station_index}: {e}")
            return False
        
        # Create subplots with clear titles
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'Boundary Layer Profiles - Station {station_index:03d}', 
                    fontsize=16, fontweight='bold')
        
        # Get eta coordinates
        eta = np.array(station_data['eta'])
        
        # Plot 1: Velocity profile (F)
        axes[0, 0].plot(station_data['F'], eta, 'black', linewidth=2, label='F')
        axes[0, 0].set_xlabel('F (Dimensionless Stream Function)')
        axes[0, 0].set_ylabel('η (Similarity Coordinate)')
        axes[0, 0].grid(True, alpha=0.3)
        axes[0, 0].set_title(f'Station {station_index:03d} - Velocity Profile (F)')
        
        # Plot 2: Temperature profile (g)
        axes[0, 1].plot(station_data['g'], eta, 'black', linewidth=2, label='g')
        axes[0, 1].set_xlabel('g (Dimensionless Enthalpy)')
        axes[0, 1].set_ylabel('η')
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_title(f'Station {station_index:03d} - Enthalpy Profile (g)')
        
        # Plot 3: Physical temperature
        axes[0, 2].plot(station_data['temperature'], eta, 'black', linewidth=2)
        axes[0, 2].set_xlabel('Temperature (K)')
        axes[0, 2].set_ylabel('η')
        axes[0, 2].grid(True, alpha=0.3)
        axes[0, 2].set_title(f'Station {station_index:03d} - Temperature Profile')
        
        # Plot 4: Major species concentrations with different line styles
        line_styles = ['-', '--', '-.', ':']  # 4 styles de base
        markers = ['', 'o', 's', '^', 'v', 'D', 'p', '*', 'h', 'H', '+', 'x']  # 12 marqueurs
        legend_added = False
        
        for i, species in enumerate(species_names):
            col_name = f'c_{species}'
            if col_name in station_data:
                # Combinaison cyclique : 4 styles × 12 marqueurs = 48 combinaisons uniques
                style_idx = i % len(line_styles)
                marker_idx = i % len(markers)
                
                axes[1, 0].plot(station_data[col_name], eta, 
                            color='black',
                            linestyle=line_styles[style_idx],
                            marker=markers[marker_idx] if markers[marker_idx] else None,
                            markevery=max(1, len(eta)//8),  # Espacer les marqueurs
                            markersize=3.5,
                            linewidth=1.5, 
                            label=species)
                legend_added = True
        
        axes[1, 0].set_xlabel('Mass Fraction')
        axes[1, 0].set_ylabel('η')
        axes[1, 0].grid(True, alpha=0.3)
        
        if legend_added:
            axes[1, 0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        axes[1, 0].set_title(f'Station {station_index:03d} - Species Concentrations')
        
        # Plot 5: Species concentrations in LOG scale
        legend_added_log = False
        min_positive_value = 1.0  # Start with max possible value
        
        for i, species in enumerate(species_names):
            col_name = f'c_{species}'
            if col_name in station_data:
                # Filter out zero/negative values for log scale
                data = station_data[col_name]
                positive_mask = data > 0
                if np.any(positive_mask):
                    # Track minimum positive value across all species
                    min_positive_value = min(min_positive_value, np.min(data[positive_mask]))
                    
                    style_idx = i % len(line_styles)
                    marker_idx = i % len(markers)
                    
                    axes[1, 1].semilogx(data[positive_mask], eta[positive_mask], 
                                color='black',
                                linestyle=line_styles[style_idx],
                                marker=markers[marker_idx] if markers[marker_idx] else None,
                                markevery=max(1, np.sum(positive_mask)//8),
                                markersize=3.5,
                                linewidth=1.5, 
                                label=species)
                    legend_added_log = True
        
        axes[1, 1].set_xlabel('Mass Fraction (log scale)')
        axes[1, 1].set_ylabel('η')
        axes[1, 1].grid(True, alpha=0.3, which='both')
        
        # Set adaptive x-axis limits with minimum bound at 1e-12
        if legend_added_log and min_positive_value < 1.0:
            # Set lower limit to one decade below the minimum value, but not below 1e-12
            lower_limit = max(min_positive_value * 0.1, 1e-12)
            axes[1, 1].set_xlim(left=lower_limit, right=1.5)
        
        if legend_added_log:
            axes[1, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        axes[1, 1].set_title(f'Station {station_index:03d} - Species Concentrations (Log Scale)')
        
        # Plot 6: V profile (moved from position 5)
        axes[1, 2].plot(station_data['V'], eta, 'black', linewidth=2)
        axes[1, 2].set_xlabel('V (Velocity Component)')
        axes[1, 2].set_ylabel('η')
        axes[1, 2].grid(True, alpha=0.3)
        axes[1, 2].set_title(f'Station {station_index:03d} - V Profile')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight')
        else:
            plt.show()
        
        plt.close()  # Free memory
        return True
    
    def plot_individual_profile_F(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot velocity profile (F) individually for stagnation point"""
        try:
            valid_stations = self.reader.get_valid_stations()
            if station_index not in valid_stations:
                print(f"Error: Station {station_index} does not have required fields for plotting.")
                return False
            
            station_data = self.reader.get_station_data(station_index)
            
        except Exception as e:
            print(f"Error loading station {station_index}: {e}")
            return False
        
        # Create individual plot with same formatting as subplot
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
        eta = np.array(station_data['eta'])
        
        ax.plot(station_data['F'], eta, 'black', linewidth=2)
        ax.set_xlabel('F (Dimensionless Stream Function)')
        ax.set_ylabel('η (Similarity Coordinate)')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Station {station_index:03d} - Velocity Profile (F)', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight', dpi=600)
        else:
            plt.show()
        
        plt.close()
        return True
    
    def plot_individual_profile_g(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot enthalpy profile (g) individually for stagnation point"""
        try:
            valid_stations = self.reader.get_valid_stations()
            if station_index not in valid_stations:
                print(f"Error: Station {station_index} does not have required fields for plotting.")
                return False
            
            station_data = self.reader.get_station_data(station_index)
            
        except Exception as e:
            print(f"Error loading station {station_index}: {e}")
            return False
        
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
        eta = np.array(station_data['eta'])
        
        ax.plot(station_data['g'], eta, 'black', linewidth=2)
        ax.set_xlabel('g (Dimensionless Enthalpy)')
        ax.set_ylabel('η')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Station {station_index:03d} - Enthalpy Profile (g)', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight', dpi=600)
        else:
            plt.show()
        
        plt.close()
        return True
    
    def plot_individual_temperature(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot temperature profile individually for stagnation point"""
        try:
            valid_stations = self.reader.get_valid_stations()
            if station_index not in valid_stations:
                print(f"Error: Station {station_index} does not have required fields for plotting.")
                return False
            
            station_data = self.reader.get_station_data(station_index)
            
        except Exception as e:
            print(f"Error loading station {station_index}: {e}")
            return False
        
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
        eta = np.array(station_data['eta'])
        
        ax.plot(station_data['temperature'], eta, 'black', linewidth=2)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('η')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Station {station_index:03d} - Temperature Profile', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight', dpi=600)
        else:
            plt.show()
        
        plt.close()
        return True
    
    def plot_individual_species_linear(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot species concentrations (linear scale) individually for stagnation point"""
        try:
            valid_stations = self.reader.get_valid_stations()
            if station_index not in valid_stations:
                print(f"Error: Station {station_index} does not have required fields for plotting.")
                return False
            
            station_data = self.reader.get_station_data(station_index)
            species_names = self.reader.get_species_names()
            
        except Exception as e:
            print(f"Error loading station {station_index}: {e}")
            return False
        
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=300)
        eta = np.array(station_data['eta'])
        
        line_styles = ['-', '--', '-.', ':']
        markers = ['', 'o', 's', '^', 'v', 'D', 'p', '*', 'h', 'H', '+', 'x']
        legend_added = False
        
        for i, species in enumerate(species_names):
            col_name = f'c_{species}'
            if col_name in station_data:
                style_idx = i % len(line_styles)
                marker_idx = i % len(markers)
                
                ax.plot(station_data[col_name], eta, 
                       color='black',
                       linestyle=line_styles[style_idx],
                       marker=markers[marker_idx] if markers[marker_idx] else None,
                       markevery=max(1, len(eta)//8),
                       markersize=3.5,
                       linewidth=1.5, 
                       label=species)
                legend_added = True
        
        ax.set_xlabel('Mass Fraction')
        ax.set_ylabel('η')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Station {station_index:03d} - Species Concentrations', fontsize=12, fontweight='bold')
        
        if legend_added:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight', dpi=600)
        else:
            plt.show()
        
        plt.close()
        return True
    
    def plot_individual_species_log(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot species concentrations (log scale) individually for stagnation point"""
        try:
            valid_stations = self.reader.get_valid_stations()
            if station_index not in valid_stations:
                print(f"Error: Station {station_index} does not have required fields for plotting.")
                return False
            
            station_data = self.reader.get_station_data(station_index)
            species_names = self.reader.get_species_names()
            
        except Exception as e:
            print(f"Error loading station {station_index}: {e}")
            return False
        
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=300)
        eta = np.array(station_data['eta'])
        
        line_styles = ['-', '--', '-.', ':']
        markers = ['', 'o', 's', '^', 'v', 'D', 'p', '*', 'h', 'H', '+', 'x']
        legend_added_log = False
        min_positive_value = 1.0
        
        for i, species in enumerate(species_names):
            col_name = f'c_{species}'
            if col_name in station_data:
                data = station_data[col_name]
                positive_mask = data > 0
                if np.any(positive_mask):
                    min_positive_value = min(min_positive_value, np.min(data[positive_mask]))
                    
                    style_idx = i % len(line_styles)
                    marker_idx = i % len(markers)
                    
                    ax.semilogx(data[positive_mask], eta[positive_mask], 
                               color='black',
                               linestyle=line_styles[style_idx],
                               marker=markers[marker_idx] if markers[marker_idx] else None,
                               markevery=max(1, np.sum(positive_mask)//8),
                               markersize=3.5,
                               linewidth=1.5, 
                               label=species)
                    legend_added_log = True
        
        ax.set_xlabel('Mass Fraction (log scale)')
        ax.set_ylabel('η')
        ax.grid(True, alpha=0.3, which='both')
        ax.set_title(f'Station {station_index:03d} - Species Concentrations (Log Scale)', fontsize=12, fontweight='bold')
        
        if legend_added_log and min_positive_value < 1.0:
            # Set lower limit to one decade below the minimum value, but not below 1e-12
            lower_limit = max(min_positive_value * 0.1, 1e-12)
            ax.set_xlim(left=lower_limit, right=1.5)
        
        if legend_added_log:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight', dpi=600)
        else:
            plt.show()
        
        plt.close()
        return True
    
    def plot_individual_profile_V(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot V profile individually for stagnation point"""
        try:
            valid_stations = self.reader.get_valid_stations()
            if station_index not in valid_stations:
                print(f"Error: Station {station_index} does not have required fields for plotting.")
                return False
            
            station_data = self.reader.get_station_data(station_index)
            
        except Exception as e:
            print(f"Error loading station {station_index}: {e}")
            return False
        
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
        eta = np.array(station_data['eta'])
        
        ax.plot(station_data['V'], eta, 'black', linewidth=2)
        ax.set_xlabel('V (Velocity Component)')
        ax.set_ylabel('η')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Station {station_index:03d} - V Profile', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight', dpi=600)
        else:
            plt.show()
        
        plt.close()
        return True
    
    def plot_individual_stagnation_all(self, output_dir: Path) -> bool:
        """Generate all individual plots for stagnation point (station 0) simultaneously"""
        station_index = 0
        
        # Check if station 0 is valid for profiles
        valid_stations = self.reader.get_valid_stations()
        valid_heat_flux_stations = self.reader.get_valid_heat_flux_stations()
        
        if station_index not in valid_stations:
            print(f"Error: Station {station_index} (stagnation point) not available for profile plotting.")
            print(f"Available stations: {sorted(valid_stations)}")
            return False
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        print(f"\nGenerating individual stagnation point plots in: {output_dir}")
        print("-" * 60)
        
        # Generate all plots simultaneously
        plots = [
            ('F_profile', self.plot_individual_profile_F),
            ('g_profile', self.plot_individual_profile_g), 
            ('temperature', self.plot_individual_temperature),
            ('species_linear', self.plot_individual_species_linear),
            ('species_log', self.plot_individual_species_log),
            ('V_profile', self.plot_individual_profile_V)
        ]
        
        # Add heat flux plots if data is available
        if station_index in valid_heat_flux_stations:
            plots.extend([
                ('heat_flux_dimensional', self.plot_individual_heat_flux_dimensional),
                ('heat_flux_species', self.plot_individual_heat_flux_species)
            ])
        
        success_count = 0
        for plot_name, plot_func in plots:
            save_path = output_dir / f'stagnation_{plot_name}'
            if plot_func(station_index, save_path):
                print(f"   [OK] {plot_name}: stagnation_{plot_name}.pdf")
                success_count += 1
            else:
                print(f"   [FAILED] {plot_name}")
        
        print("-" * 60)
        print(f"Individual stagnation point plots completed: {success_count}/{len(plots)} successful")
        
        return success_count == len(plots)

    def plot_individual_heat_flux_dimensional(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot dimensional heat flux profiles for a given station"""
        try:
            # Check if station is valid for heat flux plotting
            valid_heat_flux_stations = self.reader.get_valid_heat_flux_stations()
            if station_index not in valid_heat_flux_stations:
                return False
            
            station_data = self.reader.get_station_heat_flux_data(station_index)
            
        except Exception as e:
            print(f"Error loading heat flux data for station {station_index}: {e}")
            return False
        
        # Create single plot for dimensional flux profiles
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
        # Get eta coordinates and heat flux data
        eta = np.array(station_data['eta'])
        heat_flux = station_data['heat_flux']
        
        # Dimensional flux profiles
        q_cond_dim = np.array(heat_flux['dimensional']['q_conductive'])
        q_diff_dim = np.array(heat_flux['dimensional']['q_diffusive']) 
        q_total_dim = np.array(heat_flux['dimensional']['q_total'])
        
        ax.plot(q_cond_dim, eta, 'black', linestyle='-', linewidth=2, label='Conductive')
        ax.plot(q_diff_dim, eta, 'black', linestyle='--', linewidth=2, label='Diffusive')
        ax.plot(q_total_dim, eta, 'black', linestyle='-.', linewidth=2, label='Total')
        ax.set_xlabel('Heat Flux (W/m²)', fontsize=12)
        ax.set_ylabel('η', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10)
        ax.set_title(f'Station {station_index:03d} - Dimensional Heat Flux Profiles', fontsize=14)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        ax.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(f"{save_path}.pdf", format='pdf', bbox_inches='tight')
        else:
            plt.show()
        
        plt.close()
        return True
    
    def plot_individual_heat_flux_species(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot species diffusive heat flux contributions for a given station"""
        try:
            # Check if station is valid for heat flux plotting
            valid_heat_flux_stations = self.reader.get_valid_heat_flux_stations()
            if station_index not in valid_heat_flux_stations:
                return False
            
            station_data = self.reader.get_station_heat_flux_data(station_index)
            species_names = self.reader.get_species_names()
            
        except Exception as e:
            print(f"Error loading heat flux data for station {station_index}: {e}")
            return False
        
        # Create single plot for species contributions
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
        # Get eta coordinates and heat flux data
        eta = np.array(station_data['eta'])
        heat_flux = station_data['heat_flux']
        
        # Dimensional species diffusive contributions
        if 'q_diffusive_species' in heat_flux['dimensional']:
            q_diff_species_dim = heat_flux['dimensional']['q_diffusive_species']
            line_styles = ['-', '--', '-.', ':']
            markers = ['', 'o', 's', '^', 'v', 'D', 'p', '*', 'h', 'H', '+', 'x']
            
            for i, species in enumerate(species_names):
                if i < q_diff_species_dim.shape[0]:
                    style_idx = i % len(line_styles)
                    marker_idx = i % len(markers)
                    
                    ax.plot(q_diff_species_dim[i, :], eta,
                           color='black',
                           linestyle=line_styles[style_idx],
                           marker=markers[marker_idx] if markers[marker_idx] else None,
                           markevery=max(1, len(eta)//8),
                           markersize=4,
                           linewidth=1.5,
                           label=species)
            
            ax.legend(fontsize=10)
        
        ax.set_xlabel('Species Diffusive Heat Flux (W/m²)', fontsize=12)
        ax.set_ylabel('η', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Station {station_index:03d} - Species Heat Flux Contributions', fontsize=14)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        ax.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(f"{save_path}.pdf", format='pdf', bbox_inches='tight')
        else:
            plt.show()
        
        plt.close()
        return True
    
    def plot_heat_flux(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot heat flux analysis for a given station
        
        Returns:
            bool: True if plot was created successfully, False otherwise
        """
        
        try:
            # Check if station is valid for heat flux plotting
            valid_heat_flux_stations = self.reader.get_valid_heat_flux_stations()
            if station_index not in valid_heat_flux_stations:
                available = sorted(valid_heat_flux_stations)
                print(f"Error: Station {station_index} does not have required fields for heat flux plotting.")
                print(f"Available stations with heat flux data: {available}")
                return False
            
            station_data = self.reader.get_station_heat_flux_data(station_index)
            species_names = self.reader.get_species_names()
            
        except Exception as e:
            print(f"Error loading heat flux data for station {station_index}: {e}")
            return False
        
        # Create subplots with clear titles
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'Heat Flux Analysis - Station {station_index:03d}', 
                    fontsize=16, fontweight='bold')
        
        # Get eta coordinates and heat flux data
        eta = np.array(station_data['eta'])
        heat_flux = station_data['heat_flux']
        
        # [1,1] - Dimensional flux profiles
        q_cond_dim = np.array(heat_flux['dimensional']['q_conductive'])
        q_diff_dim = np.array(heat_flux['dimensional']['q_diffusive']) 
        q_total_dim = np.array(heat_flux['dimensional']['q_total'])
        
        axes[0, 0].plot(q_cond_dim, eta, 'black', linestyle='-', linewidth=2, label='Conductive')
        axes[0, 0].plot(q_diff_dim, eta, 'black', linestyle='--', linewidth=2, label='Diffusive')
        axes[0, 0].plot(q_total_dim, eta, 'black', linestyle='-.', linewidth=2, label='Total')
        axes[0, 0].set_xlabel('Heat Flux (W/m²)')
        axes[0, 0].set_ylabel('η')
        axes[0, 0].grid(True, alpha=0.3)
        axes[0, 0].legend()
        axes[0, 0].set_title(f'Station {station_index:03d} - Dimensional Flux Profiles')
        axes[0, 0].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        # [1,2] - Dimensional species diffusive contributions
        if 'q_diffusive_species' in heat_flux['dimensional']:
            q_diff_species_dim = heat_flux['dimensional']['q_diffusive_species']
            line_styles = ['-', '--', '-.', ':']
            markers = ['', 'o', 's', '^', 'v', 'D', 'p', '*', 'h', 'H', '+', 'x']
            
            legend_added = False
            for i, species in enumerate(species_names):
                if i < q_diff_species_dim.shape[0]:
                    style_idx = i % len(line_styles)
                    marker_idx = i % len(markers)
                    
                    axes[0, 1].plot(q_diff_species_dim[i, :], eta,
                                  color='black',
                                  linestyle=line_styles[style_idx],
                                  marker=markers[marker_idx] if markers[marker_idx] else None,
                                  markevery=max(1, len(eta)//8),
                                  markersize=3.5,
                                  linewidth=1.5,
                                  label=species)
                    legend_added = True
            
            if legend_added:
                axes[0, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        axes[0, 1].set_xlabel('Species Diffusive Heat Flux (W/m²)')
        axes[0, 1].set_ylabel('η')
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_title(f'Station {station_index:03d} - Species Contributions (Dim.)')
        
        axes[0, 1].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        axes[0, 1].tick_params(axis='x', rotation=45)
        
        # [1,3] - Focus on conductive vs diffusive (dimensional)
        axes[0, 2].plot(q_cond_dim, eta, 'black', linestyle='-', linewidth=2, label='Conductive')
        axes[0, 2].plot(q_diff_dim, eta, 'black', linestyle='--', linewidth=2, label='Diffusive')
        axes[0, 2].set_xlabel('Heat Flux (W/m²)')
        axes[0, 2].set_ylabel('η')
        axes[0, 2].grid(True, alpha=0.3)
        axes[0, 2].legend()
        axes[0, 2].set_title(f'Station {station_index:03d} - Conductive vs Diffusive (Dim.)')
        axes[0, 2].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        axes[0, 2].tick_params(axis='x', rotation=45)
        
        # [2,1] - Nondimensional flux profiles
        q_cond_nondim = np.array(heat_flux['nondimensional']['q_conductive'])
        q_diff_nondim = np.array(heat_flux['nondimensional']['q_diffusive'])
        q_total_nondim = np.array(heat_flux['nondimensional']['q_total'])
        
        axes[1, 0].plot(q_cond_nondim, eta, 'black', linestyle='-', linewidth=2, label='Conductive')
        axes[1, 0].plot(q_diff_nondim, eta, 'black', linestyle='--', linewidth=2, label='Diffusive') 
        axes[1, 0].plot(q_total_nondim, eta, 'black', linestyle='-.', linewidth=2, label='Total')
        axes[1, 0].set_xlabel('Heat Flux (-)')
        axes[1, 0].set_ylabel('η')
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].legend()
        axes[1, 0].set_title(f'Station {station_index:03d} - Nondimensional Flux Profiles')
        axes[1, 0].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        # [2,2] - Nondimensional species diffusive contributions
        if 'q_diffusive_species' in heat_flux['nondimensional']:
            q_diff_species_nondim = heat_flux['nondimensional']['q_diffusive_species']
            
            legend_added = False
            for i, species in enumerate(species_names):
                if i < q_diff_species_nondim.shape[0]:
                    style_idx = i % len(line_styles)
                    marker_idx = i % len(markers)
                    
                    axes[1, 1].plot(q_diff_species_nondim[i, :], eta,
                                  color='black',
                                  linestyle=line_styles[style_idx],
                                  marker=markers[marker_idx] if markers[marker_idx] else None,
                                  markevery=max(1, len(eta)//8),
                                  markersize=3.5,
                                  linewidth=1.5,
                                  label=species)
                    legend_added = True
            
            if legend_added:
                axes[1, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        axes[1, 1].set_xlabel('Species Diffusive Heat Flux (-)')
        axes[1, 1].set_ylabel('η')
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].set_title(f'Station {station_index:03d} - Species Contributions (Nondim.)')
        
        axes[1, 1].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        axes[1, 1].tick_params(axis='x', rotation=45)
        
        # [2,3] - Diffusive/Total ratio
        # Avoid division by zero
        q_ratio = np.where(np.abs(q_total_nondim) > 1e-12, 
                          q_diff_nondim / q_total_nondim, 
                          0.0)
        
        axes[1, 2].plot(q_ratio, eta, 'black', linewidth=2)
        axes[1, 2].set_xlabel('q_diffusive / q_total (-)')
        axes[1, 2].set_ylabel('η')
        axes[1, 2].grid(True, alpha=0.3)
        axes[1, 2].set_title(f'Station {station_index:03d} - Diffusive Contribution Ratio')
        axes[1, 2].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        axes[1, 2].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight')
        else:
            plt.show()
        
        plt.close()  # Free memory
        return True
    
    def plot_F_and_g_map(self, save_path: Optional[Path] = None) -> bool:
        """Plot 2D smoothly interpolated maps of F(η, ξ) and g(η, ξ)
        
        Returns:
            bool: True if plot was created successfully, False otherwise
        """
        
        if 'stations' not in self.reader.data:
            print("Error: No station data available for F and g mapping")
            return False
        
        valid_stations = self.reader.get_valid_stations()
        if len(valid_stations) < 2:
            print(f"Error: Need at least 2 valid stations for mapping, found {len(valid_stations)}")
            return False
        
        try:
            # Collect all F, g, eta, and xi data
            xi_vals = []
            eta_vals = []
            F_vals = []
            g_vals = []
            
            # Get xi coordinates from metadata - REQUIRED
            if not ('metadata' in self.reader.data and 'grid' in self.reader.data['metadata'] and 'xi_coordinates' in self.reader.data['metadata']['grid']):
                print("Error: xi coordinates not found in metadata. F and g mapping requires proper grid coordinates.")
                print("Required path: metadata/grid/xi_coordinates")
                return False
            
            xi_coordinates = np.array(self.reader.data['metadata']['grid']['xi_coordinates'])
            
            # Create station index to xi coordinate mapping
            station_to_xi = {}
            for i, xi_coord in enumerate(xi_coordinates):
                station_to_xi[i] = xi_coord
            
            # Collect data from valid stations only
            stations_used = []
            for station_key, station_data in self.reader.data['stations'].items():
                station_idx = int(station_key.split('_')[-1])
                
                if station_idx not in valid_stations:
                    continue
                    
                if station_idx in station_to_xi:
                    xi_station = station_to_xi[station_idx]
                    eta_station = np.array(station_data['eta'])
                    F_station = np.array(station_data['F'])
                    g_station = np.array(station_data['g'])
                    
                    stations_used.append(station_idx)
                    
                    # Add data points
                    for j in range(len(eta_station)):
                        xi_vals.append(xi_station)
                        eta_vals.append(eta_station[j])
                        F_vals.append(F_station[j])
                        g_vals.append(g_station[j])
            
            if len(stations_used) < 2:
                print("Error: Insufficient valid station data for interpolation")
                return False
            
            # Convert to numpy arrays
            xi_vals = np.array(xi_vals)
            eta_vals = np.array(eta_vals)
            F_vals = np.array(F_vals)
            g_vals = np.array(g_vals)
            
            # Create high-resolution interpolation grid
            xi_min, xi_max = xi_vals.min(), xi_vals.max()
            eta_min, eta_max = eta_vals.min(), eta_vals.max()
            
            xi_grid = np.linspace(xi_min, xi_max, 200)
            eta_grid = np.linspace(eta_min, eta_max, 200)
            Xi, Eta = np.meshgrid(xi_grid, eta_grid)
            
            # Interpolate F and g with cubic method
            points = np.column_stack((xi_vals, eta_vals))
            F_interp = griddata(points, F_vals, (Xi, Eta), method='cubic')
            g_interp = griddata(points, g_vals, (Xi, Eta), method='cubic')
            
            # Fill NaN values at boundaries using nearest neighbor
            F_interp_nearest = griddata(points, F_vals, (Xi, Eta), method='nearest')
            g_interp_nearest = griddata(points, g_vals, (Xi, Eta), method='nearest')
            
            # Replace NaN values
            F_interp = np.where(np.isnan(F_interp), F_interp_nearest, F_interp)
            g_interp = np.where(np.isnan(g_interp), g_interp_nearest, g_interp)
            
            # Apply conservative Gaussian smoothing
            F_interp = gaussian_filter(F_interp, sigma=0.6)
            g_interp = gaussian_filter(g_interp, sigma=0.6)
            
            # Create the plot
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
            fig.suptitle(f'F and g Mapping (Stations: {min(stations_used):03d}-{max(stations_used):03d})', 
                        fontsize=14, fontweight='bold')
            
            # Plot F(η, ξ) with smooth contours
            im1 = ax1.imshow(F_interp, extent=[xi_min, xi_max, eta_min, eta_max], 
                 aspect='auto', origin='lower', cmap='Greys_r', interpolation='bicubic')
            ax1.set_xlabel('ξ')
            ax1.set_ylabel('η')
            ax1.set_title('F(η, ξ)')
            cbar1 = plt.colorbar(im1, ax=ax1)
            cbar1.set_label('F')
            
            # Plot g(η, ξ) with smooth contours
            im2 = ax2.imshow(g_interp, extent=[xi_min, xi_max, eta_min, eta_max], 
                 aspect='auto', origin='lower', cmap='Greys_r', interpolation='bicubic')
            ax2.set_xlabel('ξ')
            ax2.set_ylabel('η')
            ax2.set_title('g(η, ξ)')
            cbar2 = plt.colorbar(im2, ax=ax2)
            cbar2.set_label('g')
            
            plt.tight_layout()
            
            if save_path:
                plt.savefig(save_path.with_suffix('.pdf'), format='pdf', bbox_inches='tight')
                print(f"   [OK] F and g map ({len(stations_used)} stations): {save_path.name}")
            else:
                plt.show()
            
            plt.close()
            return True
            
        except Exception as e:
            print(f"Error creating F and g map: {e}")
            return False
    
    def create_summary_report(self, output_dir: Path) -> None:
        """Create a comprehensive summary report with consistent naming"""
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        print(f"\nCreating summary report in: {output_dir}")
        print("-" * 60)
        
        valid_stations = sorted(self.reader.get_valid_stations())
        valid_heat_flux_stations = sorted(self.reader.get_valid_heat_flux_stations())
        
        if not valid_stations and not valid_heat_flux_stations:
            print("Error: No valid stations found for plotting")
            return
        
        # Track successfully plotted stations
        plotted_stations = []
        plotted_heat_flux_stations = []
        
        # Plot profiles for valid stations (limit to first few to avoid too many plots)
        stations_to_plot = valid_stations[:5]  # Limit to first 5 valid stations
        
        for station_idx in stations_to_plot:
            profile_path = output_dir / f'profiles_station_{station_idx:03d}'
            if self.plot_profiles(station_idx, save_path=profile_path):
                plotted_stations.append(station_idx)
                print(f"   [OK] Profile station {station_idx:03d}: {profile_path.name}")
        
        # Plot heat flux for valid stations (same logic: first 5)
        heat_flux_stations_to_plot = valid_heat_flux_stations[:5]
        
        for station_idx in heat_flux_stations_to_plot:
            heat_flux_path = output_dir / f'heat_flux_station_{station_idx:03d}'
            if self.plot_heat_flux(station_idx, save_path=heat_flux_path):
                plotted_heat_flux_stations.append(station_idx)
                print(f"   [OK] Heat flux station {station_idx:03d}: {heat_flux_path.name}")
        
        # Create F and g map if possible
        f_g_map_path = output_dir / 'F_g_eta_xi_map.png'
        f_g_map_created = self.plot_F_and_g_map(save_path=f_g_map_path)
        
        # Create summary data file with consistent information
        try:
            summary_data = {
                'plotted_stations': plotted_stations,
                'plotted_heat_flux_stations': plotted_heat_flux_stations,
                'f_g_map_created': f_g_map_created,
                'valid_stations_total': len(valid_stations),
                'valid_heat_flux_stations_total': len(valid_heat_flux_stations),
                'all_valid_stations': valid_stations,
                'all_valid_heat_flux_stations': valid_heat_flux_stations
            }
            self._create_summary_json(output_dir / 'summary.json', summary_data)
            print(f"   [OK] Summary JSON: summary.json")
        except Exception as e:
            print(f"Warning: Could not create summary JSON file: {e}")
        
        print("-" * 60)
        print(f"\nSummary report completed:")
        print(f"   - Profile plots: {len(plotted_stations)} stations")
        print(f"   - Heat flux plots: {len(plotted_heat_flux_stations)} stations")
        if f_g_map_created:
            print(f"   - F and g mapping: created")
    
    def _create_summary_json(self, file_path: Path, plot_summary: Dict) -> None:
        """Create JSON summary with consistent plot information"""
        
        summary = {
            'metadata': self.reader.metadata,
            'file_info': {
                'input_file': str(self.reader.input_path),
                'format': self.reader.format,
                'total_stations_in_file': len(self.reader.data['stations']) if self.reader.data else 0
            },
            'data_validation': {
                'valid_stations_for_plotting': len(self.reader.get_valid_stations()),
                'valid_heat_flux_stations_for_plotting': len(self.reader.get_valid_heat_flux_stations()),
                'required_profile_fields': list(self.reader.REQUIRED_PROFILE_FIELDS),
                'required_heat_flux_fields': list(self.reader.REQUIRED_HEAT_FLUX_FIELDS),
                'valid_station_indices': sorted(self.reader.get_valid_stations()),
                'valid_heat_flux_station_indices': sorted(self.reader.get_valid_heat_flux_stations())
            },
            'species_info': {
                'species_names': self.reader.get_species_names(),
                'species_count': len(self.reader.get_species_names())
            },
            'plots_generated': plot_summary
        }

        # Convert numpy arrays to lists for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, (bytes, np.bytes_)):
                try:
                    return obj.decode('utf-8')
                except UnicodeDecodeError:
                    return str(obj)
            elif isinstance(obj, dict):
                return {key: convert_numpy(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(item) for item in obj]
            elif isinstance(obj, set):
                return sorted(list(obj))
            return obj
        
        summary = convert_numpy(summary)
        
        with open(file_path, 'w') as f:
            json.dump(summary, f, indent=2)
        


def main():
    """Main function with command-line interface and robust error handling"""
    
    parser = argparse.ArgumentParser(description='BLAST Post-Processing Tool - Enhanced with Heat Flux Analysis')
    parser.add_argument('--input', '-i', type=str, required=True,
                       help='Input HDF5 file')
    parser.add_argument('--plots', '-p', type=str, default='all',
                       choices=['all', 'profiles', 'heat_flux', 'summary', 'f_g_map', 'stagnation_individual'],
                       help='Type of plots to generate')
    parser.add_argument('--station', '-s', type=int, default=0,
                       help='Station index for profile plots')
    parser.add_argument('--output', '-o', type=str, default='blast_plots',
                       help='Output directory for plots')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for saved figures')
    parser.add_argument('--show', action='store_true',
                       help='Show plots interactively instead of saving')
    
    args = parser.parse_args()
    
    # Initialize reader with error handling
    try:
        reader = BLASTReader(args.input)
        
        print(f"Initializing {reader.format.upper()} reader for: {reader.input_path}")
        reader.load_data()
        
        valid_stations = reader.get_valid_stations()
        valid_heat_flux_stations = reader.get_valid_heat_flux_stations()
        
        print(f"\nData loaded successfully:")
        print(f"   - Total stations: {len(reader.data['stations'])}")
        print(f"   - Valid for profile plotting: {len(valid_stations)}")
        print(f"   - Valid for heat flux plotting: {len(valid_heat_flux_stations)}")
        
        if not valid_stations and not valid_heat_flux_stations:
            print("Error: No stations have the required fields for any type of plotting")
            return 1
        
    except Exception as e:
        print(f"Error initializing data reader: {e}")
        return 1
    
    # Initialize plotter
    plotter = BLASTPlotter(reader)
    plotter.dpi = args.dpi
    
    # Create output directory
    output_path = Path(args.output) if not args.show else None
    if output_path:
        output_path.mkdir(exist_ok=True)
        print(f"\nOutput directory: {output_path}")
        print("-" * 60)
    
    # Generate plots with error handling
    try:
        if args.plots == 'all':
            if args.show:
                # Show first available plot type for the requested station
                if args.station in valid_stations:
                    success = plotter.plot_profiles(args.station)
                elif args.station in valid_heat_flux_stations:
                    success = plotter.plot_heat_flux(args.station)
                else:
                    print(f"Error: Station {args.station} not available for plotting")
                    return 1
                if not success:
                    return 1
            else:
                plotter.create_summary_report(output_path)
                
        elif args.plots == 'profiles':
            if args.station not in valid_stations:
                print(f"Error: Station {args.station} not valid for profile plotting")
                print(f"Available profile stations: {sorted(valid_stations)}")
                return 1
                
            save_path = output_path / f'profiles_station_{args.station:03d}.png' if output_path else None
            success = plotter.plot_profiles(args.station, save_path)
            if not success:
                return 1
        
        elif args.plots == 'heat_flux':
            if args.station not in valid_heat_flux_stations:
                print(f"Error: Station {args.station} not valid for heat flux plotting")
                print(f"Available heat flux stations: {sorted(valid_heat_flux_stations)}")
                return 1
                
            save_path = output_path / f'heat_flux_station_{args.station:03d}.png' if output_path else None
            success = plotter.plot_heat_flux(args.station, save_path)
            if not success:
                return 1
                
        elif args.plots == 'summary':
            if not output_path:
                output_path = Path('blast_summary')
            plotter.create_summary_report(output_path)
            
        elif args.plots == 'f_g_map':
            save_path = output_path / 'F_g_eta_xi_map.png' if output_path else None
            success = plotter.plot_F_and_g_map(save_path)
            if not success:
                return 1
        
        elif args.plots == 'stagnation_individual':
            if not output_path:
                output_path = Path('blast_stagnation_individual')
            success = plotter.plot_individual_stagnation_all(output_path)
            if not success:
                return 1
        
        print("\n" + "="*60)
        print("Post-processing completed successfully!")
        print("="*60)
        
    except Exception as e:
        print(f"Error during plotting: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
    
"""     
Usage examples:

# Generate all plots (profiles + heat flux for first 5 stations + F/g map)
python3 postprocess_blast.py --input simulation.h5 --plots all --output results

# Generate only heat flux analysis for a specific station  
python3 postprocess_blast.py --input simulation.h5 --plots heat_flux --station 2 --output results

# Generate only profiles for a specific station
python3 postprocess_blast.py --input simulation.h5 --plots profiles --station 0 --output results

# Generate individual plots for stagnation point (station 0) - 6 separate PDF files
python3 postprocess_blast.py --input simulation.h5 --plots stagnation_individual --output results

# Interactive display of heat flux for station 1
python3 postprocess_blast.py --input simulation.h5 --plots heat_flux --station 1 --show
"""