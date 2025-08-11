#!/usr/bin/env python3
"""
BLAST Boundary Layer Post-Processing Script — Temperature Focus

This script mirrors the previous heat-flux oriented tool, but focuses on temperatures:
- Physical temperature profile T(η)
- Modal temperatures (e.g., Tv, Tr, Te...) if present under stations/station_xxx/modal_temperatures/*
- Derived diagnostics (ΔT_mode = T_mode - T, ratios T_mode/T, gradients)

Usage:
    python postprocess_blast_temperatures.py --input simulation.h5 --plots all
    python postprocess_blast_temperatures.py --input simulation.h5 --plots temperatures --station 0
    python postprocess_blast_temperatures.py --input simulation.h5 --plots temp_map
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

# Set publication-quality plot style (same as the original)
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
    """Unified reader for BLAST output formats with robust error handling (temperature focus)"""

    REQUIRED_TEMP_FIELDS = {'eta', 'temperature'}  # minimal set for temperature plotting

    def __init__(self, input_path: Union[str, Path]):
        self.input_path = Path(input_path)
        self.format = self._detect_format()
        self.data = None
        self.metadata = None
        self.valid_temp_stations = set()

    def _detect_format(self) -> str:
        if not self.input_path.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_path}")
        if self.input_path.suffix.lower() in ['.h5', '.hdf5']:
            return 'hdf5'
        raise ValueError(f"Only HDF5 format is supported. Got: {self.input_path}")

    def _verify_file_integrity(self) -> None:
        try:
            with h5py.File(self.input_path, 'r') as f:
                if 'stations' not in f:
                    raise BLASTFileError("Missing 'stations' group in HDF5 file")
                if len(f['stations'].keys()) == 0:
                    raise BLASTFileError("No station data found in 'stations' group")
                if 'metadata' not in f:
                    print("Warning: No 'metadata' group found - some features may be limited")
                print(f"File integrity check passed: {len(f['stations'])} stations found")
        except Exception as e:
            if isinstance(e, BLASTFileError):
                raise
            raise BLASTFileError(f"File integrity check failed: {e}")

    def load_data(self) -> Dict:
        self._verify_file_integrity()
        if self.format == 'hdf5':
            return self._load_hdf5()
        else:
            raise NotImplementedError(f"Format {self.format} not implemented")

    def _load_hdf5(self) -> Dict:
        try:
            with h5py.File(self.input_path, 'r') as f:
                data = {}

                # metadata (optional)
                data['metadata'] = self._read_hdf5_group(f['metadata']) if 'metadata' in f else {}

                # stations
                data['stations'] = {}
                for station_name in f['stations'].keys():
                    try:
                        station_data = self._read_hdf5_group(f['stations'][station_name])
                        data['stations'][station_name] = station_data
                        idx = int(station_name.split('_')[-1])
                        if self._validate_temperature_data(station_data):
                            self.valid_temp_stations.add(idx)
                    except Exception as e:
                        print(f"Warning: Failed to load station {station_name}: {e}")

                if len(data['stations']) == 0:
                    raise BLASTFileError("No valid station data could be loaded")

                print(f"Loaded {len(data['stations'])} stations")
                print(f"- {len(self.valid_temp_stations)} valid for temperature plotting")

                self.data = data
                self.metadata = data.get('metadata', {})
                return data
        except Exception as e:
            if isinstance(e, BLASTFileError):
                raise
            raise RuntimeError(f"Failed to load HDF5 file: {e}")

    def _validate_temperature_data(self, station_data: Dict) -> bool:
        try:
            available = set(station_data.keys())
            if not self.REQUIRED_TEMP_FIELDS.issubset(available):
                return False
            eta_len = len(station_data['eta'])
            if eta_len == 0 or len(station_data['temperature']) != eta_len:
                return False
            # modal temps are optional; if present, ensure consistent lengths
            if 'modal_temperatures' in station_data and isinstance(station_data['modal_temperatures'], dict):
                for name, arr in station_data['modal_temperatures'].items():
                    # skip *_attrs keys
                    if name.endswith('_attrs'):
                        continue
                    if not hasattr(arr, '__len__'):
                        return False
                    if len(arr) != eta_len:
                        return False
            return True
        except Exception:
            return False

    def _read_hdf5_group(self, group) -> Dict:
        result = {}
        for key, item in group.items():
            if isinstance(item, h5py.Group):
                result[key] = self._read_hdf5_group(item)
            elif isinstance(item, h5py.Dataset):
                if h5py.check_string_dtype(item.dtype):
                    data = item[()]
                    if isinstance(data, np.ndarray):
                        result[key] = [s.decode('utf-8') if isinstance(s, bytes) else s for s in data]
                    else:
                        result[key] = data.decode('utf-8') if isinstance(data, bytes) else data
                else:
                    result[key] = np.array(item)
                if item.attrs:
                    result[f'{key}_attrs'] = dict(item.attrs)
        return result

    def get_valid_temperature_stations(self) -> Set[int]:
        if self.data is None:
            self.load_data()
        return self.valid_temp_stations.copy()

    def get_station_temperature_data(self, station_index: int) -> Dict:
        if self.data is None:
            self.load_data()
        key = f'station_{station_index:03d}'
        if key not in self.data['stations']:
            available = sorted([int(k.split('_')[-1]) for k in self.data['stations'].keys()])
            raise KeyError(f"Station {station_index} not found. Available stations: {available}")
        station = self.data['stations'][key].copy()
        if not self._validate_temperature_data(station):
            raise ValueError(f"Station {station_index} missing required temperature fields or lengths mismatch")
        # Build a dict of modal temps {name: array}
        modal = {}
        if 'modal_temperatures' in station and isinstance(station['modal_temperatures'], dict):
            for name, arr in station['modal_temperatures'].items():
                if name.endswith('_attrs'):
                    continue
                modal[str(name)] = np.array(arr)
        return {
            'eta': np.array(station['eta']),
            'temperature': np.array(station['temperature']),
            'modal_temperatures': modal
        }

    def get_species_names(self) -> List[str]:
        # Not central here, but kept for symmetry with the original API
        if self.metadata and 'mixture' in self.metadata and 'species_names' in self.metadata['mixture']:
            names = self.metadata['mixture']['species_names']
            if isinstance(names, (list, np.ndarray)):
                return [str(name) for name in names]
            else:
                return [str(names)]
        return []

class BLASTPlotter:
    """Publication-quality plotting for temperatures with robust error handling"""

    def __init__(self, reader: BLASTReader):
        self.reader = reader
        self.fig_size = (12, 8)
        self.dpi = 300

    def plot_temperatures(self, station_index: int = 0, save_path: Optional[Path] = None) -> bool:
        """Plot temperature profiles (physical and modal) for a given station"""
        try:
            valid = self.reader.get_valid_temperature_stations()
            if station_index not in valid:
                print(f"Error: Station {station_index} not valid for temperature plotting.")
                print(f"Available temperature stations: {sorted(valid)}")
                return False
            s = self.reader.get_station_temperature_data(station_index)
        except Exception as e:
            print(f"Error loading station {station_index}: {e}")
            return False

        eta = s['eta']
        T = s['temperature']
        modal = s['modal_temperatures']  # dict name->array

        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'Temperature Profiles - Station {station_index:03d}', fontsize=16, fontweight='bold')

        # [0,0] Physical temperature T(η)
        axes[0, 0].plot(T, eta, 'black', linewidth=2, label='T')
        axes[0, 0].set_xlabel('Temperature (K)')
        axes[0, 0].set_ylabel('η')
        axes[0, 0].grid(True, alpha=0.3)
        axes[0, 0].set_title(f'Station {station_index:03d} - Physical Temperature')
        axes[0, 0].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        axes[0, 0].tick_params(axis='x', rotation=45)

        # [0,1] Modal temperatures overlay
        line_styles = ['-', '--', '-.', ':']
        markers = ['', 'o', 's', '^', 'v', 'D', 'p', '*', 'h', 'H', '+', 'x']
        legend_added = False
        modal_names_sorted = sorted(modal.keys())  # stable order
        for i, name in enumerate(modal_names_sorted):
            arr = modal[name]
            style_idx = i % len(line_styles)
            marker_idx = i % len(markers)
            axes[0, 1].plot(arr, eta,
                            color='black',
                            linestyle=line_styles[style_idx],
                            marker=markers[marker_idx] if markers[marker_idx] else None,
                            markevery=max(1, len(eta)//8),
                            markersize=3.5,
                            linewidth=1.5,
                            label=name)
            legend_added = True
        axes[0, 1].set_xlabel('Temperature (K)')
        axes[0, 1].set_ylabel('η')
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_title(f'Station {station_index:03d} - Modal Temperatures')
        axes[0, 1].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
        axes[0, 1].tick_params(axis='x', rotation=45)
        if legend_added:
            axes[0, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        # [0,2] Differences ΔT_mode = T_mode - T
        if modal_names_sorted:
            for i, name in enumerate(modal_names_sorted):
                arr = modal[name]
                style_idx = i % len(line_styles)
                marker_idx = i % len(markers)
                axes[0, 2].plot(arr - T, eta,
                                color='black',
                                linestyle=line_styles[style_idx],
                                marker=markers[marker_idx] if markers[marker_idx] else None,
                                markevery=max(1, len(eta)//8),
                                markersize=3.5,
                                linewidth=1.5,
                                label=f'{name} - T')
            axes[0, 2].set_xlabel('ΔT (K)')
            axes[0, 2].set_ylabel('η')
            axes[0, 2].grid(True, alpha=0.3)
            axes[0, 2].set_title(f'Station {station_index:03d} - Modal Deviations')
            axes[0, 2].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
            axes[0, 2].tick_params(axis='x', rotation=45)
            axes[0, 2].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            axes[0, 2].axis('off')

        # [1,0] Gradient dT/dη (central differences)
        if len(eta) >= 3:
            dT = np.gradient(T, eta)
            axes[1, 0].plot(dT, eta, 'black', linewidth=2)
            axes[1, 0].set_xlabel('dT/dη (K)')
            axes[1, 0].set_ylabel('η')
            axes[1, 0].grid(True, alpha=0.3)
            axes[1, 0].set_title(f'Station {station_index:03d} - Temperature Gradient')
            axes[1, 0].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
            axes[1, 0].tick_params(axis='x', rotation=45)
        else:
            axes[1, 0].axis('off')

        # [1,1] Ratios T_mode / T (if T>0 to avoid division by ~0)
        if modal_names_sorted:
            safe_T = np.where(np.abs(T) > 1e-12, T, np.nan)
            for i, name in enumerate(modal_names_sorted):
                arr = modal[name]
                style_idx = i % len(line_styles)
                marker_idx = i % len(markers)
                axes[1, 1].plot(arr / safe_T, eta,
                                color='black',
                                linestyle=line_styles[style_idx],
                                marker=markers[marker_idx] if markers[marker_idx] else None,
                                markevery=max(1, len(eta)//8),
                                markersize=3.5,
                                linewidth=1.5,
                                label=f'{name}/T')
            axes[1, 1].set_xlabel('Ratio (-)')
            axes[1, 1].set_ylabel('η')
            axes[1, 1].grid(True, alpha=0.3)
            axes[1, 1].set_title(f'Station {station_index:03d} - Modal Ratios')
            axes[1, 1].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
            axes[1, 1].tick_params(axis='x', rotation=45)
            axes[1, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            axes[1, 1].axis('off')

        # [1,2] placeholder (kept off to preserve grid symmetry as in the original)
        axes[1, 2].axis('off')

        plt.tight_layout()
        if save_path:
            plt.savefig(Path(save_path).with_suffix('.pdf'), format='pdf', bbox_inches='tight')
            print(f"Temperature plots for station {station_index:03d} saved to: {save_path}")
        else:
            plt.show()
        plt.close()
        return True

    def plot_temperature_map(self, save_path: Optional[Path] = None) -> bool:
        """Plot 2D interpolated maps of T(η, ξ) and, if available, the first modal temperature"""
        if self.reader.data is None:
            print("Error: No data loaded")
            return False

        valid = sorted(self.reader.get_valid_temperature_stations())
        if len(valid) < 2:
            print(f"Error: Need at least 2 valid stations for mapping, found {len(valid)}")
            return False

        try:
            # xi coordinates from metadata
            if not ('metadata' in self.reader.data and
                    'grid' in self.reader.data['metadata'] and
                    'xi_coordinates' in self.reader.data['metadata']['grid']):
                print("Error: xi coordinates not found in metadata (metadata/grid/xi_coordinates)")
                return False
            xi_coords = np.array(self.reader.data['metadata']['grid']['xi_coordinates'])

            # Collect (xi, eta, T) and optionally first modal
            xi_vals, eta_vals, T_vals = [], [], []
            modal_name = None
            modal_vals = []

            # Determine a consistent first modal name (if any)
            # We scan valid stations and pick the first station that has modal temps
            for st_idx in valid:
                key = f'station_{st_idx:03d}'
                st = self.reader.data['stations'][key]
                if 'modal_temperatures' in st and isinstance(st['modal_temperatures'], dict):
                    candidates = sorted([n for n in st['modal_temperatures'].keys() if not n.endswith('_attrs')])
                    if candidates:
                        modal_name = candidates[0]
                        break

            stations_used = []
            for st_idx in valid:
                key = f'station_{st_idx:03d}'
                st = self.reader.data['stations'][key]
                if st_idx >= len(xi_coords):
                    continue
                xi = xi_coords[st_idx]
                eta = np.array(st['eta'])
                T = np.array(st['temperature'])
                for j in range(len(eta)):
                    xi_vals.append(xi)
                    eta_vals.append(eta[j])
                    T_vals.append(T[j])

                if modal_name and 'modal_temperatures' in st and modal_name in st['modal_temperatures']:
                    m = np.array(st['modal_temperatures'][modal_name])
                    for j in range(len(eta)):
                        modal_vals.append(m[j])
                stations_used.append(st_idx)

            if len(stations_used) < 2:
                print("Error: Insufficient valid station data for interpolation")
                return False

            xi_vals = np.array(xi_vals); eta_vals = np.array(eta_vals); T_vals = np.array(T_vals)
            points = np.column_stack((xi_vals, eta_vals))

            xi_min, xi_max = xi_vals.min(), xi_vals.max()
            eta_min, eta_max = eta_vals.min(), eta_vals.max()

            xi_grid = np.linspace(xi_min, xi_max, 200)
            eta_grid = np.linspace(eta_min, eta_max, 200)
            Xi, Eta = np.meshgrid(xi_grid, eta_grid)

            # Interpolate T
            T_interp = griddata(points, T_vals, (Xi, Eta), method='cubic')
            T_near = griddata(points, T_vals, (Xi, Eta), method='nearest')
            T_interp = np.where(np.isnan(T_interp), T_near, T_interp)
            T_interp = gaussian_filter(T_interp, sigma=0.6)

            if modal_name and len(modal_vals) == len(T_vals):
                M_vals = np.array(modal_vals)
                M_interp = griddata(points, M_vals, (Xi, Eta), method='cubic')
                M_near = griddata(points, M_vals, (Xi, Eta), method='nearest')
                M_interp = np.where(np.isnan(M_interp), M_near, M_interp)
                M_interp = gaussian_filter(M_interp, sigma=0.6)

                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
                fig.suptitle(f'Temperature Mapping (Stations: {min(stations_used):03d}-{max(stations_used):03d})',
                             fontsize=14, fontweight='bold')

                im1 = ax1.imshow(T_interp, extent=[xi_min, xi_max, eta_min, eta_max],
                                 aspect='auto', origin='lower', cmap='Greys_r', interpolation='bicubic')
                ax1.set_xlabel('ξ'); ax1.set_ylabel('η'); ax1.set_title('T(η, ξ) [K]')
                cbar1 = plt.colorbar(im1, ax=ax1); cbar1.set_label('K')

                im2 = ax2.imshow(M_interp, extent=[xi_min, xi_max, eta_min, eta_max],
                                 aspect='auto', origin='lower', cmap='Greys_r', interpolation='bicubic')
                ax2.set_xlabel('ξ'); ax2.set_ylabel('η'); ax2.set_title(f'{modal_name}(η, ξ) [K]')
                cbar2 = plt.colorbar(im2, ax=ax2); cbar2.set_label('K')

            else:
                fig, ax1 = plt.subplots(1, 1, figsize=(7, 6))
                fig.suptitle(f'Temperature Mapping (Stations: {min(stations_used):03d}-{max(stations_used):03d})',
                             fontsize=14, fontweight='bold')
                im1 = ax1.imshow(T_interp, extent=[xi_min, xi_max, eta_min, eta_max],
                                 aspect='auto', origin='lower', cmap='Greys_r', interpolation='bicubic')
                ax1.set_xlabel('ξ'); ax1.set_ylabel('η'); ax1.set_title('T(η, ξ) [K]')
                cbar1 = plt.colorbar(im1, ax=ax1); cbar1.set_label('K')

            plt.tight_layout()
            if save_path:
                plt.savefig(Path(save_path).with_suffix('.pdf'), format='pdf', bbox_inches='tight')
                print(f"Temperature map saved to: {save_path}")
            else:
                plt.show()
            plt.close()
            return True

        except Exception as e:
            print(f"Error creating temperature map: {e}")
            return False

    def create_summary_report(self, output_dir: Path) -> None:
        """Create a summary report (first 5 stations by default)"""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

        print(f"Creating summary report in: {output_dir}")

        valid = sorted(self.reader.get_valid_temperature_stations())
        if not valid:
            print("Error: No valid stations found for temperature plotting")
            return

        plotted_stations = []
        stations_to_plot = valid[:5]
        for idx in stations_to_plot:
            p = output_dir / f'temperatures_station_{idx:03d}'
            if self.plot_temperatures(idx, save_path=p):
                plotted_stations.append(idx)

        # Map
        temp_map_path = output_dir / 'temperature_map.png'
        map_created = self.plot_temperature_map(save_path=temp_map_path)

        # JSON summary
        try:
            summary_data = {
                'plotted_temperature_stations': plotted_stations,
                'temp_map_created': map_created,
                'valid_temperature_stations_total': len(valid),
                'all_valid_temperature_stations': valid
            }
            self._create_summary_json(output_dir / 'summary_temperatures.json', summary_data)
        except Exception as e:
            print(f"Warning: Could not create summary JSON file: {e}")

        print("Summary report (temperatures) completed!")
        print(f"Successfully plotted {len(plotted_stations)} stations: {plotted_stations}")
        if map_created:
            print("Temperature mapping plot created successfully")

    def _create_summary_json(self, file_path: Path, plot_summary: Dict) -> None:
        summary = {
            'metadata': self.reader.metadata,
            'file_info': {
                'input_file': str(self.reader.input_path),
                'format': self.reader.format,
                'total_stations_in_file': len(self.reader.data['stations']) if self.reader.data else 0
            },
            'data_validation': {
                'valid_temperature_stations_for_plotting': len(self.reader.get_valid_temperature_stations()),
                'required_temperature_fields': list(BLASTReader.REQUIRED_TEMP_FIELDS),
                'valid_temperature_station_indices': sorted(self.reader.get_valid_temperature_stations())
            },
            'plots_generated': plot_summary
        }

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
                return {k: convert_numpy(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(x) for x in obj]
            elif isinstance(obj, set):
                return sorted(list(obj))
            return obj

        summary = convert_numpy(summary)
        with open(file_path, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"Summary JSON saved to: {file_path}")

def main():
    parser = argparse.ArgumentParser(description='BLAST Post-Processing Tool — Temperatures')
    parser.add_argument('--input', '-i', type=str, required=True, help='Input HDF5 file')
    parser.add_argument('--plots', '-p', type=str, default='all',
                        choices=['all', 'temperatures', 'summary', 'temp_map'],
                        help='Type of plots to generate')
    parser.add_argument('--station', '-s', type=int, default=0, help='Station index for temperature plots')
    parser.add_argument('--output', '-o', type=str, default='blast_plots_temperatures',
                        help='Output directory for plots')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for saved figures')
    parser.add_argument('--show', action='store_true', help='Show plots interactively instead of saving')

    args = parser.parse_args()

    try:
        reader = BLASTReader(args.input)
        print(f"Initializing {reader.format.upper()} reader for: {reader.input_path}")
        reader.load_data()
        valid = reader.get_valid_temperature_stations()
        print("Data loaded successfully:")
        print(f"  - Total stations: {len(reader.data['stations'])}")
        print(f"  - Valid stations for temperature plotting: {len(valid)}")
        print(f"  - Valid temperature station indices: {sorted(valid)}")
        if not valid:
            print("Error: No stations have the required fields for temperature plotting")
            return 1
    except Exception as e:
        print(f"Error initializing data reader: {e}")
        return 1

    plotter = BLASTPlotter(reader)
    plotter.dpi = args.dpi

    output_path = Path(args.output) if not args.show else None
    if output_path:
        output_path.mkdir(exist_ok=True)
        print(f"Output directory: {output_path}")

    try:
        if args.plots == 'all':
            if args.show:
                if args.station in valid:
                    ok = plotter.plot_temperatures(args.station)
                    if not ok:
                        return 1
                else:
                    print(f"Error: Station {args.station} not available for temperature plotting")
                    return 1
            else:
                plotter.create_summary_report(output_path)

        elif args.plots == 'temperatures':
            if args.station not in valid:
                print(f"Error: Station {args.station} not valid for temperature plotting")
                print(f"Available temperature stations: {sorted(valid)}")
                return 1
            save_path = output_path / f'temperatures_station_{args.station:03d}.png' if output_path else None
            ok = plotter.plot_temperatures(args.station, save_path)
            if not ok:
                return 1

        elif args.plots == 'summary':
            if not output_path:
                output_path = Path('blast_temperatures_summary')
            plotter.create_summary_report(output_path)

        elif args.plots == 'temp_map':
            save_path = output_path / 'temperature_map.png' if output_path else None
            ok = plotter.plot_temperature_map(save_path)
            if not ok:
                return 1

        print("Post-processing (temperatures) completed successfully!")
    except Exception as e:
        print(f"Error during plotting: {e}")
        return 1

    return 0

if __name__ == '__main__':
    exit(main())