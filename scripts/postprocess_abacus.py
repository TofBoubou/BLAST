#!/usr/bin/env python3
"""
BLAST Abacus Post-Processing Script (map + curves)

Reads abacus HDF5 produced by BLAST with datasets at root:
- temperatures: double[N]
- catalyticity_values: double[M]
- heat_fluxes: double[M, N]   # rows: gamma index, cols: Tw index

Produces a single figure with two subplots on one row:
- Left:  2D map of q_wall vs (gamma, T_w)
- Right: line plots q_wall(T_w) for each gamma (no extra interpolation)

Usage:
    python3 postprocess_abacus.py --input abacus.h5 --output abacus_plots/abacus_map
    python3 postprocess_abacus.py --input abacus.h5 --show
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py
from pathlib import Path
from typing import Union, Dict, Optional
import warnings

warnings.filterwarnings('ignore', message='Unable to import Axes3D', category=UserWarning)
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
    'grid.linewidth': 0.3,
    'grid.alpha': 0.3,
    'figure.dpi': 300,
    'savefig.dpi': 600,
})

class AbacusError(Exception):
    """Base exception for Abacus post-processing"""
    pass

class AbacusFileError(AbacusError):
    """I/O or structure errors in HDF5 file"""
    pass

class AbacusDataError(AbacusError):
    """Content/validation errors in datasets"""
    pass

class AbacusReader:
    """Minimal robust reader for abacus HDF5"""

    REQUIRED = ('temperatures', 'catalyticity_values', 'heat_fluxes')

    def __init__(self, input_path: Union[str, Path]):
        self.input_path = Path(input_path)
        self.data: Dict[str, np.ndarray] = {}
        self.attrs: Dict[str, Union[int, float]] = {}

    def _check(self) -> None:
        if not self.input_path.exists():
            raise AbacusFileError(f"Input file not found: {self.input_path}")
        if self.input_path.suffix.lower() not in ('.h5', '.hdf5'):
            raise AbacusFileError(f"Only HDF5 is supported. Got: {self.input_path}")

    def load(self) -> Dict[str, np.ndarray]:
        self._check()
        try:
            with h5py.File(self.input_path, 'r') as f:
                for k in self.REQUIRED:
                    if k not in f:
                        raise AbacusFileError(f"Missing dataset '{k}' at root")

                temps = np.array(f['temperatures'])
                gammas = np.array(f['catalyticity_values'])
                q = np.array(f['heat_fluxes'])

                if q.ndim != 2:
                    raise AbacusDataError("'heat_fluxes' must be 2D [M, N]")
                M, N = q.shape
                if temps.ndim != 1 or temps.size != N:
                    raise AbacusDataError("'temperatures'.size must equal heat_fluxes.shape[1]")
                if gammas.ndim != 1 or gammas.size != M:
                    raise AbacusDataError("'catalyticity_values'.size must equal heat_fluxes.shape[0]")

                self.data = {
                    'temperatures': temps,
                    'catalyticity_values': gammas,
                    'heat_fluxes': q
                }

                # Optional attributes
                self.attrs = {}
                for name in ('temperature_range', 'n_temperatures', 'n_catalyticity_values'):
                    if name in f.attrs:
                        self.attrs[name] = f.attrs[name]

                print(f"Loaded abacus: M={M} γ values, N={N} temperatures")
                return self.data

        except Exception as e:
            if isinstance(e, AbacusError):
                raise
            raise AbacusFileError(f"Failed to load HDF5: {e}")

class AbacusPlotter:
    """One figure with two subplots: map (left) + curves (right)"""

    def __init__(self, reader: AbacusReader, output_format: str = 'pdf', dpi: int = 300, mode: str = 'temperature_sweep'):
        self.reader = reader
        self.output_format = output_format
        self.dpi = dpi
        self.mode = mode

    def plot_individual_map(self, save_path: Optional[Path] = None) -> bool:
        """Plot only the 2D map of heat flux vs (gamma, T_w)"""
        if not self.reader.data:
            self.reader.load()

        T = self.reader.data['temperatures']           # (N,)
        G = self.reader.data['catalyticity_values']    # (M,)
        Q = self.reader.data['heat_fluxes']            # (M, N)

        fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=self.dpi)
        
        if self.mode == 'gamma_sweep':
            fig.suptitle('Abacus: Heat Flux Map (γ-sweep)', fontsize=14, fontweight='bold')
            
            if T.size >= 2:
                # 2D map with gamma on X-axis (log scale)
                G_grid, T_grid = np.meshgrid(G, T)
                im = ax.pcolormesh(
                    G_grid, T_grid, Q.T,
                    cmap='Greys_r',
                    shading='auto'
                )
                ax.set_xscale('log')
                ax.set_xlabel('Catalyticity γ (log scale)', fontsize=12)
                ax.set_ylabel('Wall Temperature T_w (K)', fontsize=12)
                ax.grid(True, alpha=0.15, which='both')
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label('Heat Flux q_wall (W/m²)', fontsize=12)
            else:
                ax.axis('off')
                ax.text(0.5, 0.5, f'Single T_w = {float(T[0]):.1f} K',
                        ha='center', va='center', fontsize=14)
            
        else:  # temperature_sweep mode
            fig.suptitle('Abacus: Heat Flux Map (T_w-sweep)', fontsize=14, fontweight='bold')
            
            if T.size >= 2:
                # 2D map
                im = ax.imshow(
                    Q,
                    extent=[T.min(), T.max(), G.min(), G.max()],
                    aspect='auto',
                    origin='lower',
                    cmap='Greys_r',
                    interpolation='bicubic'
                )
                ax.set_xlabel('Wall Temperature T_w (K)', fontsize=12)
                ax.set_ylabel('Catalyticity γ', fontsize=12)
                ax.grid(True, alpha=0.15)
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label('Heat Flux q_wall (W/m²)', fontsize=12)
            else:
                ax.axis('off')
                ax.text(0.5, 0.5, f'Single T_w = {float(T[0]):.1f} K',
                        ha='center', va='center', fontsize=14)

        plt.tight_layout()

        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            out_path = save_path.with_suffix(f'.{self.output_format}')
            plt.savefig(out_path, format=self.output_format, bbox_inches='tight', dpi=self.dpi)
        else:
            plt.show()

        plt.close()
        return True

    def plot_individual_curves(self, save_path: Optional[Path] = None) -> bool:
        """Plot only the curves of heat flux"""
        if not self.reader.data:
            self.reader.load()

        T = self.reader.data['temperatures']           # (N,)
        G = self.reader.data['catalyticity_values']    # (M,)
        Q = self.reader.data['heat_fluxes']            # (M, N)

        fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=self.dpi)
        
        if self.mode == 'gamma_sweep':
            fig.suptitle('Abacus: Heat Flux Curves (γ-sweep)', fontsize=14, fontweight='bold')
            
            # q_wall(gamma) curves for each T_w
            colors = ['black', 'red', 'blue', 'green', 'purple', 'orange', 'brown']
            markers = ['o', 's', '^', 'v', 'D', 'p', '*']
            line_styles = ['-', '--', '-.', ':']
            
            for j in range(min(10, T.size)):  # Limit to 10 curves for clarity
                color = colors[j % len(colors)]
                marker = markers[j % len(markers)]
                style = line_styles[j % len(line_styles)]
                ax.semilogx(
                    G, Q[:, j],
                    color=color,
                    linestyle=style,
                    marker=marker,
                    markevery=max(1, len(G)//15),
                    linewidth=2,
                    markersize=6,
                    label=f'T_w = {T[j]:.0f} K'
                )
            
            ax.set_xlabel('Catalyticity γ (log scale)', fontsize=12)
            ax.set_ylabel('Heat Flux q_wall (W/m²)', fontsize=12)
            ax.grid(True, alpha=0.3, which='both')
            ax.legend(loc='best')
            
        else:  # temperature_sweep mode
            fig.suptitle('Abacus: Heat Flux Curves (T_w-sweep)', fontsize=14, fontweight='bold')
            
            if T.size >= 2:
                # q_wall(T_w) curves for each gamma
                line_styles = ['-', '--', '-.', ':']
                markers = ['', 'o', 's', '^', 'v', 'D', 'p', '*', 'h', 'H', '+', 'x']

                for i in range(min(10, Q.shape[0])):  # Limit to 10 curves for clarity
                    style = line_styles[i % len(line_styles)]
                    marker = markers[i % len(markers)]
                    ax.plot(
                        T, Q[i, :],
                        color='black',
                        linestyle=style,
                        marker=(marker if marker else None),
                        markevery=max(1, len(T)//10),
                        linewidth=1.5,
                        markersize=3.5,
                        label=f'γ = {G[i]:g}'
                    )

                ax.set_xlabel('Wall Temperature T_w (K)', fontsize=12)
                ax.set_ylabel('Heat Flux q_wall (W/m²)', fontsize=12)
                ax.grid(True, alpha=0.3)
                ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
            else:
                # Single temperature
                ax.plot(G, Q[:, 0], color='black', marker='o')
                ax.set_xlabel('Catalyticity γ', fontsize=12)
                ax.set_ylabel('Heat Flux q_wall (W/m²)', fontsize=12)
                ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            out_path = save_path.with_suffix(f'.{self.output_format}')
            plt.savefig(out_path, format=self.output_format, bbox_inches='tight', dpi=self.dpi)
        else:
            plt.show()

        plt.close()
        return True

    def plot_individual_all(self, output_dir: Path) -> bool:
        """Generate all individual plots (map and curves separately)"""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        print(f"\nGenerating individual abacus plots in: {output_dir}")
        print("-" * 60)
        
        # Generate both plots
        plots = [
            ('abacus_map', self.plot_individual_map),
            ('abacus_curves', self.plot_individual_curves)
        ]
        
        success_count = 0
        for plot_name, plot_func in plots:
            save_path = output_dir / plot_name
            if plot_func(save_path):
                print(f"   [OK] {plot_name}: {save_path.name}.{self.output_format}")
                success_count += 1
            else:
                print(f"   [FAILED] {plot_name}")
        
        print("-" * 60)
        print(f"Individual abacus plots completed: {success_count}/{len(plots)} successful")
        
        return success_count == len(plots)

    def plot_map_and_curves(self, save_path: Optional[Path] = None) -> bool:
        if not self.reader.data:
            self.reader.load()

        T = self.reader.data['temperatures']           # (N,)
        G = self.reader.data['catalyticity_values']    # (M,)
        Q = self.reader.data['heat_fluxes']            # (M, N)

        # Figure with two subplots in one row
        fig, (ax_map, ax_curves) = plt.subplots(1, 2, figsize=(16, 7), dpi=self.dpi)
        
        if self.mode == 'gamma_sweep':
            # Gamma sweep mode: plot q(γ) for each T_w
            fig.suptitle('Abacus: γ-sweep at fixed T_w', fontsize=16, fontweight='bold')
            
            if T.size >= 2:
                # LEFT: 2D map with gamma on X-axis (log scale)
                # For imshow with log scale, we need to use pcolormesh instead
                import numpy as np
                G_grid, T_grid = np.meshgrid(G, T)
                im = ax_map.pcolormesh(
                    G_grid, T_grid, Q.T,
                    cmap='Greys_r',
                    shading='auto'
                )
                ax_map.set_xscale('log')  # Set log scale for gamma axis
                ax_map.set_xlabel('Catalyticity γ (log scale)')
                ax_map.set_ylabel('Wall Temperature T_w (K)')
                ax_map.grid(True, alpha=0.15, which='both')
                cbar = plt.colorbar(im, ax=ax_map)
                cbar.set_label('Heat Flux q_wall (W/m²)')
            else:
                # Single temperature
                ax_map.axis('off')
                ax_map.text(0.5, 0.5, f'Single T_w = {float(T[0]):.1f} K',
                            ha='center', va='center', fontsize=14)
            
            # RIGHT: q_wall(gamma) curves for each T_w
            colors = ['black', 'red', 'blue', 'green', 'purple', 'orange', 'brown']
            markers = ['o', 's', '^', 'v', 'D', 'p', '*']
            line_styles = ['-', '--', '-.', ':']
            
            for j in range(T.size):
                color = colors[j % len(colors)]
                marker = markers[j % len(markers)]
                style = line_styles[j % len(line_styles)]
                ax_curves.semilogx(  # Use semilogx for log scale on x-axis
                    G, Q[:, j],
                    color=color,
                    linestyle=style,
                    marker=marker,
                    markevery=max(1, len(G)//15),
                    linewidth=2,
                    markersize=6,
                    label=f'T_w = {T[j]:.0f} K'
                )
            
            ax_curves.set_xlabel('Catalyticity γ (log scale)')
            ax_curves.set_ylabel('Heat Flux q_wall (W/m²)')
            ax_curves.grid(True, alpha=0.3, which='both')  # Grid for both major and minor ticks
            ax_curves.legend(loc='best')
            
        else:  # temperature_sweep mode
            fig.suptitle('Abacus: T_w-sweep at various γ', fontsize=16, fontweight='bold')
            
            if T.size >= 2:
                # LEFT: 2D map (imshow) with explicit extent; axis 0 = gamma
                im = ax_map.imshow(
                    Q,
                    extent=[T.min(), T.max(), G.min(), G.max()],
                    aspect='auto',
                    origin='lower',
                    cmap='Greys_r',
                    interpolation='bicubic'
                )
                ax_map.set_xlabel('Wall Temperature T_w (K)')
                ax_map.set_ylabel('Catalyticity γ')
                ax_map.grid(True, alpha=0.15)
                cbar = plt.colorbar(im, ax=ax_map)
                cbar.set_label('Heat Flux q_wall (W/m²)')

                # RIGHT: q_wall(T_w) curves for each gamma (no extra interpolation)
                # Black-only, varied linestyles/markers to differentiate series
                line_styles = ['-', '--', '-.', ':']
                markers = ['', 'o', 's', '^', 'v', 'D', 'p', '*', 'h', 'H', '+', 'x']

                for i in range(min(10, Q.shape[0])):  # Limit to 10 curves for clarity
                    style = line_styles[i % len(line_styles)]
                    marker = markers[i % len(markers)]
                    ax_curves.plot(
                        T, Q[i, :],
                        color='black',
                        linestyle=style,
                        marker=(marker if marker else None),
                        markevery=max(1, len(T)//10),
                        linewidth=1.5,
                        markersize=3.5,
                        label=f'γ = {G[i]:g}'
                    )

                ax_curves.set_xlabel('Wall Temperature T_w (K)')
                ax_curves.set_ylabel('Heat Flux q_wall (W/m²)')
                ax_curves.grid(True, alpha=0.3)
                ax_curves.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
            else:
                # Single-temperature mode: show γ → q_wall curve on the right, annotate left
                ax_map.axis('off')
                ax_map.text(0.5, 0.5, f'Single T_w = {float(T[0]):.1f} K',
                            ha='center', va='center', fontsize=14)

                ax_curves.plot(G, Q[:, 0], color='black', marker='o')
                ax_curves.set_xlabel('Catalyticity γ')
                ax_curves.set_ylabel('Heat Flux q_wall (W/m²)')
                ax_curves.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            out_path = save_path.with_suffix(f'.{self.output_format}')
            plt.savefig(out_path, format=self.output_format, bbox_inches='tight', dpi=self.dpi)
            print(f"Abacus map + curves saved to: {out_path}")
        else:
            plt.show()

        plt.close()
        return True

def main() -> int:
    parser = argparse.ArgumentParser(description='BLAST Abacus Post-Processing (map + curves)')
    parser.add_argument('--input', '-i', type=str, required=True, help='Input abacus HDF5 file')
    parser.add_argument('--output', '-o', type=str, default='result_abacus/abacus_map', help='Output path (without extension) if not --show')
    parser.add_argument('--mode', '-m', type=str, choices=['gamma_sweep', 'temperature_sweep'], required=True, 
                        help='Plotting mode: gamma_sweep (q vs γ) or temperature_sweep (q vs T_w)')
    parser.add_argument('--plots', '-p', type=str, default='combined', 
                        choices=['combined', 'individual'], 
                        help='Type of plots to generate: combined (map + curves in one figure) or individual (separate PDFs)')
    parser.add_argument('--dpi', type=int, default=300, help='Figure DPI')
    parser.add_argument('--format', '-f', type=str, default='pdf', choices=['pdf', 'png', 'svg'], help='Output format for saved figure')
    parser.add_argument('--suppress-warnings', type=str, default='none', choices=['none', 'plot', 'all'], help='Control warnings: none (default), plot, or all')
    parser.add_argument('--show', action='store_true', help='Show interactively instead of saving')
    args = parser.parse_args()

    # Configure warnings
    def _configure_warnings(mode: str):
        warnings.resetwarnings()
        if mode == 'all':
            warnings.filterwarnings('ignore')
        elif mode == 'plot':
            warnings.filterwarnings('ignore', category=UserWarning, module=r'matplotlib')
            warnings.filterwarnings('ignore', category=UserWarning, module=r'seaborn')
            warnings.filterwarnings('ignore', category=FutureWarning, module=r'h5py')

    _configure_warnings(args.suppress_warnings)

    try:
        reader = AbacusReader(args.input)
        reader.load()
    except AbacusError as e:
        print(f"Error initializing reader: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error initializing reader: {e}")
        return 1

    plotter = AbacusPlotter(reader, output_format=args.format, dpi=args.dpi, mode=args.mode)

    try:
        if args.plots == 'individual':
            if args.show:
                # Show individual plots interactively (show map only)
                ok = plotter.plot_individual_map(save_path=None)
            else:
                # Generate separate PDF files for map and curves
                output_dir = Path(args.output).parent if Path(args.output).suffix else Path(args.output)
                ok = plotter.plot_individual_all(output_dir)
        else:  # combined
            if args.show:
                ok = plotter.plot_map_and_curves(save_path=None)
            else:
                ok = plotter.plot_map_and_curves(save_path=Path(args.output))
        
        if not ok:
            return 1
        print("Post-processing completed successfully!")
        return 0
    except Exception as e:
        print(f"Error during plotting: {e}")
        return 1

if __name__ == '__main__':
    exit(main())

"""
Usage examples:

# Generate combined plot (map + curves in one figure) - default
python3 postprocess_abacus.py --input abacus.h5 --mode temperature_sweep --output result_abacus/abacus_map

# Generate individual PDF files (map and curves separately)
python3 postprocess_abacus.py --input abacus.h5 --mode temperature_sweep --plots individual --output result_abacus/

# Gamma sweep mode with individual plots
python3 postprocess_abacus.py --input abacus.h5 --mode gamma_sweep --plots individual --output result_abacus/

# Interactive display (combined)
python3 postprocess_abacus.py --input abacus.h5 --mode temperature_sweep --show

# Interactive display (individual map only)
python3 postprocess_abacus.py --input abacus.h5 --mode temperature_sweep --plots individual --show
"""
