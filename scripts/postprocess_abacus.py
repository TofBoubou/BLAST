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
warnings.filterwarnings('ignore')

# Match existing plotting style
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

class AbacusFileError(Exception):
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
            raise FileNotFoundError(f"Input file not found: {self.input_path}")
        if self.input_path.suffix.lower() not in ('.h5', '.hdf5'):
            raise ValueError(f"Only HDF5 is supported. Got: {self.input_path}")

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
                    raise AbacusFileError("'heat_fluxes' must be 2D [M, N]")
                M, N = q.shape
                if temps.ndim != 1 or temps.size != N:
                    raise AbacusFileError("'temperatures'.size must equal heat_fluxes.shape[1]")
                if gammas.ndim != 1 or gammas.size != M:
                    raise AbacusFileError("'catalyticity_values'.size must equal heat_fluxes.shape[0]")

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
            if isinstance(e, AbacusFileError):
                raise
            raise RuntimeError(f"Failed to load HDF5: {e}")

class AbacusPlotter:
    """One figure with two subplots: map (left) + curves (right)"""

    def __init__(self, reader: AbacusReader, dpi: int = 300):
        self.reader = reader
        self.dpi = dpi

    def plot_map_and_curves(self, save_path: Optional[Path] = None) -> bool:
        if not self.reader.data:
            self.reader.load()

        T = self.reader.data['temperatures']           # (N,)
        G = self.reader.data['catalyticity_values']    # (M,)
        Q = self.reader.data['heat_fluxes']            # (M, N)

        # Figure with two subplots in one row
        fig, (ax_map, ax_curves) = plt.subplots(1, 2, figsize=(16, 7), dpi=self.dpi)
        fig.suptitle('Abacus: q_wall(γ, T_w)', fontsize=16, fontweight='bold')

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

            for i in range(Q.shape[0]):
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
            out_pdf = save_path.with_suffix('.pdf')
            plt.savefig(out_pdf, format='pdf', bbox_inches='tight')
            print(f"Abacus map + curves saved to: {out_pdf}")
        else:
            plt.show()

        plt.close()
        return True

def main() -> int:
    parser = argparse.ArgumentParser(description='BLAST Abacus Post-Processing (map + curves)')
    parser.add_argument('--input', '-i', type=str, required=True, help='Input abacus HDF5 file')
    parser.add_argument('--output', '-o', type=str, default='result_abacus/abacus_map', help='Output path (without extension) if not --show')
    parser.add_argument('--dpi', type=int, default=300, help='Figure DPI')
    parser.add_argument('--show', action='store_true', help='Show interactively instead of saving')
    args = parser.parse_args()

    try:
        reader = AbacusReader(args.input)
        reader.load()
    except Exception as e:
        print(f"Error initializing reader: {e}")
        return 1

    plotter = AbacusPlotter(reader, dpi=args.dpi)

    try:
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
