#!/usr/bin/env python3
"""
Script: FINITE SIZE COMPARISON WITH ANALYTICAL SOLUTION
Expected folders: "100", "400", "800", "1600"
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# --- CONFIGURATION ---
FOLDERS = ["100", "400", "800", "1600"]
LABELS  = ["100", "400", "800", "1600"]
COLORS = ['blue', 'green', 'red', "black"]
MARKERS = ['o', 's', '^', "<"]
ALPHA_SHADE = 0.1      # very light shadow for finite size error
SOL_COLOR = 'orange'    # color for analytical solution
SOL_LINEWIDTH = 2.5     # thicker line for analytical solution

# -------------------- FUNCTIONS --------------------
def read_data_file(filename, is_csv=False):
    """
    Read data file:
    - Simulation files (.dat): space/tab-delimited
    - Analytical solution (.dat CSV-like): comma-delimited
    Returns x, y, err (err=0 if not present)
    """
    if not os.path.exists(filename):
        print(f"Warning: {filename} not found.")
        return None, None, None
    try:
        if is_csv:
            # Skip comment lines starting with '#', take first 3 columns
            data = np.genfromtxt(filename, delimiter=',', comments='#', usecols=(0,1,2))
        else:
            # Space/tab-delimited simulation data
            data = np.loadtxt(filename, comments='#')
            if data.ndim == 1:
                data = data[np.newaxis, :]
        x = data[:, 0]
        y = data[:, 1]
        err = data[:, 2] if data.shape[1] >= 3 else np.zeros_like(y)
        return x, y, err
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None, None, None

def normalize_x(x):
    """Normalize x to range 0-1"""
    if x is None: return None
    return (x - x.min()) / (x.max() - x.min())

def subsample(x, y, err, n=25):
    """Reduce number of points for clear visualization"""
    if x is None: return None, None, None
    if len(x) <= n: return x, y, err
    idx = np.linspace(0, len(x)-1, n, dtype=int)
    return x[idx], y[idx], err[idx]

# -------------------- MAIN --------------------
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Load analytical solution (CSV)
    x_sol, y_sol, _ = read_data_file('solution.dat', is_csv=True)
    if x_sol is not None:
        x_sol = normalize_x(x_sol)
    
    print("--- Generating Finite Size Comparison ---")
    
    # ==========================================================
    # ==================== DENSITY FIGURE ======================
    # ==========================================================
    fig_rho, ax_rho = plt.subplots(figsize=(7, 6))
    
    for folder, label, color, mark in zip(FOLDERS, LABELS, COLORS, MARKERS):
        path_rho = os.path.join(folder, 'density_profile.dat')
        
        x, y, err = read_data_file(path_rho)
        if x is not None:
            x = normalize_x(x)
            xs, ys, es = subsample(x, y, err)
            
            ax_rho.fill_between(x, y - err, y + err,
                                color=color, alpha=ALPHA_SHADE, zorder=1)
            
            ax_rho.plot(x, y, color=color,
                        alpha=0.6, linewidth=1, zorder=2)
            
            ax_rho.errorbar(xs, ys, yerr=es,
                            fmt=mark, color=color, label=label,
                            capsize=3, markersize=5,
                            alpha=0.9, zorder=3)
    
    if x_sol is not None:
        ax_rho.plot(x_sol, y_sol,
                    color=SOL_COLOR, linestyle='--',
                    linewidth=SOL_LINEWIDTH,
                    label='Analytical', zorder=4)
    
    ax_rho.set_title("Density Profile ($\\rho$)", fontsize=14)
    ax_rho.set_xlabel("Position x/L")
    ax_rho.set_ylabel(r"$\rho(x)$")
    ax_rho.grid(True, alpha=0.3)
    ax_rho.legend()
    
    plt.tight_layout()
    plt.savefig('density_comparison.png', dpi=300)
    print("Saved: density_comparison.png")
    plt.show()
    
    
    # ==========================================================
    # ================== TEMPERATURE FIGURE ====================
    # ==========================================================
    fig_temp, ax_temp = plt.subplots(figsize=(7, 6))
    
    for folder, label, color, mark in zip(FOLDERS, LABELS, COLORS, MARKERS):
        path_temp = os.path.join(folder, 'temperature_profile.dat')
        
        x, y, err = read_data_file(path_temp)
        if x is not None:
            x = normalize_x(x)
            xs, ys, es = subsample(x, y, err)
            
            ax_temp.fill_between(x, y - err, y + err,
                                 color=color, alpha=ALPHA_SHADE, zorder=1)
            
            ax_temp.plot(x, y, color=color,
                         alpha=0.6, linewidth=1, zorder=2)
            
            ax_temp.errorbar(xs, ys, yerr=es,
                             fmt=mark, color=color, label=label,
                             capsize=3, markersize=5,
                             alpha=0.9, zorder=3)
    
    if x_sol is not None:
        ax_temp.plot(x_sol, y_sol,
                     color=SOL_COLOR, linestyle='--',
                     linewidth=SOL_LINEWIDTH,
                     label='Analytical', zorder=4)
    
    ax_temp.set_title("Temperature Profile ($T$)", fontsize=14)
    ax_temp.set_xlabel("Position x/L")
    ax_temp.set_ylabel(r"$T(x)$")
    ax_temp.grid(True, alpha=0.3)
    ax_temp.legend()
    
    plt.tight_layout()
    plt.savefig('temperature_comparison.png', dpi=300)
    print("Saved: temperature_comparison.png")
    plt.show()

# -------------------- RUN --------------------
if __name__ == "__main__":
    main()
