#!/usr/bin/env python3
"""
Script 3: COMPARACION DE GRADIENTE DE DENSIDAD
Carpetas esperadas: "rho=0.2", "rho=0.5", "rho=0.9"
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# --- CONFIGURACION ---
FOLDERS = ["rho=0.2", "rho=0.5", "rho=0.9"]
LABELS = [r"$\rho_L=0.2$ ", r"$\rho_L=0.5$", r"$\rho_L=0.9$ "]
COLORS = ['#88CCEE', '#3377AA', '#114477'] # Escala de Azules
MARKERS = ['o', 's', 'D']

def read_data_file(filename):
    try:
        if not os.path.exists(filename): return None, None, None
        data = np.loadtxt(filename, skiprows=1)
        if data.ndim == 1 or data.shape[1] < 2: return None, None, None
        x = data[:, 0]; y = data[:, 1]
        err = data[:, 2] if data.shape[1] >= 3 else np.zeros_like(y)
        return x, y, err
    except: return None, None, None

def normalize_x(x): return (x - x.min()) / (x.max() - x.min()) if x is not None else None

def subsample(x, y, err, n=25):
    if x is None: return None, None, None
    if len(x) <= n: return x, y, err
    idx = np.linspace(0, len(x)-1, n, dtype=int)
    return x[idx], y[idx], err[idx]

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax_rho, ax_temp = axes[0], axes[1]
    
    print("--- Generando Comparativa de Densidad ---")
    
    for folder, label, color, mark in zip(FOLDERS, LABELS, COLORS, MARKERS):
        path_rho = os.path.join(folder, 'density_profile.dat')
        path_temp = os.path.join(folder, 'temperature_profile.dat')
        
        # 1. DENSIDAD
        x, y, err = read_data_file(path_rho)
        if x is not None:
            x = normalize_x(x)
            xs, ys, es = subsample(x, y, err)
            
            # --- SOLO PUNTOS (Lineas eliminadas) ---
            ax_rho.errorbar(xs, ys, yerr=es, fmt=mark, color=color, label=label,
                            capsize=3, markersize=5, linestyle='none') 
            # 'linestyle=none' asegura que no haya lineas uniendo los marcadores
        
        # 2. TEMPERATURA (Buscamos Efecto Dufour)
        x, y, err = read_data_file(path_temp)
        if x is not None:
            x = normalize_x(x)
            xs, ys, es = subsample(x, y, err)
            
            # --- SOLO PUNTOS (Lineas eliminadas) ---
            ax_temp.errorbar(xs, ys, yerr=es, fmt=mark, color=color, label=label,
                             capsize=3, markersize=5, linestyle='none')
    
    # Zoom en temperatura si quieres ver mejor el efecto Dufour (Opcional)
    # ax_temp.set_ylim(1.95, 2.05) 

    ax_rho.set_xlabel("Normalized position x/L")
    ax_rho.set_ylabel(r"Density $\rho$")
    ax_rho.grid(True, alpha=0.3)
    ax_rho.legend()
    
    ax_temp.set_title("", fontsize=14)
    ax_temp.set_xlabel("Normalized position x/L")
    ax_temp.set_ylabel(r"Temperature $T$")
    ax_temp.set_ylim(1.5, 2.5) 
    ax_temp.grid(True, alpha=0.3)
    ax_temp.legend()
    
    plt.tight_layout()
    plt.savefig('comparacion_density.png', dpi=300)
    print("Guardado: comparacion_density.png")
    plt.show()

if __name__ == "__main__":
    main()