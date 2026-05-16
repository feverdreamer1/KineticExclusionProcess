#!/usr/bin/env python3
"""
Script 2: COMPARACION DE GRADIENTE TERMICO
Carpetas esperadas: "T=2", "T=5", "T=10"
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# --- CONFIGURACION ---
FOLDERS = ["T=2", "T=5", "T=10"]
LABELS = [r"$T_L=2$", r"$T_L=5$", r"$T_L=10$"]
# Usamos mapa de color 'plasma' o 'coolwarm' manual
COLORS = ['#3498db', '#e67e22', '#c0392b'] # Azul, Naranja, Rojo Oscuro
MARKERS = ['o', 'D', 'v']

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
    
    print("--- Generando Comparativa de Gradiente T ---")
    
    for folder, label, color, mark in zip(FOLDERS, LABELS, COLORS, MARKERS):
        path_rho = os.path.join(folder, 'density_profile.dat')
        path_temp = os.path.join(folder, 'temperature_profile.dat')
        
        # 1. DENSIDAD (Buscamos Efecto Soret)
        x, y, err = read_data_file(path_rho)
        if x is not None:
            x = normalize_x(x)
            xs, ys, es = subsample(x, y, err)
            ax_rho.plot(x, y, color=color, alpha=0.4)
            ax_rho.errorbar(xs, ys, yerr=es, fmt=mark, color=color, label=label,
                            capsize=3, markersize=5)
        
        # 2. TEMPERATURA
        x, y, err = read_data_file(path_temp)
        if x is not None:
            x = normalize_x(x)
            xs, ys, es = subsample(x, y, err)
            ax_temp.plot(x, y, color=color, alpha=0.4)
            ax_temp.errorbar(xs, ys, yerr=es, fmt=mark, color=color, label=label,
                             capsize=3, markersize=5)

    ax_rho.set_xlabel("Normalized position x/L")
    ax_rho.set_ylabel(r"Density $\rho$")
    ax_rho.grid(True, alpha=0.3)
    ax_rho.legend()
    
    ax_temp.set_title("", fontsize=14)
    ax_temp.set_xlabel("Normalized position x/L")
    ax_temp.set_ylabel(r"Temperature $T$")
    ax_temp.grid(True, alpha=0.3)
    ax_temp.legend()
    
    plt.tight_layout()
    plt.savefig('comparacion_thermal.png', dpi=300)
    print("Guardado: comparacion_thermal.png")
    plt.show()

if __name__ == "__main__":
    main()