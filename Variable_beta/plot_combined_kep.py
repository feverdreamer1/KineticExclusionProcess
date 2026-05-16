#!/usr/bin/env python3
"""
Script 1: COMPARACION DE EFECTO BETA
Carpetas esperadas: "beta=0", "beta=0.5", "beta=1"
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# --- CONFIGURACION ---
FOLDERS = ["beta=0", "beta=0.5", "beta=1"]
LABELS = [r"$\beta=0$ ", r"$\beta=0.5$ ", r"$\beta=1$ "]
COLORS = ['blue', 'green', 'red']
MARKERS = ['o', 's', '^']

def read_data_file(filename):
    """Lee datos: Posicion, Valor, Error"""
    try:
        if not os.path.exists(filename):
            print(f"Warning: {filename} no encontrado.")
            return None, None, None
        data = np.loadtxt(filename, skiprows=1)
        if data.ndim == 1 or data.shape[1] < 2: return None, None, None
        x = data[:, 0]
        y = data[:, 1]
        err = data[:, 2] if data.shape[1] >= 3 else np.zeros_like(y)
        return x, y, err
    except Exception as e:
        print(f"Error leyendo {filename}: {e}")
        return None, None, None

def normalize_x(x):
    if x is None: return None
    return (x - x.min()) / (x.max() - x.min())

def subsample(x, y, err, n=25):
    """Reduce puntos para visualizacion clara"""
    if x is None: return None, None, None
    if len(x) <= n: return x, y, err
    idx = np.linspace(0, len(x)-1, n, dtype=int)
    return x[idx], y[idx], err[idx]

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Preparamos figura (2 paneles: Densidad y Temperatura)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax_rho, ax_temp = axes[0], axes[1]
    
    print("--- Generando Comparativa de BETA ---")
    
    for folder, label, color, mark in zip(FOLDERS, LABELS, COLORS, MARKERS):
        path_rho = os.path.join(folder, 'density_profile.dat')
        path_temp = os.path.join(folder, 'temperature_profile.dat')
        
        # 1. DENSIDAD
        x, y, err = read_data_file(path_rho)
        if x is not None:
            x = normalize_x(x)
            xs, ys, es = subsample(x, y, err)
            # Linea suave completa
            ax_rho.plot(x, y, color=color, alpha=0.4, linewidth=1)
            # Puntos con error
            ax_rho.errorbar(xs, ys, yerr=es, fmt=mark, color=color, label=label,
                            capsize=3, markersize=5, alpha=0.9)
        
        # 2. TEMPERATURA
        x, y, err = read_data_file(path_temp)
        if x is not None:
            x = normalize_x(x)
            xs, ys, es = subsample(x, y, err)
            ax_temp.plot(x, y, color=color, alpha=0.4, linewidth=1)
            ax_temp.errorbar(xs, ys, yerr=es, fmt=mark, color=color, label=label,
                             capsize=3, markersize=5, alpha=0.9)

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
    plt.savefig('comparacion_beta.png', dpi=300)
    print("Guardado: comparacion_beta.png")
    plt.show()

if __name__ == "__main__":
    main()