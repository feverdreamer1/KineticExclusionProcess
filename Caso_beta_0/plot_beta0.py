

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

base_dir =os.path.dirname(os.path.abspath(__file__))

cases = [
    ("N=100 T=2,1 b=0, rho=0.4 , 0.8", 2.0, 1.0, 0.4, 0.8, "T=2.0, ρ=0.4→0.8"),
    ("N=100 T=4,1 b=0,rho= 0.3 0.8", 4.0, 1.0, 0.3, 0.8, "T=4.0, ρ=0.3→0.8"),
    ("N=100 T=6,1 b=0 rho=0.4 0.6", 6.0, 1.0, 0.4, 0.6, "T=6.0, ρ=0.4→0.6"),
]

def read_theoretical_data(folder_path):
    """Read theoretical data from solution.csv"""
    file_path = os.path.join(folder_path, "solution.csv")
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found")
        return None, None, None
    
    try:
        data = np.loadtxt(file_path, delimiter=',', skiprows=1)
        x_theo = data[:, 0]  # Position (already normalized)
        rho_theo = data[:, 1]  # Theoretical density
        T_theo = data[:, 2]  # Theoretical temperature
        
        return x_theo, rho_theo, T_theo
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None, None, None

def read_density_file(folder_path):
    """Read density profile from density_profile.dat"""
    file_path = os.path.join(folder_path, "density_profile.dat")
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found")
        return None, None, None
    
    try:
        data = np.loadtxt(file_path, skiprows=1)
        if data.ndim == 1:
            positions = np.arange(1, len(data) + 1)
            density = data
            errors = None
        else:
            positions = data[:, 0]
            density = data[:, 1]
            if data.shape[1] > 2:
                errors = data[:, 2]
            else:
                errors = None
        
        # Normalize positions to [0, 1]
        L = len(positions)
        x_normalized = (positions - 1) / (L - 1) if L > 1 else positions
        
        return x_normalized, density, errors
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None, None, None

def read_temperature_file(folder_path):
    """Read temperature profile from temperature_profile.dat"""
    file_path = os.path.join(folder_path, "temperature_profile.dat")
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found")
        return None, None
    
    try:
        data = np.loadtxt(file_path, skiprows=1)
        if data.ndim == 1:
            positions = np.arange(1, len(data) + 1)
            temperature = data
        else:
            positions = data[:, 0]
            temperature = data[:, 1]
        
        # Normalize positions to [0, 1]
        L = len(positions)
        x_normalized = (positions - 1) / (L - 1) if L > 1 else positions
        
        return x_normalized, temperature
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None, None

def plot_density_comparison():
    """Plot density profiles for all cases with theoretical comparison"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    
    colors = ['blue', 'green', 'red']
    markers = ['o', 's', '^']  # circle, square, triangle
    
    for idx, (folder, T_left, T_right, rho_left, rho_right, label) in enumerate(cases):
        folder_path = os.path.join(base_dir, folder)
        
        x, density, errors = read_density_file(folder_path)
        x_theo, rho_theo, _ = read_theoretical_data(folder_path)
        
        if x is not None and density is not None:
            # Plot simulation data - only markers, no line
            if errors is not None:
                ax.errorbar(x[::12], density[::12], yerr=errors[::12], fmt=markers[idx], color=colors[idx], 
                           markersize=12, capsize=2, capthick=1, elinewidth=1,
                           label=f'Sim: {label}', alpha=0.7, linestyle='None')
            else:
                ax.plot(x[::12], density[::12], marker=markers[idx], color=colors[idx], 
                       markersize=12, linestyle='None', label=f'Sim: {label}', alpha=0.7)
            
            # Plot theoretical profile - only line, no markers
            if x_theo is not None and rho_theo is not None:
                ax.plot(x_theo, rho_theo, '-', color=colors[idx], linewidth=2, 
                       label=f'Theoretical: {label}', alpha=0.8)
    
    ax.set_xlabel('Normalized position (x/L)', fontsize=12)
    ax.set_ylabel('Density ρ', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc='best')
    ax.set_xlim(0, 1)
    
    plt.tight_layout()
    plt.savefig('density_comparison_beta0.png', dpi=300, bbox_inches='tight')
    print("Saved: density_comparison_beta0.png")
    return fig

def plot_temperature_comparison():
    """Plot temperature profiles for all cases with theoretical comparison"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    
    colors = ['blue', 'green', 'red']
    markers = ['o', 's', '^']  # circle, square, triangle
    
    for idx, (folder, T_left, T_right, rho_left, rho_right, label) in enumerate(cases):
        folder_path = os.path.join(base_dir, folder)
        
        x, temperature = read_temperature_file(folder_path)
        x_theo, _, T_theo = read_theoretical_data(folder_path)
        
        if x is not None and temperature is not None:
            # Plot simulation data - only markers, no line
            ax.plot(x[::12], temperature[::12], marker=markers[idx], color=colors[idx], 
                   markersize=12, linestyle='None', label=f'Sim: {label}', alpha=0.7)
            
            # Plot theoretical profile - only line, no markers
            if x_theo is not None and T_theo is not None:
                ax.plot(x_theo, T_theo, '-', color=colors[idx], linewidth=2, 
                       label=f'Theoretical: {label}', alpha=0.8)
    
    ax.set_xlabel('Normalized position (x/L)', fontsize=12)
    ax.set_ylabel('Temperature T', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc='best')
    ax.set_xlim(0, 1)
    
    plt.tight_layout()
    plt.savefig('temperature_comparison_beta0.png', dpi=300, bbox_inches='tight')
    print("Saved: temperature_comparison_beta0.png")
    return fig

def main():
    """Main function"""
    print("=" * 60)
    print("Generando gráficas para β=0")
    print("=" * 60)
    
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    print(f"Directorio de trabajo: {os.getcwd()}")
    print(f"Directorio base de datos: {base_dir}\n")
    
    # Check if data folders exist
    print("Verificando carpetas de datos...")
    for folder, _, _, _, _, label in cases:
        folder_path = os.path.join(base_dir, folder)
        if os.path.exists(folder_path):
            print(f"  ✓ {folder}")
        else:
            print(f"  ✗ {folder} NO ENCONTRADA")
    
    print("\nGenerando gráfica de densidad...")
    plot_density_comparison()
    
    print("Generando gráfica de temperatura...")
    plot_temperature_comparison()
    
    print("\n" + "=" * 60)
    print("¡Gráficas generadas exitosamente!")
    print("=" * 60)
    print("\nArchivos generados:")
    print("  - density_comparison_beta0.png")
    print("  - temperature_comparison_beta0.png")
    
    # Show plots
    plt.show()

if __name__ == "__main__":
    main()
