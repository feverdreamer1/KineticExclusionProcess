#!/usr/bin/env python3
"""
KMP Validation Plotting Script
Compares KEP simulation (beta=0) with theoretical KMP solution
Theoretical KMP: T(x) = T_L + (x/(L+1)) * (T_R - T_L)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters (should match kmp_sim.f90)
L = 200
T_left = 4.0
T_right = 1.0

def theoretical_kmp_temperature(x):
    """
    Theoretical KMP temperature profile
    T(x) = T_L + (x/(L+1)) * (T_R - T_L)
    
    Parameters:
    -----------
    x : array-like
        Position normalized to [0, 1] (x/L)
    
    Returns:
    --------
    T : array-like
        Theoretical temperature at each position
    """
    # Convert normalized position to actual position in [0, L+1]
    x_actual = x * (L + 1)
    # Apply theoretical formula
    T = T_left + (x_actual / (L + 1)) * (T_right - T_left)
    return T

def read_temperature_file(filename='temperature_profile.dat'):
    """
    Read temperature profile from Fortran output
    
    Parameters:
    -----------
    filename : str
        Name of the temperature profile file
    
    Returns:
    --------
    x : numpy array
        Normalized positions [0, 1]
    T : numpy array
        Temperature values
    T_err : numpy array or None
        Temperature errors (if available)
    """
    if not os.path.exists(filename):
        print(f"Error: {filename} not found!")
        return None, None, None
    
    try:
        # Read data (skip header line)
        data = np.loadtxt(filename, skiprows=1)
        
        if data.ndim == 1:
            # Single column (shouldn't happen with our format)
            positions = np.arange(1, len(data) + 1)
            temperature = data
            errors = None
        else:
            positions = data[:, 0]
            temperature = data[:, 1]
            if data.shape[1] > 2:
                errors = data[:, 2]
            else:
                errors = None
        
        # Normalize positions to [0, 1]
        x_normalized = (positions - 1) / (L - 1) if L > 1 else positions
        
        return x_normalized, temperature, errors
    
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None, None, None

def plot_validation():
    """
    Create validation plot comparing simulation with theoretical KMP solution
    """
    # Read simulation data
    x_sim, T_sim, T_err = read_temperature_file('temperature_profile.dat')
    
    if x_sim is None or T_sim is None:
        print("Could not read simulation data. Make sure temperature_profile.dat exists.")
        return
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    
    # Plot simulation data with error bars if available
    if T_err is not None:
        ax.errorbar(x_sim[::5], T_sim[::5], yerr=T_err, fmt='o', color='blue', 
                   markersize=6, capsize=2, capthick=1, elinewidth=1,
                   label='Simulation', alpha=0.7)
    else:
        ax.plot(x_sim[::5], T_sim[::5], 'o', color='blue', markersize=6,
               label='Simulation', alpha=0.7)
    
    # Plot theoretical KMP solution
    x_theo = np.linspace(0, 1, 1000)
    T_theo = theoretical_kmp_temperature(x_theo)
    ax.plot(x_theo, T_theo, '-', color='red', linewidth=2,
           label='KMP theory', alpha=0.8)
    
    # Labels and formatting
    ax.set_xlabel('Normalized position (x/L)', fontsize=12)
    ax.set_ylabel('Temperature T', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    ax.set_xlim(0, 1)
    
    # # Add text box with parameters
    # textstr = f'p:\nT_L = {T_left}\nT_R = {T_right}\nL = {L}'
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
    #        verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig('validation_kmp.png', dpi=300, bbox_inches='tight')
    print("Gráfica guardada: validation_kmp.png")
    
    # Show plot
    plt.show()

def main():
    """Main function"""
    print("=" * 60)
    print("KMP Validation Plotting")
    print("=" * 60)
    
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    print(f"Directorio de trabajo: {os.getcwd()}\n")
    
    # Check if data file exists
    if not os.path.exists('temperature_profile.dat'):
        print("ADVERTENCIA: temperature_profile.dat no encontrado.")
        print("Por favor ejecuta primero kmp_sim.f90 para generar los datos.")
        return
    
    # Create plot
    plot_validation()
    
    print("\n" + "=" * 60)
    print("Validación completada!")
    print("=" * 60)

if __name__ == "__main__":
    main()
