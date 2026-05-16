#!/usr/bin/env python3
"""
SSEP Validation Plotting Script
Compares KEP simulation (beta=0) with theoretical SSEP solution
Theoretical SSEP: rho(x) = rho_L + (x/(L+1)) * (rho_R - rho_L)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters (should match ssep_sim.f90)
L = 200
rho_left = 0.8
rho_right = 0.2

def theoretical_ssep_density(x):
    """
    Theoretical SSEP density profile
    rho(x) = rho_L + (x/(L+1)) * (rho_R - rho_L)
    
    Parameters:
    -----------
    x : array-like
        Position normalized to [0, 1] (x/L)
    
    Returns:
    --------
    rho : array-like
        Theoretical density at each position
    """
    # Convert normalized position to actual position in [0, L+1]
    x_actual = x * (L + 1)
    # Apply theoretical formula
    rho = rho_left + (x_actual / (L + 1)) * (rho_right - rho_left)
    return rho

def read_density_file(filename='density_profile.dat'):
    """
    Read density profile from Fortran output
    
    Parameters:
    -----------
    filename : str
        Name of the density profile file
    
    Returns:
    --------
    x : numpy array
        Normalized positions [0, 1]
    rho : numpy array
        Density values
    rho_err : numpy array or None
        Density errors (if available)
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
        x_normalized = (positions - 1) / (L - 1) if L > 1 else positions
        
        return x_normalized, density, errors
    
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None, None, None

def plot_validation():
    """
    Create validation plot comparing simulation with theoretical SSEP solution
    """
    # Read simulation data
    x_sim, rho_sim, rho_err = read_density_file('density_profile.dat')
    
    if x_sim is None or rho_sim is None:
        print("Could not read simulation data. Make sure density_profile.dat exists.")
        return
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    
    # Plot simulation data with error bars if available
    if rho_err is not None:
        ax.errorbar(x_sim[::5], rho_sim[::5], yerr=rho_err[::5], fmt='o', color='blue', 
                   markersize=6, capsize=2, capthick=1, elinewidth=1,
                   label='Simulation', alpha=0.7)
    else:
        ax.plot(x_sim[::5], rho_sim[::5], 'o', color='blue', markersize=6,
               label='Simulation', alpha=0.7)
    
    # Plot theoretical SSEP solution
    x_theo = np.linspace(0, 1, 1000)
    rho_theo = theoretical_ssep_density(x_theo)
    ax.plot(x_theo, rho_theo, '-', color='red', linewidth=2,
           label='SSEP theory', alpha=0.8)
    
    # Labels and formatting
    ax.set_xlabel('Normalized position (x/L)', fontsize=12)
    ax.set_ylabel(r'Density $\rho$', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    ax.set_xlim(0, 1)
    
    # # Add text box with parameters (Commented out to match KMP script)
    # textstr = f'Parameters:\nrho_L = {rho_left}\nrho_R = {rho_right}\nL = {L}'
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
    #         verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig('validation_ssep.png', dpi=300, bbox_inches='tight')
    print("Plot saved: validation_ssep.png")
    
    # Show plot
    plt.show()

def main():
    """Main function"""
    print("=" * 60)
    print("SSEP Validation Plotting")
    print("=" * 60)
    
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    print(f"Working directory: {os.getcwd()}\n")
    
    # Check if data file exists
    if not os.path.exists('density_profile.dat'):
        print("WARNING: density_profile.dat not found.")
        print("Please run ssep_sim.f90 first to generate the data.")
        return
    
    # Create plot
    plot_validation()
    
    print("\n" + "=" * 60)
    print("Validation completed!")
    print("=" * 60)

if __name__ == "__main__":
    main()