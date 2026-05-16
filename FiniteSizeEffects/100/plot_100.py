#!/usr/bin/env python3
"""
Combined KEP Plotting Script
Plots both simulation results (with error bars) and numerical algorithm solution on the same graphs
Parameters: T_left=6.0, T_right=1.0, beta=1.0, rho_left=0.4, rho_right=0.8
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def read_data_file(filename):
    """Read data from a file and return position, values, and errors (if available)"""
    try:
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found in current directory: {os.getcwd()}")
            return None, None, None
        
        data = np.loadtxt(filename, skiprows=1)  # Skip header line
        if data.ndim == 1:
            # Single column data (old format, no errors)
            return np.arange(1, len(data) + 1), data, None
        elif data.shape[1] == 2:
            # Two columns: position and value (old format, no errors)
            return data[:, 0], data[:, 1], None
        elif data.shape[1] >= 3:
            # Three or more columns: position, value, error
            return data[:, 0], data[:, 1], data[:, 2]
        else:
            return None, None, None
    except (FileNotFoundError, OSError, IOError) as e:
        print(f"Warning: {filename} not found ({type(e).__name__})")
        return None, None, None
    except Exception as e:
        print(f"Error reading {filename}: {type(e).__name__}: {e}")
        return None, None, None

def read_numerical_solution(filename):
    """Read numerical algorithm solution from solution.dat"""
    try:
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found in current directory: {os.getcwd()}")
            return None, None, None
        
        # Read CSV format, skip header lines (first 3 lines are comments)
        data = np.loadtxt(filename, delimiter=',', skiprows=3, comments='#')
        x = data[:, 0]
        rho = data[:, 1]
        T = data[:, 2]
        return x, rho, T
    except (FileNotFoundError, OSError, IOError) as e:
        print(f"Warning: {filename} not found ({type(e).__name__})")
        return None, None, None
    except Exception as e:
        print(f"Error reading {filename}: {type(e).__name__}: {e}")
        return None, None, None

def normalize_x_axis(positions):
    """Normalize positions to [0, 1] range"""
    if positions is None or len(positions) == 0:
        return None
    pos_min = np.min(positions)
    pos_max = np.max(positions)
    if pos_max == pos_min:
        return np.zeros_like(positions)
    return (positions - pos_min) / (pos_max - pos_min)

def subsample_data(x, y, yerr=None, n_points=25):
    """Subsample data to show only n_points evenly spaced points"""
    if x is None or y is None:
        return None, None, None
    if len(x) <= n_points:
        return x, y, yerr
    
    # Create evenly spaced indices
    indices = np.linspace(0, len(x) - 1, n_points, dtype=int)
    if yerr is not None:
        return x[indices], y[indices], yerr[indices]
    else:
        return x[indices], y[indices], None

def create_combined_plots():
    """Create combined plots showing simulation results and numerical algorithm solution"""
    
    # Get script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    print("=" * 60)
    print("KEP Combined Plotting (β ≠ 0)")
    print("=" * 60)
    
    # Read simulation data (with errors)
    print("\nReading simulation data...")
    pos_density, density_data, density_error = read_data_file('density_profile.dat')
    pos_temp, temp_data, temp_error = read_data_file('temperature_profile.dat')
    
    # Read numerical algorithm solution
    print("Reading numerical algorithm solution...")
    x_num, rho_num, T_num = read_numerical_solution('solution.dat')
    
    # Normalize simulation x-axis to [0, 1]
    if pos_density is not None:
        x_density = normalize_x_axis(pos_density)
    else:
        x_density = None
    
    if pos_temp is not None:
        x_temp = normalize_x_axis(pos_temp)
    else:
        x_temp = None
    
    # Subsample simulation data to show only 25 points (with errors)
    x_density_sub, density_data_sub, density_error_sub = subsample_data(x_density, density_data, density_error, n_points=25)
    x_temp_sub, temp_data_sub, temp_error_sub = subsample_data(x_temp, temp_data, temp_error, n_points=20)
    
    # Create figure with 2 subplots side by side
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('KEP Simulation Results vs Numerical Algorithm Solution (β ≠ 0)', 
                 fontsize=16, fontweight='bold')
    
    # Plot 1: Density Profile
    ax1 = axes[0]
    
    # Plot numerical algorithm solution
    if x_num is not None and rho_num is not None:
        ax1.plot(x_num, rho_num, 'r-', linewidth=2.5, label='Numerical Algorithm', alpha=0.8)
    
    # Plot simulation data (subsampled) with error bars
    if x_density_sub is not None and density_data_sub is not None:
        if density_error_sub is not None:
            ax1.errorbar(x_density_sub, density_data_sub, yerr=density_error_sub, 
                        fmt='bo-', linewidth=1.5, markersize=5, capsize=3, capthick=1.5,
                        label='Simulation Data', alpha=0.7, markeredgewidth=0.5, elinewidth=1.5)
        else:
            ax1.plot(x_density_sub, density_data_sub, 'bo-', linewidth=1.5, markersize=5, 
                    label='Simulation Data', alpha=0.7, markeredgewidth=0.5)
    
    ax1.set_xlabel('Normalized Position (x/L)', fontsize=13)
    ax1.set_ylabel('Density ρ', fontsize=13)
    ax1.set_title('Density Profile', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    ax1.legend(fontsize=11, loc='best')
    
    # Add statistics if simulation data available
    if density_data is not None:
        mean_density = np.mean(density_data)
        ax1.axhline(y=mean_density, color='g', linestyle='--', alpha=0.5, linewidth=1,
                   label=f'Mean (Sim): {mean_density:.3f}')
        ax1.legend(fontsize=11, loc='best')
    
    # Plot 2: Temperature Profile
    ax2 = axes[1]
    
    # Plot numerical algorithm solution
    if x_num is not None and T_num is not None:
        ax2.plot(x_num, T_num, 'r-', linewidth=2.5, label='Numerical Algorithm', alpha=0.8)
    
    # Plot simulation data (subsampled) with error bars
    if x_temp_sub is not None and temp_data_sub is not None:
        if temp_error_sub is not None:
            ax2.errorbar(x_temp_sub, temp_data_sub, yerr=temp_error_sub, 
                        fmt='bo-', linewidth=1.5, markersize=8, capsize=3, capthick=1.5,
                        label='Simulation Data', alpha=0.7, markeredgewidth=0.5, elinewidth=1.5)
        else:
            ax2.plot(x_temp_sub, temp_data_sub, 'bo-', linewidth=1.5, markersize=8, 
                    label='Simulation Data', alpha=0.7, markeredgewidth=0.5)
    
    ax2.set_xlabel('Normalized Position (x/L)', fontsize=13)
    ax2.set_ylabel('Temperature T', fontsize=13)
    ax2.set_title('Temperature Profile', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 1)
    ax2.legend(fontsize=11, loc='best')
    
    # Auto-adjust y-limits based on data (use original temp_data for limits)
    if temp_data is not None and T_num is not None:
        y_min = min(np.min(temp_data), np.min(T_num)) * 0.95
        y_max = max(np.max(temp_data), np.max(T_num)) * 1.05
        ax2.set_ylim(y_min, y_max)
    
    plt.tight_layout()
    
    # Save the combined plot
    output_file = 'kep_combined_results.png'
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✅ Combined plot saved: {output_file}")
    
    return fig

def create_individual_combined_plots():
    """Create individual plots for density and temperature with both datasets"""
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Read data (with errors)
    pos_density, density_data, density_error = read_data_file('density_profile.dat')
    pos_temp, temp_data, temp_error = read_data_file('temperature_profile.dat')
    x_num, rho_num, T_num = read_numerical_solution('solution.dat')
    
    # Normalize simulation positions
    if pos_density is not None:
        x_density = normalize_x_axis(pos_density)
    if pos_temp is not None:
        x_temp = normalize_x_axis(pos_temp)
    
    # Subsample simulation data to show only 25 points (with errors)
    x_density_sub, density_data_sub, density_error_sub = subsample_data(x_density, density_data, density_error, n_points=25)
    x_temp_sub, temp_data_sub, temp_error_sub = subsample_data(x_temp, temp_data, temp_error, n_points=25)
    
    # Density Profile - Individual Plot
    fig1, ax1 = plt.subplots(1, 1, figsize=(10, 7))
    
    if x_num is not None and rho_num is not None:
        ax1.plot(x_num, rho_num, 'r-', linewidth=2.5, label='Numerical Algorithm', alpha=0.8)
    
    if x_density_sub is not None and density_data_sub is not None:
        if density_error_sub is not None:
            ax1.errorbar(x_density_sub, density_data_sub, yerr=density_error_sub, 
                       fmt='bo-', linewidth=1.5, markersize=5, capsize=3, capthick=1.5,
                       label='Simulation Data', alpha=0.7, markeredgewidth=0.5, elinewidth=1.5)
        else:
            ax1.plot(x_density_sub, density_data_sub, 'bo-', linewidth=1.5, markersize=5, 
                    label='Simulation Data', alpha=0.7, markeredgewidth=0.5)
        mean_density = np.mean(density_data)
        ax1.axhline(y=mean_density, color='g', linestyle='--', alpha=0.5, linewidth=1,
                   label=f'Mean (Sim): {mean_density:.3f}')
    
    ax1.set_xlabel('Normalized Position (x/L)', fontsize=13)
    ax1.set_ylabel('Density ρ', fontsize=13)
    ax1.set_title(f'KEP Density Profile: Simulation vs Numerical Algorithm (β = 1.0)', 
                 fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    ax1.legend(fontsize=11, loc='best')
    
    plt.tight_layout()
    plt.savefig('density_combined.png', dpi=300, bbox_inches='tight')
    print("✅ Saved: density_combined.png")
    
    # Temperature Profile - Individual Plot
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 7))
    
    if x_num is not None and T_num is not None:
        ax2.plot(x_num, T_num, 'r-', linewidth=2.5, label='Numerical Algorithm', alpha=0.8)
    
    if x_temp_sub is not None and temp_data_sub is not None:
        if temp_error_sub is not None:
            ax2.errorbar(x_temp_sub, temp_data_sub, yerr=temp_error_sub, 
                        fmt='bo-', linewidth=1.5, markersize=5, capsize=3, capthick=1.5,
                        label='Simulation Data', alpha=0.7, markeredgewidth=0.5, elinewidth=1.5)
        else:
            ax2.plot(x_temp_sub, temp_data_sub, 'bo-', linewidth=1.5, markersize=5, 
                    label='Simulation Data', alpha=0.7, markeredgewidth=0.5)
    
    ax2.set_xlabel('Normalized Position (x/L)', fontsize=13)
    ax2.set_ylabel('Temperature T', fontsize=13)
    ax2.set_title(f'KEP Temperature Profile: Simulation vs Numerical Algorithm (β = 1.0)', 
                 fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 1)
    
    if temp_data is not None and T_num is not None:
        y_min = min(np.min(temp_data), np.min(T_num)) * 0.95
        y_max = max(np.max(temp_data), np.max(T_num)) * 1.05
        ax2.set_ylim(y_min, y_max)
    
    ax2.legend(fontsize=11, loc='best')
    plt.tight_layout()
    plt.savefig('temperature_combined.png', dpi=300, bbox_inches='tight')
    print("✅ Saved: temperature_combined.png")

def main():
    """Main function to create all combined plots"""
    
    print("\nStarting combined plotting process...\n")
    
    # Create main combined dashboard
    fig = create_combined_plots()
    
    # Create individual combined plots
    print("\nCreating individual combined plots...")
    create_individual_combined_plots()
    
    # Show plots
    print("\nDisplaying plots...")
    plt.show()
    
    print("\n" + "=" * 60)
    print("Plotting completed successfully!")
    print("=" * 60)
    print("\nGenerated files:")
    print("  - kep_combined_results.png (combined dashboard)")
    print("  - density_combined.png (density comparison)")
    print("  - temperature_combined.png (temperature comparison)")

if __name__ == "__main__":
    main()

