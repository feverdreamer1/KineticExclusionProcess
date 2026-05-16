
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

base_dir = os.path.dirname(os.path.abspath(__file__))

cases = [
    ("N=100 T=2,1 b=1,rho=0.5 0.7", 2.0, 1.0, 0.5, 0.7, "T=2.0, ρ=0.5→0.7"),
    ("N=100 T=4,1 b=1 rho=0.5 0.8", 4.0, 1.0, 0.5, 0.8, "T=4.0, ρ=0.5→0.8"),
    ("N=100 T=6,1 b=1 rho=0.6 0.9", 6.0, 1.0, 0.6, 0.9, "T=6.0, ρ=0.6→0.9"),
]

def read_theoretical_data(folder_path):
    file_path = os.path.join(folder_path, "solution.csv")
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found")
        return None, None, None

    try:
        data = np.loadtxt(file_path, delimiter=',', skiprows=1)
        x_theo = data[:, 0]
        rho_theo = data[:, 1]
        T_theo = data[:, 2]
        return x_theo, rho_theo, T_theo
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None, None, None


def read_density_file(folder_path):
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
            errors = data[:, 2] if data.shape[1] > 2 else None

        L = len(positions)
        x_normalized = (positions - 1) / (L - 1) if L > 1 else positions

        return x_normalized, density, errors

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None, None, None


def read_temperature_file(folder_path):
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

        L = len(positions)
        x_normalized = (positions - 1) / (L - 1) if L > 1 else positions

        return x_normalized, temperature

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None, None


def plot_density_comparison():
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    colors = ['blue', 'green', 'red']
    markers = ['o', 's', '^']

    for idx, (folder, T_left, T_right, rho_left, rho_right, label) in enumerate(cases):
        folder_path = os.path.join(base_dir, folder)

        x, density, errors = read_density_file(folder_path)
        x_theo, rho_theo, _ = read_theoretical_data(folder_path)

        if x is not None and density is not None:

            if errors is not None:
                ax.errorbar(
                    x[::40], density[::40], yerr=errors[::40],
                    fmt=markers[idx], color=colors[idx],
                    markersize=12, capsize=2,
                    label=f'Sim: {label}', alpha=0.7,
                    linestyle='None'
                )
            else:
                ax.plot(
                    x[::40], density[::40],
                    marker=markers[idx], color=colors[idx],
                    markersize=12, linestyle='None',
                    label=f'Sim: {label}', alpha=0.7
                )

            if x_theo is not None and rho_theo is not None:
                ax.plot(
                    x_theo, rho_theo,
                    '-', color=colors[idx],
                    linewidth=2,
                    label=f'Theoretical: {label}',
                    alpha=0.8
                )

    ax.set_xlabel('Normalized position (x/L)', fontsize=12)
    ax.set_ylabel('Density ρ', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.45, 0.9)   # ← AQUÍ
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig('density_comparison_beta1.png', dpi=300, bbox_inches='tight')
    print("Saved: density_comparison_beta1.png")
    return fig


def plot_temperature_comparison():
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    colors = ['blue', 'green', 'red']
    markers = ['o', 's', '^']

    for idx, (folder, T_left, T_right, rho_left, rho_right, label) in enumerate(cases):
        folder_path = os.path.join(base_dir, folder)

        x, temperature = read_temperature_file(folder_path)
        x_theo, _, T_theo = read_theoretical_data(folder_path)

        if x is not None and temperature is not None:

            ax.plot(
                x[::40], temperature[::40],
                marker=markers[idx], color=colors[idx],
                markersize=12, linestyle='None',
                label=f'Sim: {label}', alpha=0.7
            )

            if x_theo is not None and T_theo is not None:
                ax.plot(
                    x_theo, T_theo,
                    '-', color=colors[idx],
                    linewidth=2,
                    label=f'Theoretical: {label}',
                    alpha=0.8
                )

    ax.set_xlabel('Normalized position (x/L)', fontsize=12)
    ax.set_ylabel('Temperature T', fontsize=12)
    ax.set_xlim(0, 1)  # ← AQUÍ
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig('temperature_comparison_beta1.png', dpi=300, bbox_inches='tight')
    print("Saved: temperature_comparison_beta1.png")
    return fig


def main():
    print("=" * 60)
    print("Generando gráficas para β=1")
    print("=" * 60)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    plot_density_comparison()
    plot_temperature_comparison()

    print("\nGráficas generadas correctamente.")
    plt.show()


if __name__ == "__main__":
    main()
