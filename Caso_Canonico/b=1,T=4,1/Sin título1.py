import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
import scipy.special as sp
import os

# ==========================================
# 1. PARÁMETROS DEL SISTEMA
# ==========================================
beta = 1.0
T_L = 4.0
T_R = 1.0
rho_bar = 0.8  # Densidad media (N=1600 / L=2000)
L_sim = 2000   # Tamaño de la red KMC

# ==========================================
# 2. SOLVER ANALÍTICO EXACTO (MFT)
# ==========================================
print("Resolviendo las ecuaciones hidrodinámicas de MFT asumiendo LEA...")

# Constantes Gamma
g1 = sp.gamma(beta + 1.0)
g2 = sp.gamma(beta + 2.0)

def coeficientes_D(rho, T):
    """Matriz de transporte exacta derivada por Hurtado et al. (2019)"""
    D11 = g1 * T**beta
    D12 = beta * g1 * rho * (1.0 - rho) * T**(beta - 1.0)
    D21 = g2 * T**(beta + 1.0)
    # Expresión exacta para D22
    D22 = g2 * T**beta * rho * ( (1.0 + beta) + rho * ((beta**2 - 7.0*beta - 6.0)/12.0) )
    return D11, D12, D21, D22

def sistema_EDO(x, y, p):
    """
    y[0] = T(x)
    y[1] = I_rho(x) (Integral acumulada de densidad para forzar rho_bar)
    p[0] = A (Constante de estratificación de Soret)
    p[1] = c2 (Corriente de energía estacionaria)
    """
    T = y[0]
    A = p[0]
    c2 = p[1]
    
    # Ecuación de estado asumiendo LEA (c1 = 0)
    rho = A / (T**beta + A)
    
    D11, D12, D21, D22 = coeficientes_D(rho, T)
    
    # Conductividad térmica efectiva
    kappa_eff = D22 - D21 * (D12 / D11)
    
    dT_dx = c2 / kappa_eff
    dIrho_dx = rho  
    
    return np.vstack((dT_dx, dIrho_dx))

def condiciones_contorno(ya, yb, p):
    return np.array([
        ya[0] - T_L,        # Frontera caliente
        yb[0] - T_R,        # Frontera fría
        ya[1] - 0.0,        # Integral acumulada empieza en 0
        yb[1] - rho_bar     # Conservación de masa total
    ])

# Malla y estimación inicial (Guess)
x_mesh = np.linspace(0, 1, 500)
y_guess = np.vstack((np.linspace(T_L, T_R, 500), np.linspace(0, rho_bar, 500)))
p_guess = np.array([10.0, -5.0]) 

sol = solve_bvp(sistema_EDO, condiciones_contorno, x_mesh, y_guess, p=p_guess)

if sol.success:
    print(f"✓ Convergencia MFT exitosa (A={sol.p[0]:.4f}, Corriente calor={sol.p[1]:.4f})")
    x_th = np.linspace(0, 1, 1000)
    T_th = sol.sol(x_th)[0]
    rho_th = sol.p[0] / (T_th**beta + sol.p[0])
else:
    print("X Fallo en el solver MFT. Revisa los parámetros.")
    exit()

# ==========================================
# 3. LECTURA DE DATOS NUMÉRICOS KMC
# ==========================================
print("Cargando datos de la simulación KMC...")

def cargar_dat(filename):
    if os.path.exists(filename):
        data = np.loadtxt(filename, comments='#')
        x_num = data[:, 0] / L_sim
        y_num = data[:, 1]
        err_num = data[:, 2] if data.shape[1] > 2 else None
        return x_num, y_num, err_num
    return None, None, None

x_rho, rho_num, err_rho = cargar_dat('density_profile.dat')
x_T, T_num, err_T = cargar_dat('temperature_profile.dat')
x_C, C_num, err_C = cargar_dat('correlation_profile.dat')

# ==========================================
# 4. REPRESENTACIÓN GRÁFICA TIPO TFG
# ==========================================
plt.rcParams.update({'font.size': 12, 'axes.labelsize': 14, 'lines.markersize': 2})
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

skip = 15 # Dibuja 1 de cada 15 puntos para no saturar

# PANEL 1: Densidad
axs[0].plot(x_th, rho_th, 'k-', linewidth=2, label='MFT (Límite $\infty$)')
if rho_num is not None:
    axs[0].errorbar(x_rho[::skip], rho_num[::skip], yerr=(err_rho[::skip] if err_rho is not None else None), 
                    fmt='o', color='tab:blue', alpha=0.7, label='KMC (L=2000)')
axs[0].set_title(r'Perfil de Densidad $\rho(x)$')
axs[0].set_xlabel('Posición $x/L$')
axs[0].set_ylabel(r'Densidad $\rho$')
axs[0].legend()
axs[0].grid(True, alpha=0.3)

# PANEL 2: Temperatura
axs[1].plot(x_th, T_th, 'k-', linewidth=2, label='MFT (Límite $\infty$)')
if T_num is not None:
    axs[1].errorbar(x_T[::skip], T_num[::skip], yerr=(err_T[::skip] if err_T is not None else None), 
                    fmt='o', color='tab:red', alpha=0.7, label='KMC (L=2000)')
axs[1].set_title(r'Perfil de Temperatura $T(x)$')
axs[1].set_xlabel('Posición $x/L$')
axs[1].set_ylabel(r'Temperatura $T$')
axs[1].legend()
axs[1].grid(True, alpha=0.3)

# PANEL 3: Correlación
if C_num is not None:
    axs[2].plot(x_C[::skip], C_num[::skip], 'o', color='tab:purple', alpha=0.7, label='KMC')
axs[2].axhline(0, color='k', linestyle='--', linewidth=2, label='Predicción LEA (=0)')
axs[2].set_title(r'Correlación $\langle n_i E_i \rangle - \langle n_i \rangle\langle E_i \rangle$')
axs[2].set_xlabel('Posición $x/L$')
axs[2].set_ylabel('Covarianza Local')
axs[2].legend()
axs[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('analisis_ruptura_LEA.png', dpi=300)
print("-> Gráfico final guardado como 'analisis_ruptura_LEA.png'")
plt.show()