import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from kinetics import calculate_rate, calculate_Kp
from scipy.optimize import root_scalar

def calculate_x_eq(T, P0, inerts_ratio):
    """
    Solves for the equilibrium conversion X_eq at a given Temperature and Pressure.
    """
    Kp = calculate_Kp(T)
    
    def eq_residual(X):
        # Prevent math domain errors from invalid X
        if X >= 0.999:
            return 1e6
        if X <= 1e-5:
            return -1e6
            
        n_total = 4.0 - 2.0 * X + inerts_ratio
        P_N2 = ((1.0 - X) / n_total) * P0
        P_H2 = ((3.0 - 3.0 * X) / n_total) * P0
        P_NH3 = ((2.0 * X) / n_total) * P0
        
        return (P_NH3**2 / (P_N2 * P_H2**3)) - Kp
        
    try:
        sol = root_scalar(eq_residual, bracket=[1e-5, 0.99], method='brentq')
        if sol.converged:
            return sol.root
        return 0.0
    except ValueError:
        return 0.0

def generate_analytical_plots(W_all, X_all, T_all, P0, inerts_ratio):
    """
    Generates crucial analytical plots for reactor analysis.
    """
    # 1. Reactor Profiles (Conversion & Temp)
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    color_x = 'tab:blue'
    ax1.set_xlabel('Catalyst Weight, W (kg)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Conversion (X)', color=color_x, fontsize=12, fontweight='bold')
    ax1.plot(W_all, X_all, color=color_x, linewidth=2.5, label='Conversion')
    ax1.tick_params(axis='y', labelcolor=color_x)
    ax1.grid(True, alpha=0.3)
    
    ax2 = ax1.twinx()
    color_t = 'tab:red'
    ax2.set_ylabel('Temperature (K)', color=color_t, fontsize=12, fontweight='bold')
    ax2.plot(W_all, T_all, color=color_t, linewidth=2.5, linestyle='--', label='Temperature')
    ax2.tick_params(axis='y', labelcolor=color_t)
    
    plt.title('Multi-Bed PFR Profiles: Conversion and Temperature', fontsize=14, fontweight='bold')
    fig1.tight_layout()
    fig1.savefig('reactor_profiles.png', dpi=300)
    plt.close(fig1)

    # 2. Sawtooth Trajectory with Equilibrium Curve
    # Generate range of temperatures to plot the equilibrium curve
    T_min = min(T_all) - 20
    T_max = max(T_all) + 50
    T_range = np.linspace(T_min, T_max, 100)
    X_eq = [calculate_x_eq(T, P0, inerts_ratio) for T in T_range]
    
    fig2, ax3 = plt.subplots(figsize=(9, 6))
    ax3.plot(T_all, X_all, 'k-', linewidth=2.5, label='Actual Operating Trajectory')
    ax3.plot(T_range, X_eq, 'r--', linewidth=2, label='Equilibrium Limit (X_eq)')
    
    ax3.set_xlabel('Temperature (K)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Conversion (X)', fontsize=12, fontweight='bold')
    ax3.set_title('Sawtooth Operating Trajectory vs. Equilibrium Curve', fontsize=14, fontweight='bold')
    ax3.set_xlim(T_min, T_max)
    ax3.set_ylim(0, max(max(X_all)+0.1, 0.7))
    ax3.grid(True, alpha=0.5)
    ax3.legend(loc='lower right', fontsize=11)
    
    fig2.tight_layout()
    fig2.savefig('sawtooth_equilibrium.png', dpi=300)
    plt.close(fig2)
    
    # 3. Mole Fractions vs Catalyst Weight
    y_N2, y_H2, y_NH3, y_inerts = [], [], [], []
    for X in X_all:
        n_total = 4.0 - 2.0 * X + inerts_ratio
        y_N2.append((1.0 - X) / n_total)
        y_H2.append((3.0 - 3.0 * X) / n_total)
        y_NH3.append((2.0 * X) / n_total)
        y_inerts.append(inerts_ratio / n_total)
        
    fig3, ax4 = plt.subplots(figsize=(10, 6))
    ax4.plot(W_all, y_H2, label='Hydrogen (H2)', linewidth=2.5, color='tab:orange')
    ax4.plot(W_all, y_N2, label='Nitrogen (N2)', linewidth=2.5, color='tab:blue')
    ax4.plot(W_all, y_NH3, label='Ammonia (NH3)', linewidth=2.5, color='tab:green')
    ax4.plot(W_all, y_inerts, label='Inerts (Ar, CH4)', linewidth=2, color='gray', linestyle=':')
    
    ax4.set_xlabel('Catalyst Weight, W (kg)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Mole Fraction (y_i)', fontsize=12, fontweight='bold')
    ax4.set_title('Gas Composition along the Reactor', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=11)
    fig3.tight_layout()
    fig3.savefig('composition_profile.png', dpi=300)
    plt.close(fig3)
    
    # 4. Reaction Rate Profile
    rates = [calculate_rate(X, T, P0, inerts_ratio) for X, T in zip(X_all, T_all)]
    fig4, ax5 = plt.subplots(figsize=(10, 6))
    ax5.plot(W_all, rates, color='purple', linewidth=2.5)
    ax5.set_xlabel('Catalyst Weight, W (kg)', fontsize=12, fontweight='bold')
    ax5.set_ylabel("Reaction Rate ($-r'_{N2}$, mol/kg·s)", fontsize=12, fontweight='bold')
    ax5.set_title('Reaction Rate vs. Catalyst Weight', fontsize=14, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    fig4.tight_layout()
    fig4.savefig('rate_profile.png', dpi=300)
    plt.close(fig4)
