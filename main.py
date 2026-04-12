import matplotlib
matplotlib.use('Agg')  # Critical GUI Fix
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from reactor_model import pfr_balances

def main():
    # Operating Conditions
    T0 = 673.15      # K
    P0 = 200.0       # atm
    F_N2_0 = 100.0   # mol/s
    X0 = 0.001       # Critical Math Fix: Start > 0 to prevent div by zero
    
    # Reactor Sizing & Cooling
    num_beds = 3
    W_bed = 2500.0   # kg of catalyst per bed
    delta_T_cooling = 60.0  # K
    
    W_all = []
    X_all = []
    T_all = []
    
    current_X = X0
    current_T = T0
    
    for i in range(num_beds):
        W_start = i * W_bed
        W_end = (i + 1) * W_bed
        
        if i > 0:
            # Inter-stage cooling
            current_T -= delta_T_cooling
            
            # Store the point representing the vertical drop in temperature
            # at the same catalyst weight between beds
            W_all.append(W_start)
            X_all.append(current_X)
            T_all.append(current_T)
        
        print(f"\n[{'='*30}]")
        print(f"BED {i+1} CONDITIONS:")
        print(f"  -> Entrance: Conversion X = {current_X:.4f}, Temperature T = {current_T:.2f} K")
            
        y0 = [current_X, current_T]
        
        # Solver Setup: Must use Radau for stiff chemical ODEs
        sol = solve_ivp(
            fun=lambda W, y: pfr_balances(W, y, P0, F_N2_0),
            t_span=[W_start, W_end],
            y0=y0,
            method='Radau',
            dense_output=True
        )
        
        # Evaluate using dense output for smooth plotting
        W_eval = np.linspace(W_start, W_end, 200)
        y_eval = sol.sol(W_eval)
        
        if i == 0:
            W_all.extend(W_eval)
            X_all.extend(y_eval[0])
            T_all.extend(y_eval[1])
        else:
            # Skip the first point to smoothly transition over intervals 
            W_all.extend(W_eval[1:])
            X_all.extend(y_eval[0][1:])
            T_all.extend(y_eval[1][1:])
            
        # Update states for the next bed
        current_X = sol.y[0][-1]
        current_T = sol.y[1][-1]
        
        print(f"  -> Exit:     Conversion X = {current_X:.4f}, Temperature T = {current_T:.2f} K")

    print(f"\n[{'='*30}]")
    print("FINAL REACTOR METRICS:")
    print(f"  -> Total Catalyst Weight: {num_beds * W_bed} kg")
    print(f"  -> Final Ammonia Production Rate: {2.0 * F_N2_0 * current_X:.2f} mol/s\n")

    # Convert lists to Numpy arrays
    W_all = np.array(W_all)
    X_all = np.array(X_all)
    T_all = np.array(T_all)
    
    # --- Plot 1: Reactor Profiles ---
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    
    color_x = 'tab:blue'
    ax1.set_xlabel('Catalyst Weight, W (kg)')
    ax1.set_ylabel('Conversion (X)', color=color_x)
    ax1.plot(W_all, X_all, color=color_x, linewidth=2, label='Conversion')
    ax1.tick_params(axis='y', labelcolor=color_x)
    
    ax2 = ax1.twinx()
    color_t = 'tab:red'
    ax2.set_ylabel('Temperature (K)', color=color_t)
    ax2.plot(W_all, T_all, color=color_t, linewidth=2, label='Temperature')
    ax2.tick_params(axis='y', labelcolor=color_t)
    
    plt.title('Multi-Bed PFR Profiles: Conversion and Temperature vs Catalyst Weight')
    fig1.tight_layout()
    fig1.savefig('reactor_profiles.png')
    plt.close(fig1)
    
    # --- Plot 2: Sawtooth Trajectory ---
    fig2, ax3 = plt.subplots(figsize=(8, 6))
    
    ax3.plot(T_all, X_all, 'k-', linewidth=2)
    ax3.set_xlabel('Temperature (K)')
    ax3.set_ylabel('Conversion (X)')
    ax3.set_title('Sawtooth Operating Trajectory: Conversion vs. Temperature')
    ax3.grid(True)
    
    fig2.tight_layout()
    fig2.savefig('sawtooth_trajectory.png')
    plt.close(fig2)

if __name__ == "__main__":
    main()