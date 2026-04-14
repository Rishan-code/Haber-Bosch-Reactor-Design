import matplotlib
matplotlib.use('Agg')  # Critical GUI Fix
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from reactor_model import pfr_balances
from plotting import generate_analytical_plots

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
        # 4% inerts (Argon, Methane) in feed, given basic 3:1 stoich.
        # N2: 24%, H2: 72%, Inerts: 4% -> inerts_ratio = 4.0 / 24.0
        sol = solve_ivp(
            fun=lambda W, y: pfr_balances(W, y, P0, F_N2_0, inerts_ratio=4.0/24.0),
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
    # Generate analytical visuals from our new module
    print("Generating analysis plots...")
    generate_analytical_plots(W_all, X_all, T_all, P0, inerts_ratio=4.0/24.0)
    print("Done! Plots saved safely.")

if __name__ == "__main__":
    main()