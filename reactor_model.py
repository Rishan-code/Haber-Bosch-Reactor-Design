from kinetics import calculate_rate
from thermo import calculate_heat_of_rxn, calculate_thermal_mass

def pfr_balances(W, y, P0, F_N2_0, inerts_ratio=0.0):
    """
    ODE function for the plug flow reactor.
    y = [X, T]
    Returns [dX/dW, dT/dW]
    """
    X, T = y
    
    # Clamp variables to prevent solver from evaluating unphysical states
    X = min(max(X, 1e-5), 0.999)
    T = max(T, 100.0)
    
    # Get the kinetic rate (-r'_N2)
    minus_r_prime_N2 = calculate_rate(X, T, P0, inerts_ratio)
    
    # Mole balance
    dX_dW = minus_r_prime_N2 / F_N2_0
    
    # Energy balance
    delta_H_rxn = calculate_heat_of_rxn(T)
    thermal_mass = calculate_thermal_mass(X, T, F_N2_0, inerts_ratio)
    
    # dT/dW = (-r'_N2) * (-Delta H_rxn) / (Sum(F_i * Cp_i))
    dT_dW = minus_r_prime_N2 * (-delta_H_rxn) / thermal_mass
    
    return [dX_dW, dT_dW]