import math

def calculate_Kp(T):
    """
    Empirical equilibrium constant for N2 + 3H2 <-> 2NH3.
    """
    # Gillespie-Beattie correlation for Kp (atm^-2)
    log10_Keq = -2.691122 * math.log10(T) - 5.519265e-5 * T + 1.848863e-7 * T**2 + 2001.6 / T + 2.6899
    # Clamp to prevent overflow warning during solver exploration
    log10_Keq = min(log10_Keq, 300.0)
    return 10**log10_Keq

def calculate_rate(X, T, P0, inerts_ratio=0.0):
    """
    Calculate the rate of reaction using the Temkin-Pyzhev kinetics.
    Returns -r'_N2.
    """
    R = 8.314  # J/(mol K)
    
    # Arrhenius Parameters
    E1 = 170660.0  # J/mol
    A1 = 1.79e14   # kmol/m³·h·atm^0.5 roughly, typical for models
    alpha = 0.5
    
    k1 = A1 * math.exp(-E1 / (R * T))
    
    # Use equilibrium constant Kp instead of Arrhenius for backward rate
    Kp = calculate_Kp(T)
    k2 = k1 / Kp
    
    # Molar flow rates proportional calculations
    n_N2 = 1.0 - X
    n_H2 = 3.0 * (1.0 - X)
    n_NH3 = 2.0 * X
    n_inerts = inerts_ratio
    n_total = n_N2 + n_H2 + n_NH3 + n_inerts
    
    # Mole fractions
    y_N2 = n_N2 / n_total
    y_H2 = n_H2 / n_total
    y_NH3 = n_NH3 / n_total
    
    # Partial pressures
    P_N2 = y_N2 * P0
    P_H2 = y_H2 * P0
    P_NH3 = y_NH3 * P0
    
    # Critical Fix: Clamp P_NH3 to a minimum of 1e-5 to prevent division by zero
    P_NH3 = max(P_NH3, 1e-5)
    
    # Temkin-Pyzhev Rate Equation
    term1 = k1 * P_N2 * ((P_H2**3) / (P_NH3**2))**alpha
    term2 = k2 * ((P_NH3**2) / (P_H2**3))**(1.0 - alpha)
    
    return term1 - term2