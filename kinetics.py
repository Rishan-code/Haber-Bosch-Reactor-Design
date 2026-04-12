import math

def calculate_rate(X, T, P0):
    """
    Calculate the rate of reaction using the Temkin-Pyzhev kinetics.
    Returns -r'_N2.
    """
    R = 8.314  # J/(mol K)
    
    # Arrhenius Parameters
    E1 = 170660.0  # J/mol
    E2 = 150420.0  # J/mol
    A1 = 1.79e14
    A2 = 2.5e10
    alpha = 0.5
    
    k1 = A1 * math.exp(-E1 / (R * T))
    k2 = A2 * math.exp(-E2 / (R * T))
    
    # Molar flow rates proportional calculations
    n_N2 = 1.0 - X
    n_H2 = 3.0 * (1.0 - X)
    n_NH3 = 2.0 * X
    n_total = n_N2 + n_H2 + n_NH3
    
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