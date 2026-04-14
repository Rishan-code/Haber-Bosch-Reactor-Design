def calculate_Cp(T, component):
    """
    Returns the heat capacity in J/(mol K) for a given component at temperature T (K).
    """
    if component == 'N2':
        A, B, C, D = 29.00, 0.2199e-2, 0.5723e-5, -2.871e-9
    elif component == 'H2':
        A, B, C, D = 28.84, 0.00765e-2, 0.3288e-5, -0.8698e-9
    elif component == 'NH3':
        A, B, C, D = 27.31, 2.383e-2, 1.707e-5, -1.185e-8
    elif component == 'Inerts':
        A, B, C, D = 29.00, 0.2199e-2, 0.5723e-5, -2.871e-9  # Roughly approximate to N2/Argon
    else:
        A, B, C, D = 29.0, 0.0, 0.0, 0.0
        
    return A + B*T + C*T**2 + D*T**3

def calculate_heat_of_rxn(T):
    """
    Returns the heat of reaction in J/mol for the synthesis of ammonia
    using Kirchhoff's equation.
    """
    A_N2, B_N2, C_N2, D_N2 = 29.00, 0.2199e-2, 0.5723e-5, -2.871e-9
    A_H2, B_H2, C_H2, D_H2 = 28.84, 0.00765e-2, 0.3288e-5, -0.8698e-9
    A_NH3, B_NH3, C_NH3, D_NH3 = 27.31, 2.383e-2, 1.707e-5, -1.185e-8
    
    Delta_A = 2 * A_NH3 - (A_N2 + 3 * A_H2)
    Delta_B = 2 * B_NH3 - (B_N2 + 3 * B_H2)
    Delta_C = 2 * C_NH3 - (C_N2 + 3 * C_H2)
    Delta_D = 2 * D_NH3 - (D_N2 + 3 * D_H2)
    
    T_ref = 298.15
    Delta_H_298 = -92220.0  # J/mol N2
    
    integral_Cp = (
        Delta_A * (T - T_ref) + 
        (Delta_B / 2.0) * (T**2 - T_ref**2) + 
        (Delta_C / 3.0) * (T**3 - T_ref**3) + 
        (Delta_D / 4.0) * (T**4 - T_ref**4)
    )
    
    return Delta_H_298 + integral_Cp

def calculate_thermal_mass(X, T, F_N2_0, inerts_ratio=0.0):
    """
    Calculates sum(F_i * C_p,i) given the conversion X, temperature T,
    initial molar flow rate of N2, and the molar ratio of inerts to N2 feed.
    """
    Cp_N2 = calculate_Cp(T, 'N2')
    Cp_H2 = calculate_Cp(T, 'H2')
    Cp_NH3 = calculate_Cp(T, 'NH3')
    Cp_Inerts = calculate_Cp(T, 'Inerts')
    
    # Molar flow rates
    F_N2 = F_N2_0 * (1.0 - X)
    F_H2 = 3.0 * F_N2_0 * (1.0 - X)
    F_NH3 = 2.0 * F_N2_0 * X
    F_Inerts = F_N2_0 * inerts_ratio
    
    # Thermal mass = sum(F_i * Cp_i)
    return F_N2 * Cp_N2 + F_H2 * Cp_H2 + F_NH3 * Cp_NH3 + F_Inerts * Cp_Inerts