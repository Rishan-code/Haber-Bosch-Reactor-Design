def calculate_heat_of_rxn(T):
    """
    Returns the heat of reaction in J/mol for the synthesis of ammonia.
    Approximated as constant across temperature.
    """
    return -92220.0

def calculate_thermal_mass(X, T, F_N2_0):
    """
    Calculates sum(F_i * C_p,i) given the conversion X, temperature T,
    and initial molar flow rate of N2.
    """
    # Constant average heat capacities in J/(mol K)
    Cp_N2 = 29.1
    Cp_H2 = 28.8
    Cp_NH3 = 35.1
    
    # Molar flow rates
    F_N2 = F_N2_0 * (1.0 - X)
    F_H2 = 3.0 * F_N2_0 * (1.0 - X)
    F_NH3 = 2.0 * F_N2_0 * X
    
    # Thermal mass = sum(F_i * Cp_i)
    return F_N2 * Cp_N2 + F_H2 * Cp_H2 + F_NH3 * Cp_NH3