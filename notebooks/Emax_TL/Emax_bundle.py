# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
# Code developed in ChatGPT using equation (2.9) from "High Voltage Engineering Fundamentals" by E. Kuffel, W.S. Zaengl, J. Kuffel
# to calculate the voltage gradient Ei [kV/cm] for a given set of parameters.


import math

def voltage_gradient(Ci, n2, r, s, U):
    """
    Calculate the voltage gradient Ei [kV/cm] .
    
    Parameters:
        Ci : float
            Capacitance per unit length [F/m]
        n2 : int
            Number of conductors per bundle
        r  : float
            Radius of sub-conductor [m]
        s  : float
            Distance between sub-conductors [m]
        U  : float
            Rated voltage [kV]
            
    Returns:
        Ei : float
            Voltage gradient [kV/cm]
    """
    
    # Constants
    epsilon_0 = 8.854e-12  # F/m
    
    # Compute the term inside brackets
    bracket_term = 1 + 2 * (r / s) * (n2 - 1) * math.sin(math.pi / n2)
    
    # Equation (2.9)
    Ei = (Ci / (2 * math.pi * epsilon_0 * n2 * r)) * bracket_term * (U / (math.sqrt(3) * 100))
    
    return Ei


# Example usage
if __name__ == "__main__":
    # Example parameters
    Ci = 1e-11      # F/m (example value)
    n2 = 1          # conductors per bundle
    r = 0.01       # m (15 mm)
    s = .457        # m (400 mm)
    U = 115         # kV

    Ei = voltage_gradient(Ci, n2, r, s, U)
    print(f"Voltage gradient Ei = {Ei:.4f} kV/cm")
