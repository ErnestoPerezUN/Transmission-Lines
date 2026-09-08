# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
import numpy as np
import matplotlib.pyplot as plt

# ==============================
# Simulación de Cargas
# Ernesto Pérez (refactored to Python)
# ==============================

# ------------------------------
# Input Data
# ------------------------------

num_cargas = 6          # Number of simulation charges per conductor
scl = 30                # Plotting scale
det = 10                # Grid detail (points per diameter)

# Define conductors (list of dictionaries)
conductors = [
    {
        "Pos": np.array([0.0, 15.0]),      # [x,y] in meters
        "Pot": 220/np.sqrt(3),             # Potential [kV]
        "diametro": 0.010                  # Diameter [m]
    },
        {
        "Pos": np.array([0.0, 20.0]),      # [x,y] in meters
        "Pot": 220/np.sqrt(3),             # Potential [kV]
        "diametro": 0.010                  # Diameter [m]
    }
]

# ------------------------------
# Initialization
# ------------------------------

num_con = len(conductors)
angc = 2 * np.pi / num_cargas
k = 1 / (2 * np.pi * 8.85e-12)

carx = []
cary = []
Pcx = []
Pcy = []
Perrx = []
Perry = []
Pot = []

for i, cond in enumerate(conductors):

    if num_cargas == 1:
        r_carga = 0
    else:
        r_carga = cond["diametro"] / 4

    for j in range(1, num_cargas + 1):
        angle = j * angc

        # Charge positions
        carx.append(cond["Pos"][0] + r_carga * np.cos(angle))
        cary.append(cond["Pos"][1] + r_carga * np.sin(angle))

        # Boundary potential points
        Pcx.append(cond["Pos"][0] + cond["diametro"] * np.cos(angle) / 2)
        Pcy.append(cond["Pos"][1] + cond["diametro"] * np.sin(angle) / 2)

        # Error evaluation points
        Perrx.append(cond["Pos"][0] + cond["diametro"] * np.cos(1.5 * angle) / 2)
        Perry.append(cond["Pos"][1] + cond["diametro"] * np.sin(1.5 * angle) / 2)

        Pot.append(cond["Pot"])

carx = np.array(carx)
cary = np.array(cary)
Pcx = np.array(Pcx)
Pcy = np.array(Pcy)
Perrx = np.array(Perrx)
Perry = np.array(Perry)
Pot = np.array(Pot) 

# ------------------------------
# Simulation of Charges Method
# ------------------------------

n = len(Pcy)
M1 = np.zeros((n, n))
M1err = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        r1 = np.sqrt((Pcy[i] - cary[j])**2 + (Pcx[i] - carx[j])**2)
        r2 = np.sqrt((Pcy[i] + cary[j])**2 + (Pcx[i] - carx[j])**2)

        r1err = np.sqrt((Perry[i] - cary[j])**2 + (Perrx[i] - carx[j])**2)
        r2err = np.sqrt((Perry[i] + cary[j])**2 + (Perrx[i] - carx[j])**2)

        M1[i, j] = k * np.log(r2 / r1)
        M1err[i, j] = k * np.log(r2err / r1err)

# Solve for charge densities
ro = np.linalg.solve(M1, Pot)

Poterr = M1err @ ro
error = (Pot - Poterr) / Pot

print("Relative Error (max): ", np.max(np.abs(error)))

# ------------------------------
# Field and Potential Calculation
# ------------------------------

for kk, cond in enumerate(conductors):

    Ax = np.arange(
        cond["Pos"][0] - cond["diametro"] * scl,
        cond["Pos"][0] + cond["diametro"] * scl,
        cond["diametro"] / det
    )

    Ay = np.arange(
        cond["Pos"][1] - cond["diametro"] * scl,
        cond["Pos"][1] + cond["diametro"] * scl,
        cond["diametro"] / det
    )

    Po = np.zeros((len(Ay), len(Ax)))

    for l, y in enumerate(Ay):
        for i, x in enumerate(Ax):

            inside = False
            M2 = np.zeros(len(carx))

            for j in range(len(carx)):

                r1 = np.sqrt((carx[j] - x)**2 + (cary[j] - y)**2)

                if r1 <= cond["diametro"] / 2:
                    inside = True
                    P = Pot[j]
                else:
                    r2 = np.sqrt((carx[j] - x)**2 + (cary[j] + y)**2)
                    M2[j] = k * np.log(r2 / r1)

            if inside:
                Po[l, i] = P
            else:
                Po[l, i] = M2 @ ro

    # ------------------------------
    # Electric Field
    # ------------------------------

    dx = cond["diametro"] / det
    dy = cond["diametro"] / det

    Ey, Ex = np.gradient(Po, dy, dx)
    Ex = -Ex
    Ey = -Ey

    E_mag = np.sqrt(Ex**2 + Ey**2)

    print(f"Conductor {kk+1}")
    print("Emax [kV/m] =", np.max(E_mag))

    X, Y = np.meshgrid(Ax, Ay)

    # ------------------------------
    # Plots
    # ------------------------------

    plt.figure()
    plt.contourf(X, Y, Po, levels=15)
    plt.colorbar(label="Potential [kV]")
    plt.title(f"Potential Distribution - Conductor {kk+1}")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.show()

    plt.figure()
    plt.quiver(X, Y, Ex, Ey)
    plt.title(f"Electric Field Vectors - Conductor {kk+1}")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.show()

    plt.figure()
    plt.contourf(X, Y, E_mag, levels=15)
    plt.colorbar(label="|E| [kV/m]")
    plt.title(f"Electric Field Magnitude - Conductor {kk+1}")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.show()