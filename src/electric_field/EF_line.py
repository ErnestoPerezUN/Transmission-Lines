# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
# Original script by Ernesto Pérez (port of a MATLAB code); refactored into functions with the help of
# Claude (Anthropic). The solver now lives in charge_simulation.py. See NOTICE.md.
"""Worked example: 500 kV line with 4-bundles and two grounded shield wires.

Run:  ``python src/electric_field/EF_line.py``           (shows the figures)
      ``python src/electric_field/EF_line.py --save DIR`` (writes PNG files instead)

It builds the geometry and applies the **charge simulation** method (``charge_simulation.solve``) with the
shield wires at 0 V. It compares three shield wire strategies (grounded, insulated, none), checks the
surface gradient given by the simulation against the critical gradient given by **Peek's formula**
gradient and draws the layout, the surface gradient and the field map.
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np                                                    # noqa: E402
import corona                                                         # noqa: E402
from charge_simulation import (Wire, add_bundle, charge_per_group, plot_field_map,   # noqa: E402
                               plot_layout, plot_surface_gradient, solve,
                               three_phase_potentials, V_PER_M_TO_KV_PER_CM)

# ---------------------------------------------------------------------------
# Input data: change these numbers to study another line
# ---------------------------------------------------------------------------
V_LL_KV = 500.0            # line-to-line rms voltage [kV]
PHASE_X = [-11.0, 0.0, 11.0]   # horizontal position of phases a, b, c [m]
PHASE_Y = 25.0             # height of the bundle centres [m] (the sag is ignored)
N_SUB, SPACING, R_SUB = 4, 0.45, 0.0141   # bundle: sub-conductors, spacing [m], radius [m]
GUARD_X = [-7.0, 7.0]      # shield wire positions [m]
GUARD_Y, R_GUARD = 35.0, 0.0055
ALTITUDE_M, TEMP_C, M_SURFACE = 1500.0, 20.0, 0.8


def build_line(with_guard=True):
    """Return (wires, group_potentials). Groups: 0, 1, 2 are the phases, 3 and 4 the shield wires."""
    wires = []
    for i, x in enumerate(PHASE_X):
        add_bundle(wires, x, PHASE_Y, N_SUB, SPACING, R_SUB, group=i, label="abc"[i])
    va, vb, vc = three_phase_potentials(V_LL_KV)
    pot = {0: va, 1: vb, 2: vc}
    if with_guard:
        for j, x in enumerate(GUARD_X):
            wires.append(Wire(x, GUARD_Y, R_GUARD, group=3 + j, label=f"G{j + 1}"))
    return wires, pot


def main(save_dir=None):
    delta = corona.air_density(ALTITUDE_M, TEMP_C)
    e_c = corona.peek_critical_gradient(R_SUB * 100, M_SURFACE, delta)
    print(f"delta = {delta:.3f}, E_c from Peek formula = {e_c:.2f} kV/cm (peak) for m = {M_SURFACE}")

    results = {}
    for name in ("grounded", "insulated", "none"):
        wires, pot = build_line(with_guard=(name != "none"))
        if name == "grounded":
            pot.update({3: 0.0, 4: 0.0})       # <- shield wires kept in the system at 0 V
        elif name == "insulated":
            pot.update({3: None, 4: None})     # <- floating: potential unknown, net charge 0
        sol = solve(wires, pot, n_charges=12)
        emax = sol.max_gradient_by_group()
        results[name] = (wires, sol)
        line = ", ".join(f"{'abc'[g]}: {emax[g] * V_PER_M_TO_KV_PER_CM:.2f}" for g in range(3))
        print(f"{name:9s} E_max from charge simulation [kV/cm peak] {line}   simulation error {sol.bc_error:.1e}")
        if name == "insulated":
            print("          insulated wire potential [kV peak]:",
                  round(abs(sol.group_potentials[3]) / 1e3, 1))

    wires, sol = results["grounded"]
    emax = sol.max_gradient_by_group()
    worst = max(emax[g] for g in range(3)) * V_PER_M_TO_KV_PER_CM
    ratio, ok = corona.gradient_margin(worst, e_c)
    print(f"Worst phase: E_max (simulation) / E_c (Peek) = {ratio:.2f} -> {'meets' if ok else 'does NOT meet'} the 0.95 criterion")
    print("Induced charge on the shield wires [uC/m, peak]:",
          [round(abs(charge_per_group(sol)[g]) * 1e6, 3) for g in (3, 4)])

    import matplotlib
    if save_dir:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    guards = (3, 4)
    figs = {
        "layout": plot_layout(wires, guard_groups=guards, zoom_group=1,
                              title="500 kV line: phases (colour) and grounded shield wires (grey)"),
        "surface_gradient": plot_surface_gradient(sol, guards, e_crit_kv_cm=e_c),
        "field_map": plot_field_map(sol, (-25, 25), (0, 45), guards, vmax_kv_m=40),
    }
    if save_dir:
        out = pathlib.Path(save_dir)
        out.mkdir(parents=True, exist_ok=True)
        for k, f in figs.items():
            f.savefig(out / f"ef_line_{k}.png", dpi=150)
        print("Figures written to", out)
    else:
        plt.show()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--save", metavar="DIR", help="write PNG files to DIR instead of showing figures")
    main(ap.parse_args().save)
