# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
# Written with the help of Claude (Anthropic); see NOTICE.md.
"""A three-phase line described by a few parameters, and helpers to vary them.

Built for workshops: describe a line once with ``LineCase``, then ask for the results of a variation
with ``analyze(case, spacing=0.30)`` or for a whole list of values with ``sweep``. Nothing is changed
in the original case, so a variation is always compared with the base line.

    case = LineCase()                                   # base line (500 kV, 3 sub-conductors, grounded guards)
    analyze(case)                                       # results of the base line
    analyze(case, n_sub=4)                              # the same line with 4 sub-conductors
    sweep(case, "n_sub", [2, 4])                        # one row per value

Which method produced each number
---------------------------------
Every result is labelled with the method that produced it:

* ``(simulación)``: **charge simulation** (``charge_simulation.py``), a numerical method with the exact
  geometry of all the conductors, the guards and the ground. The surface gradients ``E_max`` and the
  ground field come from it.
* ``(Kuffel)``: **formula** of Kuffel / CIGRE for the maximum gradient of a bundle, fed with the phase
  charge of the Maxwell potential matrix (``corona.phase_charges``). It is the approximate way of getting
  ``E_max`` without the simulation, and the analysis compares it with the simulation.
* ``(Peek)``: **formula** of Peek for the critical gradient ``E_c``.

Gradients are **peak** values in kV/cm and the ground field is **rms** in kV/m.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np

import corona
from charge_simulation import (V_PER_M_TO_KV_PER_CM, Wire, add_bundle, solve, three_phase_potentials)

PHASE_NAMES = "abc"


@dataclass(frozen=True)
class LineCase:
    """Parameters of a line. Positions are ``(x, y)`` in metres; y is the height above ground."""
    v_ll_kv: float = 500.0                                                  # line-to-line rms voltage [kV]
    phases: tuple = ((-11.0, 25.0), (0.0, 25.0), (11.0, 25.0))              # bundle centre of phases a, b, c
    n_sub: int = 3                                                          # sub-conductors per phase
    spacing: float = 0.45                                                   # distance between adjacent sub-conductors [m]
    r_sub: float = 0.0141                                                   # radius of ONE sub-conductor [m]
    guards: tuple = ((-7.0, 35.0), (7.0, 35.0))                             # shield wires; () for none
    r_guard: float = 0.0055                                                 # shield wire radius [m]
    guard_mode: str = "grounded"                                            # "grounded" (0 V), "insulated" or "none"
    altitude_m: float = 1500.0                                              # [m]
    temp_c: float = 20.0                                                    # [degC]
    m_surface: float = 0.8                                                  # Peek surface factor
    sequence: str = "abc"                                                   # phase sequence, "abc" or "acb"
    n_charges: int = 16                                                     # fictitious line charges per conductor (charge simulation)

    def replace(self, **changes) -> "LineCase":
        """Copy of the case with some parameters changed. An unknown name raises ``TypeError``."""
        return dataclasses.replace(self, **changes)


def build(case: LineCase):
    """Wires and group potentials [V, complex peak] of the case. Groups 0-2 are the phases, 3+ the guards."""
    wires = []
    for i, (x, y) in enumerate(case.phases):
        add_bundle(wires, x, y, case.n_sub, case.spacing, case.r_sub, group=i, label=PHASE_NAMES[i])
    pot = dict(enumerate(three_phase_potentials(case.v_ll_kv, case.sequence)))
    if case.guard_mode not in ("grounded", "insulated", "none"):
        raise ValueError("guard_mode must be 'grounded', 'insulated' or 'none'")
    if case.guard_mode != "none":
        for j, (x, y) in enumerate(case.guards):
            wires.append(Wire(x, y, case.r_guard, group=3 + j, label=f"G{j + 1}"))
            pot[3 + j] = 0.0 if case.guard_mode == "grounded" else None
    return wires, pot


def solve_case(case: LineCase):
    """Charge simulation of the case. Returns ``(wires, solution)``."""
    wires, pot = build(case)
    return wires, solve(wires, pot, n_charges=case.n_charges)


def kuffel_gradients(case: LineCase):
    """Maximum surface gradient of each phase by the **formula** of Kuffel / CIGRE [kV/cm peak].

    The phase charge comes from the Maxwell potential matrix with the grounded shield wires kept in the
    system at 0 V (``corona.phase_charges``) and the gradient from ``corona.bundle_max_gradient_v_m``.
    An insulated shield wire is ignored (its effect on the phases is negligible, see the tests).
    """
    xy = [tuple(p) for p in case.phases]
    r_eq = [corona.equivalent_radius(case.n_sub, case.r_sub, case.spacing)] * 3
    v = list(three_phase_potentials(case.v_ll_kv, case.sequence))
    if case.guard_mode == "grounded":
        xy += [tuple(g) for g in case.guards]
        r_eq += [case.r_guard] * len(case.guards)
        v += [0.0] * len(case.guards)
    q = corona.phase_charges(np.array(xy, float), r_eq, np.array(v))
    return [corona.bundle_max_gradient_v_m(abs(q[i]), case.n_sub, case.r_sub, case.spacing) * V_PER_M_TO_KV_PER_CM
            for i in range(3)]


def ground_profile(case: LineCase, half_width: float = None, n: int = 241, height: float = 1.0, **changes):
    """Vertical field at ``height`` above ground [kV/m rms] along x, by charge simulation. Returns ``(x, E)``.

    ``half_width`` defaults to the outermost conductor plus 15 m. Keyword ``changes`` replace case
    parameters for this call only."""
    case = case.replace(**changes) if changes else case
    _, sol = solve_case(case)
    xs = [p[0] for p in case.phases] + [g[0] for g in case.guards]
    hw = half_width or (max(abs(v) for v in xs) + 15.0)
    x = np.linspace(-hw, hw, n)
    _, ey = sol.field(x, np.full_like(x, height))
    return x, np.abs(ey) / np.sqrt(2.0) / 1e3


def analyze(case: LineCase, **changes) -> dict:
    """Results of a case, optionally with some parameters changed for this call only.

    Every key says which method produced the number: ``(simulación)`` is the charge simulation,
    ``(Kuffel)`` and ``(Peek)`` are formulas. Returns the peak surface gradient of each phase by both
    methods [kV/cm], the Peek critical gradient, the ratio ``E_max / E_c`` of the worst phase by each
    method, whether it meets 0.95, the difference between the Kuffel formula and the simulation, and
    the highest ground field of the simulation (value and its distance |x| from the line centre; the
    profile is symmetric).
    """
    case = case.replace(**changes) if changes else case
    wires, sol = solve_case(case)
    e = sol.max_gradient_by_group()
    sim = [e[k] * V_PER_M_TO_KV_PER_CM for k in range(3)]
    kuf = kuffel_gradients(case)
    ec = corona.peek_critical_gradient(case.r_sub * 100.0, case.m_surface,
                                       corona.air_density(case.altitude_m, case.temp_c))
    worst = int(np.argmax(sim))
    x, eg = ground_profile(case)
    out = {}
    for k in range(3):
        out[f"E_max {PHASE_NAMES[k]} (simulación) [kV/cm]"] = sim[k]
    for k in range(3):
        out[f"E_max {PHASE_NAMES[k]} (Kuffel) [kV/cm]"] = kuf[k]
    out.update({
        "E_c (Peek) [kV/cm]": ec,
        "peor fase (simulación)": PHASE_NAMES[worst],
        "E_max / E_c, peor fase (simulación)": sim[worst] / ec,
        "E_max / E_c, peor fase (Kuffel)": max(kuf) / ec,
        "cumple < 0,95 (simulación)": bool(sim[worst] / ec < 0.95),
        "cumple < 0,95 (Kuffel)": bool(max(kuf) / ec < 0.95),
        "Kuffel vs simulación, peor fase [%]": 100.0 * (kuf[worst] / sim[worst] - 1.0),
        "campo en el suelo, máx. (simulación) [kV/m rms]": float(eg.max()),
        "|x| del máx. en el suelo (simulación) [m]": float(abs(x[eg.argmax()])),
    })
    return out


def sweep(case: LineCase, label: str, values, make=None):
    """One ``analyze`` row per value. Returns a DataFrame indexed by the values.

    ``label`` is the column name of the varied quantity and, if ``make`` is not given, the name of the
    ``LineCase`` parameter that takes each value. ``make(value)`` can instead return a dict of
    parameters, for variations that change several of them at once, for example
    ``make=lambda h: {"guards": ((-7.0, h), (7.0, h))}``.
    """
    import pandas as pd

    rows = []
    for v in values:
        changes = make(v) if make else {label: v}
        rows.append({label: v, **analyze(case, **changes)})
    return pd.DataFrame(rows).set_index(label)


def plot_sweep(df, columns=("E_max b (simulación) [kV/cm]",), title=None, ylabel="Gradiente superficial pico [kV/cm]"):
    """Plot some columns of a ``sweep`` table against its index.

    Columns of the charge simulation are drawn with solid lines and columns of the Kuffel formula with
    dashed lines. The Peek critical gradient and the 0.95 limit are added when the table has them."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    for c in columns:
        ax.plot(df.index, df[c], "s--" if "Kuffel" in c else "o-", label=c)
    ec_col = "E_c (Peek) [kV/cm]"
    if ec_col in df and any("E_max" in c for c in columns):
        ax.axhline(df[ec_col].iloc[0], color="k", ls="--", label="E_c (fórmula de Peek)")
        ax.axhline(0.95 * df[ec_col].iloc[0], color="k", ls=":", label="0,95 E_c")
    ax.set_xlabel(df.index.name)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig
