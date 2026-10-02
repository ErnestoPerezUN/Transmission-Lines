# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
# First version generated with ChatGPT (Kuffel, Zaengl, Kuffel, eq. 2.9); rewritten with the help of
# Claude (Anthropic) so that the physics lives in src/electric_field/corona.py. See NOTICE.md.
"""Surface gradient of a conductor bundle and corona margin for a 3-phase line.

Run it from anywhere:  ``python notebooks/Emax_TL/Emax_bundle.py``

What this script teaches
------------------------
1. ``voltage_gradient`` is the classic hand calculation (Kuffel eq. 2.9): if you know the
   capacitance per unit length of the phase, the charge is ``q = C * V_phase`` and the highest
   surface gradient follows from the bundle geometry.
2. The capacitance is not a free input. It comes from the geometry of the whole line, and it is
   affected by the other phases and by the shield wires. ``phase_charges`` (in ``corona.py``)
   computes it with the shield wires kept in the system at 0 V.
3. The gradient is compared with Peek's critical gradient at the site altitude.
   The design criterion is E_max / E_c < 0.95.

The gradient formula itself is not repeated here: it lives in ``corona.py`` (one place, one source).
"""
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "src" / "electric_field"))

import numpy as np                                    # noqa: E402
import corona                                         # noqa: E402
from charge_simulation import three_phase_potentials  # noqa: E402


def voltage_gradient(Ci, n2, r, s, U):
    """Maximum surface gradient of a bundle, in kV/cm (rms), from the phase capacitance.

    Parameters
    ----------
    Ci : float
        Capacitance per unit length of the phase to ground [F/m]. It already includes the effect of
        the other phases and of the shield wires.
    n2 : int
        Number of sub-conductors per bundle.
    r : float
        Radius of one sub-conductor [m]. (Radius, not diameter.)
    s : float
        Distance between adjacent sub-conductors [m].
    U : float
        Line-to-line rms voltage [kV]. The phase-to-ground voltage is U / sqrt(3).

    Returns
    -------
    float
        Maximum surface gradient [kV/cm], rms. Multiply by sqrt(2) for the peak value that is
        compared with Peek's critical gradient.
    """
    q_rms = Ci * U * 1e3 / math.sqrt(3.0)             # line charge [C/m], rms
    return corona.bundle_max_gradient_v_m(q_rms, n2, r, s) * 1e-5


def demo_230kv_line():
    """230 kV horizontal line with a 2-bundle and two shield wires, at 1000 m altitude."""
    v_ll = 230.0                     # kV, line to line, rms
    n, r, s = 2, 0.0141, 0.40        # ACSR Bluejay-like sub-conductor: r = 14.1 mm, 40 cm spacing
    xy = np.array([[-6.5, 20.0], [0.0, 20.0], [6.5, 20.0],     # phases a, b, c
                   [-4.0, 27.0], [4.0, 27.0]])                 # two shield wires
    r_eq = [corona.equivalent_radius(n, r, s)] * 3 + [0.0055] * 2
    v = np.r_[three_phase_potentials(v_ll), 0.0, 0.0]          # shield wires at 0 V
    q = corona.phase_charges(xy, r_eq, v)

    delta = corona.air_density(1000.0, 20.0)
    e_c = corona.peek_critical_gradient(r * 100, m=0.8, delta=delta)
    print(f"Relative air density at 1000 m, 20 degC: {delta:.3f}")
    print(f"Peek critical gradient (peak): {e_c:.1f} kV/cm")
    for name, qi in zip("abc", q[:3]):
        e_max = corona.bundle_max_gradient_v_m(qi, n, r, s) * 1e-5
        ratio, ok = corona.gradient_margin(e_max, e_c)
        print(f"Phase {name}: E_max = {e_max:5.2f} kV/cm (peak), E_max/E_c = {ratio:.2f}  "
              f"{'OK' if ok else 'above 0.95: corona expected'}")


if __name__ == "__main__":
    # The original example of this file: single conductor, 115 kV, Ci = 10 pF/m
    Ei = voltage_gradient(Ci=1e-11, n2=1, r=0.01, s=0.457, U=115)
    print(f"Voltage gradient Ei = {Ei:.4f} kV/cm (rms), Kuffel formula")
    print()
    demo_230kv_line()
