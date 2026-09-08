# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Lumped-circuit approximation of a transmission line: a cascade of N L-C sections.

This is the model a circuits course reaches for when the field problem is set
aside.  The line of length l is chopped into N pieces, and each piece becomes a
series inductor followed by a shunt capacitor:

    Vs(t) o--[ L ]--+--[ L ]--+-- ... --[ L ]--+---o  open end
                    |         |                |
                   ===       ===              ===
                    C         C                C
                    |         |                |
    ----------------+---------+----------------+------  ground

    L = L' l/N        C = C' l/N

State equations, with v0 = Vs(t) and no current past the last node:

    L di_k/dt = v_(k-1) - v_k         k = 1 .. N
    C dv_k/dt = i_k - i_(k+1)         i_(N+1) = 0  (open end)

N = 1 is the plain lumped model: the whole line collapses to one L and one C,
and every notion of position along the line disappears with it.  Raising N
recovers propagation, but only in the limit -- which is the point the
accompanying demo makes.

The system is integrated with solve_ivp rather than a leapfrog on purpose.  A
leapfrog would add its own discretisation error on top of the modelling error,
and at N = 1 a student could then dismiss the disagreement with the FDTD result
as "just the integrator".  With an accurate integrator, what is left is the
model being wrong.
"""
import numpy as np
from scipy.integrate import solve_ivp


def solve_ladder(n_sections, L_per_m, C_per_m, length, source, t_eval,
                 load=None, end_capacitance=0.0, rtol=1e-8, atol=1e-10):
    """Integrate an N-section L-C ladder driven by an ideal voltage source.

    Parameters
    ----------
    n_sections : int
        Number of L-C sections, N.  N = 1 is the lumped approximation.
    L_per_m, C_per_m : float
        Per-unit-length inductance [H/m] and capacitance [F/m].
    length : float
        Line length [m].
    source : callable
        Ideal source voltage Vs(t) [V], zero internal impedance.
    t_eval : array_like
        Times at which to report the solution [s].
    load : float, optional
        Resistive load across the far end [ohm].  Default is an open circuit.
    end_capacitance : float, optional
        Extra shunt capacitance at the far node [F], modelling the fringing
        field that spills out of an open end.  Refining the ladder cannot
        produce this term: it is a property of the geometry at the end, not of
        the line, so its value has to come from a field solution (or from
        measurement).  See the open-end discussion in the demo.

    Returns
    -------
    (currents, voltages) : ndarray, ndarray
        Shapes (N, len(t_eval)).  currents[k] is the current through inductor
        k+1; voltages[k] is the voltage on capacitor k+1, so voltages[-1] is
        the far-end (open-circuit) voltage.
    """
    if n_sections < 1:
        raise ValueError("n_sections must be at least 1")

    n = int(n_sections)
    dl = length / n
    L = L_per_m * dl
    C = C_per_m * dl
    g_load = 0.0 if load is None else 1.0 / load

    t_eval = np.asarray(t_eval, dtype=float)

    def rhs(t, y):
        i = y[:n]
        v = y[n:]
        di = np.empty(n)
        dv = np.empty(n)

        di[0] = (source(t) - v[0]) / L
        if n > 1:
            di[1:] = (v[:-1] - v[1:]) / L
            dv[:-1] = (i[:-1] - i[1:]) / C
        dv[-1] = (i[-1] - g_load * v[-1]) / (C + end_capacitance)

        return np.concatenate((di, dv))

    sol = solve_ivp(
        rhs,
        (t_eval[0], t_eval[-1]),
        np.zeros(2 * n),
        t_eval=t_eval,
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"ladder integration failed: {sol.message}")

    return sol.y[:n], sol.y[n:]


def lossless_line_step_response(t, transit_time, amplitude=1.0):
    """Exact open-circuit far-end voltage of an ideal line driven by a step.

    The wave needs one transit time T to arrive, doubles at the open end, comes
    back, inverts at the ideal (zero-impedance) source and cancels: the result
    is a square wave alternating between 0 and 2*amplitude with period 4T.

    This is the N -> infinity limit that the ladder converges to, useful as a
    reference curve that involves no simulation at all.
    """
    t = np.asarray(t, dtype=float)
    phase = np.floor_divide(t - transit_time, 2.0 * transit_time)
    return np.where(
        t < transit_time,
        0.0,
        2.0 * amplitude * np.where(phase % 2 == 0, 1.0, 0.0),
    )
