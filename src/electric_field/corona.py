# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
# Written with the help of Claude (Anthropic); see NOTICE.md.
"""Corona inception on overhead lines: critical gradient, bundle gradient and design margin.

Corona starts when the electric field at the conductor surface exceeds the value at which the air
around it ionises. The design question is therefore a comparison of two numbers:

* ``E_max``: the highest surface gradient of the conductor at operating voltage. It depends on the
  geometry only (heights, phase spacing, bundle) and grows linearly with voltage.
* ``E_c``: the critical gradient of the air (Peek's law). It depends on the air density
  (altitude and temperature), on the conductor radius and on the surface condition.

The design criterion used in the course is ``E_max / E_c < 0.95``.

Everything here uses **peak** values, because Peek's 30 kV/cm is a peak value (21.2 kV/cm rms).

Units
-----
* Radii and spacings passed in metres unless the name ends in ``_cm``.
* Line charges in C/m, potentials in volts (complex peak phasors when three-phase).
* Gradients in V/m from the ``*_v_m`` functions and in kV/cm from the ``*_kv_cm`` functions.

Sources
-------
* Peek's law: F. W. Peek, *Dielectric Phenomena in High Voltage Engineering* (1929). Written here
  in its usual form for stranded conductors (surface constant 0.308; 0.301 for a smooth
  cylinder), as used in the course notes (Cardona Correa, 2022) and CIGRE TB 638 section 2.6.
* Bundle gradient: E. Kuffel, W. S. Zaengl, J. Kuffel, *High Voltage Engineering Fundamentals*,
  eq. (2.9), the same formula that ``notebooks/Emax_TL/Emax_bundle.py`` used before.
* Air density versus altitude: Peek's barometric relation ``b = 76 * 10**(-h/18336)`` cmHg.

Everything in this module is a **formula** evaluated directly. The phase charges use the Maxwell
potential coefficients with the image method, with the ground wires at zero potential (see
``phase_charges``). ``charge_simulation.py`` solves the same problem with the **charge simulation**
method and the exact conductor geometry, and is the reference for the tests. Results shown to
students always say which of the two produced them.
"""
from __future__ import annotations

import numpy as np

EPS0 = 8.8541878128e-12          # vacuum permittivity [F/m]
K_COULOMB = 1.0 / (2.0 * np.pi * EPS0)


# ---------------------------------------------------------------------------
# Air and critical gradient
# ---------------------------------------------------------------------------

def air_density(altitude_m: float, temp_c: float = 25.0) -> float:
    """Relative air density ``delta`` (1.0 at 25 degC and 76 cmHg).

    ``delta = 3.92 * b / (273 + T)`` with the barometric pressure ``b = 76 * 10**(-h / 18336)``
    in cmHg (Peek). At sea level and 25 degC it gives 0.9997; at 2600 m and 15 degC, 0.746.
    """
    b = 76.0 * 10.0 ** (-altitude_m / 18336.0)
    return 3.92 * b / (273.0 + temp_c)


def peek_critical_gradient(r_cm: float, m: float = 0.8, delta: float = 1.0,
                           stranded: bool = True) -> float:
    """Peek critical (inception) gradient, **peak** value in kV/cm.

    ``E_c = 30 * m * delta * (1 + k / sqrt(delta * r))`` with ``r`` in cm and ``k`` = 0.308 for
    a stranded conductor (0.301 for a smooth cylinder).

    ``m`` is the surface factor: 1.0 smooth clean cylinder, 0.8 to 0.9 clean stranded conductor,
    0.6 to 0.7 weathered or dusty, 0.3 to 0.6 with rain drops (the drops act as sharp points).
    ``r_cm`` is the radius of one **sub-conductor**, not of the bundle.
    """
    k = 0.308 if stranded else 0.301
    return 30.0 * m * delta * (1.0 + k / np.sqrt(delta * r_cm))


def rms_from_peak(x_peak):
    """Peak to rms of a sinusoid."""
    return x_peak / np.sqrt(2.0)


# ---------------------------------------------------------------------------
# Bundle geometry and gradient
# ---------------------------------------------------------------------------

def bundle_radius(n: int, spacing: float) -> float:
    """Radius ``R`` of the circle through the centres of ``n`` sub-conductors.

    ``R = spacing / (2 sin(pi/n))`` where ``spacing`` is the distance between adjacent
    sub-conductors. For ``n = 1`` it returns 0."""
    if n == 1:
        return 0.0
    return spacing / (2.0 * np.sin(np.pi / n))


def equivalent_radius(n: int, r: float, spacing: float) -> float:
    """Equivalent radius ``r_eq = R * (n r / R)**(1/n)`` of the bundle (geometric mean distance).

    A bundle of ``n`` sub-conductors of radius ``r`` has the same capacitance to the far world as
    a single conductor of radius ``r_eq``. For ``n = 1`` it is just ``r``."""
    if n == 1:
        return r
    R = bundle_radius(n, spacing)
    return R * (n * r / R) ** (1.0 / n)


def bundle_average_gradient_v_m(q: float, n: int, r: float) -> float:
    """Average surface gradient of a sub-conductor, ``q / (2 pi eps0 n r)`` [V/m].

    ``q`` is the total line charge of the phase [C/m] (all sub-conductors together)."""
    return abs(q) * K_COULOMB / (n * r)


def bundle_max_gradient_v_m(q: float, n: int, r: float, spacing: float) -> float:
    """Maximum surface gradient of a bundle, Kuffel eq. (2.9) [V/m].

    ``E_max = q / (2 pi eps0 n r) * [1 + (n - 1) r / R]`` with ``R = spacing / (2 sin(pi/n))``.
    The bracket is the non-uniformity: the field is highest on the outer side of a sub-conductor
    (away from the bundle centre) and lowest on the inner side, and the average is the term
    before the bracket. It assumes the field of the other phases is negligible.
    """
    if n == 1:
        return bundle_average_gradient_v_m(q, n, r)
    R = bundle_radius(n, spacing)
    return bundle_average_gradient_v_m(q, n, r) * (1.0 + (n - 1) * r / R)


# ---------------------------------------------------------------------------
# Phase charges with ground wires at zero potential
# ---------------------------------------------------------------------------

def phase_charges(xy: np.ndarray, r_eq: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Line charges [C/m] of every conductor when some of them are at a prescribed potential.

    Parameters
    ----------
    xy : (N, 2) array with the (x, y) of every conductor **centre** [m] (bundle centre for a
        bundle), phases first and ground wires last.
    r_eq : (N,) equivalent radius of each conductor [m] (see ``equivalent_radius``).
    v : (N,) potentials [V]; **use 0 for a grounded shield wire**.

    Strategy for the ground wires
    -----------------------------
    The potential of every conductor is ``V = P q`` with ``P_ii = ln(2 y_i / r_i) / (2 pi eps0)``
    and ``P_ij = ln(D'_ij / D_ij) / (2 pi eps0)`` (``D`` the distance between conductors and
    ``D'`` the distance from one to the image of the other). The ground wires are **not**
    dropped: they stay in ``P`` with ``V = 0`` and ``q = P^-1 V`` is solved for all conductors.
    The charges of the phases are then the phase rows of ``P^-1 V``, which is the same as using
    the reduced capacitance matrix ``C_pp = (P_pp - P_pg P_gg^-1 P_gp)^-1``. Inverting only the
    phase block ``P_pp`` would ignore the shield wire and overestimate the gradient.
    """
    xy = np.asarray(xy, float)
    n = len(xy)
    dx = xy[:, 0][:, None] - xy[:, 0][None, :]
    d = np.hypot(dx, xy[:, 1][:, None] - xy[:, 1][None, :])
    dp = np.hypot(dx, xy[:, 1][:, None] + xy[:, 1][None, :])
    P = np.zeros((n, n))
    off = ~np.eye(n, dtype=bool)
    P[off] = np.log(dp[off] / d[off])
    P[np.eye(n, dtype=bool)] = np.log(2.0 * xy[:, 1] / np.asarray(r_eq, float))
    P *= K_COULOMB
    return np.linalg.solve(P, np.asarray(v))


# ---------------------------------------------------------------------------
# Design check
# ---------------------------------------------------------------------------

def gradient_margin(e_max_kv_cm: float, e_crit_kv_cm: float, limit: float = 0.95) -> tuple:
    """Return ``(ratio, ok)`` with ``ratio = E_max / E_c`` and ``ok = ratio < limit``.

    Both gradients must be in the same units and both peak (or both rms)."""
    ratio = e_max_kv_cm / e_crit_kv_cm
    return ratio, bool(ratio < limit)


def onset_voltage_ratio(e_max_kv_cm: float, e_crit_kv_cm: float) -> float:
    """Factor by which the line voltage can rise before ``E_max`` reaches ``E_c``.

    The gradient is proportional to voltage, so the onset voltage is ``V_operating * E_c / E_max``."""
    return e_crit_kv_cm / e_max_kv_cm


# ---------------------------------------------------------------------------
# Nolasco et al. (2017), sec. 4.19: average bundle gradient, corona loss, radio interference, audible noise
# ---------------------------------------------------------------------------
# Source: J. F. Nolasco et al., chapter "Electrical Design" of CIGRE Green Books, *Overhead Lines*
# (Springer, 2017), sec. 4.19, eqs. 4.208-4.211 and 4.223-4.226. The base case of that section
# (+-500 kV bipolar HVDC, 3 x Lapwing) reproduces in tests/test_corona.py.
#
# WARNING: the loss, RI and AN equations are empirical formulas for BIPOLAR HVDC lines (the book says
# so explicitly). They are not valid for an AC line. ``g`` is the maximum bundle gradient in kV/cm.

def markt_mengele_average_gradient_kv_cm(v_kv: float, n: int, r_cm: float, h_cm: float, s_cm: float,
                                         a_cm: float) -> float:
    """Average bundle gradient of a bipolar line by Markt and Mengele's method [kV/cm], eq. 4.208.

    ``E_a = V / (n r ln( 2H / (r_eq sqrt((2H/S)^2 + 1)) ))`` with ``V`` the pole voltage [kV], ``r`` the
    sub-conductor radius, ``H`` the conductor height, ``S`` the pole spacing, all in cm, and ``r_eq`` the
    equivalent bundle radius (eqs. 4.210-4.211). The maximum gradient is
    ``E_m = E_a (1 + (n-1) r/R)`` (eq. 4.209, see ``bundle_max_from_average``).
    """
    r_eq = equivalent_radius(n, r_cm, a_cm)
    return v_kv / (n * r_cm * np.log(2.0 * h_cm / (r_eq * np.sqrt((2.0 * h_cm / s_cm) ** 2 + 1.0))))


def bundle_max_from_average(e_avg, n: int, r: float, spacing: float):
    """Maximum bundle gradient from the average one, eq. 4.209: ``E_m = E_a (1 + (n-1) r/R)``.

    ``r`` and ``spacing`` in the same unit; the result has the unit of ``e_avg``."""
    if n == 1:
        return e_avg
    return e_avg * (1.0 + (n - 1) * r / bundle_radius(n, spacing))


def hvdc_bipolar_corona_loss_w_m(g_kv_cm: float, d_cm: float, n: int, h_m: float, s_m: float,
                                 weather: str = "fair") -> float:
    """Corona loss of a bipolar HVDC line [W/m], eqs. 4.223 and 4.224.

    ``P[dB] = P0 + k_g log(g/25) + k_d log(d/3.05) + k_n log(n/3) - 10 log(H S / (15 15))`` and
    ``P[W/m] = 10**(P[dB]/10)``, with (P0, k_g, k_d, k_n) = (2.9, 50, 30, 20) in fair weather and
    (11, 40, 20, 15) in foul weather. ``d`` is the sub-conductor diameter [cm], ``H`` the height and
    ``S`` the pole spacing [m]. The book weights 80 % fair and 20 % foul weather in the economic study.
    """
    if weather == "fair":
        p0, kg, kd, kn = 2.9, 50.0, 30.0, 20.0
    elif weather == "foul":
        p0, kg, kd, kn = 11.0, 40.0, 20.0, 15.0
    else:
        raise ValueError("weather must be 'fair' or 'foul'")
    db = (p0 + kg * np.log10(g_kv_cm / 25.0) + kd * np.log10(d_cm / 3.05) + kn * np.log10(n / 3.0)
          - 10.0 * np.log10(h_m * s_m / (15.0 * 15.0)))
    return 10.0 ** (db / 10.0)


def radio_interference_db(g_kv_cm: float, d_cm: float, dist_m: float, freq_mhz: float = 1.0,
                          altitude_m: float = 0.0) -> float:
    """Fair-weather radio interference of a bipolar HVDC line [dB(uV/m)], eq. 4.225.

    ``dist_m`` is the radial distance from the positive pole to the observation point (not the lateral
    distance: in the book's base case it is sqrt(12.5^2 + 23.5^2) = 26.6 m, although the printed
    expression shows 30 m).

    ``RI = 51.7 + 86 log10(g/25.6) + 40 log10(d/4.62) + 10 (1 - log10(10 f)^2) + 40 log10(19.9/D) + q/300``
    with ``d`` the sub-conductor diameter [cm], ``f`` the frequency [MHz] and ``q`` the altitude [m].
    """
    return (51.7 + 86.0 * np.log10(g_kv_cm / 25.6) + 40.0 * np.log10(d_cm / 4.62)
            + 10.0 * (1.0 - np.log10(10.0 * freq_mhz) ** 2) + 40.0 * np.log10(19.9 / dist_m)
            + altitude_m / 300.0)


def audible_noise_dba(g_kv_cm: float, n: int, d_cm: float, dist_m: float, altitude_m: float = 0.0) -> float:
    """Fair-weather audible noise of a bipolar HVDC line [dBA], eq. 4.226.

    ``dist_m`` is the radial distance from the positive pole to the observation point.

    ``AN = AN0 + 86 log10(g) + k log10(n) + 40 log10(d) - 11.4 log10(R) + q/300`` with
    (k, AN0) = (0, -93.40) for n <= 2 and (25.6, -100.62) for n >= 3.
    """
    k, an0 = (0.0, -93.40) if n <= 2 else (25.6, -100.62)
    return (an0 + 86.0 * np.log10(g_kv_cm) + k * np.log10(n) + 40.0 * np.log10(d_cm)
            - 11.4 * np.log10(dist_m) + altitude_m / 300.0)


def ri_cispr_20m_db(g_kv_cm: float, r_cm: float) -> float:
    """Simplified radio interference at 20 m [dB(uV/m)]: ``3.5 g + 12 r - 30``.

    Taken from the course summary (attributed there to CISPR TR 18-3); not checked against that source."""
    return 3.5 * g_kv_cm + 12.0 * r_cm - 30.0
