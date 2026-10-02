# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
# Written with the help of Claude (Anthropic); see NOTICE.md.
"""Physics tests for the charge simulation, the bundle gradient and Peek's law.

Each test pins down one claim made by ``src/electric_field``:

  1. the simulation reproduces the exact single-cylinder-over-ground solution,
  2. the boundary condition holds on the conductor surface (small error, even n_charges too),
  3. Kuffel eq. (2.9) agrees with the simulation for 2-, 3- and 4-bundles,
  4. a grounded shield wire is at 0 V, carries induced charge and raises the phase gradient,
  5. an insulated wire has net charge zero and floats between 0 V and the phase potentials,
  6. the Maxwell-matrix charges with the shield wires kept in the system match the simulation,
  7. Peek's law and the air density reproduce hand-computed values.
"""
import numpy as np
import pytest

import corona
from charge_simulation import (EPS0, K_COULOMB, Wire, add_bundle, charge_per_group, solve,
                               three_phase_potentials)

V = 1.0e5   # arbitrary test potential [V]


def line_500kv(guard):
    """Three 4-bundles at 25 m and two shield wires at 35 m. ``guard``: 0, None or 'absent'."""
    wires = []
    for i, x in enumerate([-11.0, 0.0, 11.0]):
        add_bundle(wires, x, 25.0, 4, 0.45, 0.0141, group=i, label="abc"[i])
    pot = dict(enumerate(three_phase_potentials(500.0)))
    if guard != "absent":
        for j, x in enumerate([-7.0, 7.0]):
            wires.append(Wire(x, 35.0, 0.0055, group=3 + j, label=f"G{j}"))
            pot[3 + j] = guard
    return wires, pot


# --- 1 and 2: the solver ------------------------------------------------------------------------

def test_single_cylinder_charge_matches_exact_solution():
    h, r = 15.0, 0.0125
    sol = solve([Wire(0.0, h, r, 0)], {0: V})
    q_exact = 2 * np.pi * EPS0 * V / np.arccosh(h / r)      # exact for a cylinder over ground
    assert charge_per_group(sol)[0].real == pytest.approx(q_exact, rel=1e-4)


@pytest.mark.parametrize("n_charges", [6, 12, 24])
def test_boundary_condition_holds_also_for_even_number_of_charges(n_charges):
    # Regression: interleaving control points and charges made P singular for even n_charges.
    wires = add_bundle([], 0.0, 25.0, 4, 0.45, 0.0141, group=0)
    sol = solve(wires, {0: V}, n_charges=n_charges)
    assert sol.bc_error < 1e-3


def test_gradient_is_linear_in_voltage():
    wires = add_bundle([], 0.0, 25.0, 3, 0.4, 0.0141, group=0)
    e1 = solve(wires, {0: V}).max_gradient_by_group()[0]
    e2 = solve(wires, {0: 2 * V}).max_gradient_by_group()[0]
    assert e2 == pytest.approx(2 * e1, rel=1e-9)


# --- 3: Kuffel eq. (2.9) ------------------------------------------------------------------------

@pytest.mark.parametrize("n, spacing", [(2, 0.40), (3, 0.40), (4, 0.45)])
def test_kuffel_bundle_gradient_matches_simulation(n, spacing):
    r = 0.0141
    wires = add_bundle([], 0.0, 25.0, n, spacing, r, group=0)
    sol = solve(wires, {0: V}, n_charges=24)
    q = charge_per_group(sol)[0].real
    e_sim = sol.max_gradient_by_group()[0]
    e_kuffel = corona.bundle_max_gradient_v_m(q, n, r, spacing)
    assert e_kuffel == pytest.approx(e_sim, rel=0.01)


def test_bundle_max_is_above_average_by_the_kuffel_bracket():
    q, n, r, a = 1e-6, 4, 0.0141, 0.45
    R = corona.bundle_radius(n, a)
    ratio = corona.bundle_max_gradient_v_m(q, n, r, a) / corona.bundle_average_gradient_v_m(q, n, r)
    assert ratio == pytest.approx(1 + (n - 1) * r / R)
    assert R == pytest.approx(a / np.sqrt(2))          # square bundle: R = a / sqrt(2)


# --- 4 and 5: shield wires ----------------------------------------------------------------------

def test_grounded_shield_wire_is_at_zero_volts_and_collects_charge():
    wires, pot = line_500kv(guard=0.0)
    sol = solve(wires, pot)
    for g in (3, 4):
        assert abs(charge_per_group(sol)[g]) > 1e-8          # induced charge, not ignored
    # Potential on the surface of the shield wires (check at points off the control points)
    for w in wires[12:]:
        phi = np.linspace(0, 2 * np.pi, 37)
        v = sol.potential(w.x + w.r * np.cos(phi), w.y + w.r * np.sin(phi))
        assert np.abs(v).max() < 1e-4 * np.abs(pot[0])


def test_grounded_shield_wire_raises_phase_gradient_slightly():
    e_none = solve(*line_500kv("absent")).max_gradient_by_group()[1]
    e_gnd = solve(*line_500kv(0.0)).max_gradient_by_group()[1]
    assert 1.0 < e_gnd / e_none < 1.02


def test_insulated_shield_wire_has_zero_net_charge_and_floats_between_0_and_phase_peak():
    wires, pot = line_500kv(guard=None)
    sol = solve(wires, pot)
    assert sol.floating == (3, 4)
    for g in (3, 4):
        assert abs(charge_per_group(sol)[g]) < 1e-9 * abs(charge_per_group(sol)[1])
        vg = abs(sol.group_potentials[g])
        assert 0.0 < vg < abs(pot[0])


def test_insulated_wire_has_almost_no_effect_on_phase_gradient():
    e_none = solve(*line_500kv("absent")).max_gradient_by_group()[1]
    e_ins = solve(*line_500kv(None)).max_gradient_by_group()[1]
    assert e_ins == pytest.approx(e_none, rel=2e-3)


# --- 6: Maxwell matrix with the shield wires kept in the system ---------------------------------

def test_maxwell_charges_with_grounded_guard_match_simulation():
    wires, pot = line_500kv(guard=0.0)
    sol = solve(wires, pot, n_charges=12)
    q_sim = charge_per_group(sol)
    xy = np.array([[-11, 25], [0, 25], [11, 25], [-7, 35], [7, 35]], float)
    r_eq = [corona.equivalent_radius(4, 0.0141, 0.45)] * 3 + [0.0055] * 2
    v = np.array([pot[0], pot[1], pot[2], 0.0, 0.0])
    q = corona.phase_charges(xy, r_eq, v)
    for i in range(5):
        assert abs(q[i]) == pytest.approx(abs(q_sim[i]), rel=2e-3)


def test_dropping_the_guard_from_the_matrix_gives_a_different_charge():
    # Inverting only the phase block P_pp ignores the shield wire: the strategy is to keep it.
    xy = np.array([[-11, 25], [0, 25], [11, 25], [-7, 35], [7, 35]], float)
    r_eq = np.array([corona.equivalent_radius(4, 0.0141, 0.45)] * 3 + [0.0055] * 2)
    vp = three_phase_potentials(500.0)
    q_full = corona.phase_charges(xy, r_eq, np.r_[vp, 0.0, 0.0])[:3]
    q_phases_only = corona.phase_charges(xy[:3], r_eq[:3], vp)
    assert np.abs(q_full[1]) > np.abs(q_phases_only[1])


# --- 7: Peek and air density --------------------------------------------------------------------

def test_air_density_hand_values():
    assert corona.air_density(0.0, 25.0) == pytest.approx(0.9997, abs=5e-4)
    b = 76 * 10 ** (-2600 / 18336)
    assert corona.air_density(2600.0, 15.0) == pytest.approx(3.92 * b / 288.0)
    assert corona.air_density(2600.0, 15.0) < corona.air_density(0.0, 15.0)


def test_peek_critical_gradient_hand_values():
    # r = 1 cm, delta = 1, m = 0.8, stranded: 30 * 0.8 * (1 + 0.308) = 31.392 kV/cm peak
    assert corona.peek_critical_gradient(1.0, 0.8, 1.0) == pytest.approx(31.392)
    assert corona.peek_critical_gradient(1.0, 1.0, 1.0, stranded=False) == pytest.approx(30 * 1.301)
    # E_c grows as the radius shrinks (thin wires tolerate a higher surface gradient)
    assert corona.peek_critical_gradient(0.5, 0.8, 1.0) > corona.peek_critical_gradient(2.0, 0.8, 1.0)
    assert corona.rms_from_peak(30.0) == pytest.approx(21.213, abs=1e-3)


def test_lower_air_density_lowers_critical_gradient():
    assert (corona.peek_critical_gradient(1.41, 0.8, corona.air_density(2600, 15))
            < corona.peek_critical_gradient(1.41, 0.8, corona.air_density(0, 15)))


def test_margin_criterion():
    ratio, ok = corona.gradient_margin(20.0, 25.0)
    assert ratio == pytest.approx(0.8) and ok
    assert not corona.gradient_margin(24.0, 25.0)[1]
    assert corona.onset_voltage_ratio(20.0, 25.0) == pytest.approx(1.25)


def test_emax_bundle_script_keeps_the_original_hand_calculation():
    import Emax_bundle
    ei = Emax_bundle.voltage_gradient(Ci=1e-11, n2=1, r=0.01, s=0.457, U=115)
    old = 1e-11 / (2 * np.pi * 8.854e-12 * 1 * 0.01) * (115 / (np.sqrt(3) * 100))   # the old formula
    assert ei == pytest.approx(old, rel=1e-3)


# --- 8: empirical RI / AN formulas (hand-computed anchors; constants unverified, see corona.py) ------

def test_radio_interference_reference_point_and_gradient_slope():
    # At the reference values of the formula (g0, d0, f = 1 MHz, D = 19.9 m, sea level) only 51.7 remains
    assert corona.radio_interference_db(25.6, 4.62, 19.9) == pytest.approx(51.7)
    # Doubling g adds 86 log10(2) dB
    assert (corona.radio_interference_db(51.2, 4.62, 19.9) - 51.7) == pytest.approx(86 * np.log10(2))
    # 40 log10 decay with distance: doubling D removes 40 log10(2) dB
    assert (corona.radio_interference_db(25.6, 4.62, 39.8) - 51.7) == pytest.approx(-40 * np.log10(2))


def test_audible_noise_constants_switch_at_three_subconductors():
    assert corona.audible_noise_dba(1.0, 2, 1.0, 1.0) == pytest.approx(-93.40)
    assert corona.audible_noise_dba(1.0, 4, 1.0, 1.0) == pytest.approx(-100.62 + 25.6 * np.log10(4))
    assert (corona.audible_noise_dba(1.0, 2, 1.0, 10.0) - corona.audible_noise_dba(1.0, 2, 1.0, 1.0)
            == pytest.approx(-11.4))


def test_cispr_20m_formula():
    assert corona.ri_cispr_20m_db(15.0, 1.41) == pytest.approx(3.5 * 15 + 12 * 1.41 - 30)


# --- 9: time domain (one 60 Hz cycle) ------------------------------------------------------------

def _line_solution():
    wires, pot = line_500kv(guard=0.0)
    return wires, pot, solve(wires, pot, n_charges=12)


@pytest.mark.parametrize("k", [0, 1, 2])
def test_bundle_potentials_follow_the_three_phase_sequence(k):
    # At wt = 120 deg * k, phase a, b, c is at +Vp in turn and the other two are at -Vp/2 (positive sequence)
    wires, pot, sol = _line_solution()
    vp = np.sqrt(2 / 3) * 500e3
    wt = 2 * np.pi * k / 3
    measured = []
    for g in range(3):
        w = next(x for x in wires if x.group == g)
        measured.append(float(sol.potential_at(w.x + w.r, w.y, wt)))
    expected = [-vp / 2] * 3
    expected[k] = vp
    assert measured == pytest.approx(expected, rel=2e-4)


def test_field_is_periodic_and_half_wave_antisymmetric():
    _, _, sol = _line_solution()
    x, y = np.array([-3.0, 4.0]), np.array([20.0, 30.0])
    ex0, ey0 = sol.field_at(x, y, 0.7)
    ex1, ey1 = sol.field_at(x, y, 0.7 + np.pi)
    ex2, ey2 = sol.field_at(x, y, 0.7 + 2 * np.pi)
    assert ex1 == pytest.approx(-ex0) and ey1 == pytest.approx(-ey0)
    assert ex2 == pytest.approx(ex0) and ey2 == pytest.approx(ey0)


def test_peak_of_the_instantaneous_gradient_equals_the_phasor_amplitude():
    _, _, sol = _line_solution()
    wts, res = sol.max_gradient_vs_time(n_times=180, n_points=180)
    amp = sol.max_gradient_by_group(180)
    for g in range(3):
        assert res[g].max() == pytest.approx(amp[g], rel=1e-3)


def test_each_phase_reaches_its_gradient_peak_at_a_different_instant():
    _, _, sol = _line_solution()
    wts, res = sol.max_gradient_vs_time(n_times=180, n_points=90)
    # the centre phase (b) and the side phase (c) peak at different angles; |E| peaks every 180 deg
    db = np.degrees(wts[res[1].argmax()]) % 180
    dc = np.degrees(wts[res[2].argmax()]) % 180
    assert abs(db - dc) > 30


def test_bundle_zoom_draws_equipotentials_and_field_lines():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from charge_simulation import draw_bundle_zoom
    _, _, sol = _line_solution()
    fig, ax = plt.subplots()
    draw_bundle_zoom(ax, sol, group=1, wt=0.0)
    assert len(ax.patches) >= 4                      # the four sub-conductors
    assert len(ax.collections) > 0 and len(ax.patches) >= 4
    plt.close(fig)


# --- 10: base case of Nolasco et al. (CIGRE Green Books, Overhead Lines, sec. 4.19) --------------------
# +-500 kV bipolar HVDC, 3 x Lapwing: d = 3.822 cm, bundle spacing 45 cm, pole spacing 13.0 m,
# minimum clearance 12.5 m. The numbers below are the ones printed in eqs. 4.212-4.215, 4.221 and
# the loss, RI and AN examples of the same section.

BOOK = dict(v_kv=500.0, n=3, r_cm=1.911, a_cm=45.0, h_cm=1250.0, s_cm=1300.0)


def test_book_bundle_geometry_eqs_4_210_to_4_213():
    assert corona.bundle_radius(3, 45.0) == pytest.approx(26.0, abs=0.05)                       # eq. 4.212
    assert corona.equivalent_radius(3, 1.911, 45.0) == pytest.approx(15.7, abs=0.05)             # eq. 4.213


def test_book_average_and_maximum_gradient_eqs_4_214_and_4_215():
    ea = corona.markt_mengele_average_gradient_kv_cm(**BOOK)
    assert ea == pytest.approx(20.297, abs=2e-3)                                                 # eq. 4.214
    em = corona.bundle_max_from_average(ea, 3, 1.911, 45.0)
    assert em == pytest.approx(23.28, abs=0.01)                                                  # eq. 4.215


def test_book_critical_gradient_eq_4_221():
    # The book uses the smooth-cylinder constant 0.301 with m = 0.82 and delta = 0.915
    ec = corona.peek_critical_gradient(1.911, m=0.82, delta=0.915, stranded=False)
    assert ec == pytest.approx(27.6, abs=0.05)


def test_book_corona_loss_fair_foul_and_weighted():
    fair = corona.hvdc_bipolar_corona_loss_w_m(23.28, 3.822, 3, 12.5, 13.0, "fair")
    foul = corona.hvdc_bipolar_corona_loss_w_m(23.28, 3.822, 3, 12.5, 13.0, "foul")
    assert fair == pytest.approx(3.7, abs=0.05)
    assert foul == pytest.approx(20.6, abs=0.1)
    assert 0.8 * fair + 0.2 * foul == pytest.approx(7.1, abs=0.05)


def test_book_audible_noise_and_radio_interference():
    radial = np.hypot(12.5, 30.0 - 6.5)                                  # 26.6 m from the positive pole
    assert corona.audible_noise_dba(23.28, 3, 3.822, radial, altitude_m=600.0) == pytest.approx(38.2, abs=0.05)
    # The printed expression shows D = 30 m, which gives 39.7 dB; the printed result (41.8 dB) matches
    # the radial distance used in the AN example.
    ri = corona.radio_interference_db(23.28, 3.822, radial, freq_mhz=1.0, altitude_m=600.0)
    assert ri == pytest.approx(41.8, abs=0.1)
    assert corona.radio_interference_db(23.28, 3.822, 30.0, 1.0, 600.0) == pytest.approx(39.7, abs=0.05)


def test_charge_simulation_agrees_with_the_book_equations_for_the_base_case():
    wires = []
    add_bundle(wires, -6.5, 12.5, 3, 0.45, 0.01911, group=0)
    add_bundle(wires, 6.5, 12.5, 3, 0.45, 0.01911, group=1)
    sol = solve(wires, {0: 500e3, 1: -500e3}, n_charges=16)
    q = abs(charge_per_group(sol)[0])
    e_avg = corona.bundle_average_gradient_v_m(q, 3, 0.01911) * 1e-5
    e_max = sol.max_gradient_by_group()[0] * 1e-5
    assert e_avg == pytest.approx(20.297, rel=2e-3)          # average gradient, eq. 4.214
    assert e_max == pytest.approx(23.28, rel=0.02)           # the book says the equations are "reasonably accurate"
    assert e_max > 23.28                                      # the simulation includes the other pole and the ground


def test_shield_wire_gradient_eq_4_220_matches_the_simulation():
    # E_sw = Q_sw / (2 pi eps r): the average gradient of a grounded shield wire from its induced charge
    wires, pot = line_500kv(guard=0.0)
    sol = solve(wires, pot, n_charges=16)
    q = abs(charge_per_group(sol)[3])
    e_avg = q * K_COULOMB / 0.0055
    e_max = sol.max_gradient_by_group()[3]
    assert e_avg == pytest.approx(e_max, rel=0.05)
