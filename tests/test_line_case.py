# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
# Written with the help of Claude (Anthropic); see NOTICE.md.
"""Tests of the workshop helper ``line_case``: variations must be physical, must not leak, and every
result must say which method (charge simulation or formula) produced it."""
import numpy as np
import pytest

import corona
from line_case import LineCase, analyze, build, ground_profile, kuffel_gradients, solve_case, sweep

BASE = LineCase()
SIM_B = "E_max b (simulación) [kV/cm]"
KUF_B = "E_max b (Kuffel) [kV/cm]"


def test_a_variation_does_not_change_the_base_case():
    before = analyze(BASE)
    analyze(BASE, n_sub=5, spacing=0.3)
    assert BASE.n_sub == 3 and BASE.spacing == 0.45
    assert analyze(BASE) == before


def test_unknown_parameter_is_rejected():
    with pytest.raises(TypeError):
        analyze(BASE, number_of_conductors=4)


def test_every_result_names_its_method():
    keys = analyze(BASE).keys()
    for k in keys:
        if k.startswith("E_max") or "campo en el suelo" in k or "|x|" in k:
            assert "(simulación)" in k or "(Kuffel)" in k, k
        if k.startswith("E_c"):
            assert "(Peek)" in k, k


def test_charge_simulation_gradient_agrees_with_the_kuffel_formula_for_the_centre_phase():
    from charge_simulation import charge_per_group
    _, sol = solve_case(BASE)
    q = abs(charge_per_group(sol)[1])
    e_kuffel_with_simulated_charge = corona.bundle_max_gradient_v_m(q, 3, 0.0141, 0.45) * 1e-5
    res = analyze(BASE)
    assert res[SIM_B] == pytest.approx(e_kuffel_with_simulated_charge, rel=0.01)
    # the formula fed with the Maxwell-matrix charge (no simulation at all) is also close
    assert res[KUF_B] == pytest.approx(res[SIM_B], rel=0.01)
    assert abs(res["Kuffel vs simulación, peor fase [%]"]) < 1.0


def test_kuffel_formula_ignores_an_insulated_guard_and_includes_a_grounded_one():
    grounded = kuffel_gradients(BASE)
    none = kuffel_gradients(BASE.replace(guard_mode="none"))
    insulated = kuffel_gradients(BASE.replace(guard_mode="insulated"))
    assert grounded[1] > none[1]
    assert insulated == none


def test_more_sub_conductors_lower_the_gradient_and_the_centre_phase_is_the_worst():
    df = sweep(BASE, "n_sub", [1, 2, 3, 4])
    col = df[SIM_B].to_numpy()
    assert np.all(np.diff(col) < 0)
    assert np.all(np.diff(df[KUF_B].to_numpy()) < 0)
    assert (df["peor fase (simulación)"] == "b").all()
    assert not df.loc[3, "cumple < 0,95 (simulación)"] and df.loc[4, "cumple < 0,95 (simulación)"]
    # the formula leads to the same design decision as the simulation for every case of the sweep
    assert (df["cumple < 0,95 (simulación)"] == df["cumple < 0,95 (Kuffel)"]).all()


def test_larger_phase_distance_lowers_the_gradient():
    make = lambda d: {"phases": ((-d, 25.0), (0.0, 25.0), (d, 25.0))}
    df = sweep(BASE, "distancia entre fases [m]", [8.0, 13.0], make=make)
    assert df[SIM_B].iloc[1] < df[SIM_B].iloc[0]


def test_single_conductor_per_phase_does_not_depend_on_the_spacing():
    a = analyze(BASE, n_sub=1, spacing=0.3)[SIM_B]
    b = analyze(BASE, n_sub=1, spacing=0.6)[SIM_B]
    assert a == pytest.approx(b, rel=1e-9)


def test_ground_field_is_symmetric_and_highest_outside_the_outer_phases():
    x, e = ground_profile(BASE)
    assert e == pytest.approx(e[::-1], rel=1e-6)
    assert abs(x[e.argmax()]) > 11.0                  # outside the outer phase (x = 11 m)
    assert e.max() > 3 * e[np.abs(x).argmin()]        # much higher than under the centre phase


def test_guard_modes():
    none = analyze(BASE, guard_mode="none")
    grounded = analyze(BASE)
    insulated = analyze(BASE, guard_mode="insulated")
    assert grounded[SIM_B] > none[SIM_B]
    assert insulated[SIM_B] == pytest.approx(none[SIM_B], rel=2e-3)
    wires, pot = build(BASE.replace(guard_mode="none"))
    assert len(wires) == 9 and set(pot) == {0, 1, 2}
    with pytest.raises(ValueError):
        build(BASE.replace(guard_mode="floating"))
