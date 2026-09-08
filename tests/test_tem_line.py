# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Tests for the parallel-plate TEM line: FDTD ground truth vs circuit models.

These are physics assertions, not smoke tests. They pin down the four claims
the demo in ``examples/voltage_source_on_line.py`` makes:

  1. the FDTD wave travels at c and has Z0 = eta0 * h / w,
  2. an open end doubles the voltage,
  3. a single lumped L-C section is *wrong*, not merely inaccurate,
  4. a cascade of many sections converges onto the FDTD answer.
"""
import numpy as np
import pytest

from ladder_circuit import solve_ladder
from tem_line_fdtd import C0, ETA0, LineGeometry, TEMLineFDTD
from waveforms import step_source

# One shared FDTD run: it is the ground truth for most of the assertions below.
GEOM = LineGeometry()
T_RISE = GEOM.transit_time / 5.0
AMPLITUDE = 1.0
SOURCE = step_source(amplitude=AMPLITUDE, t_rise=T_RISE)
N_TRANSITS = 4.0


@pytest.fixture(scope="module")
def fdtd():
    sim = TEMLineFDTD(GEOM, SOURCE, snapshot_every=None)
    return sim.run(t_max=N_TRANSITS * GEOM.transit_time)


# --------------------------------------------------------------------------
# Geometry: the per-unit-length parameters must be exact, not approximate.
# --------------------------------------------------------------------------

def test_parallel_plate_parameters_are_exact():
    assert GEOM.velocity == pytest.approx(C0, rel=1e-12)
    assert GEOM.impedance == pytest.approx(ETA0 * GEOM.h / GEOM.width, rel=1e-12)
    assert 1.0 / np.sqrt(GEOM.L_per_m * GEOM.C_per_m) == pytest.approx(C0, rel=1e-12)


def test_operating_bandwidth_stays_below_parallel_plate_cutoff():
    """Above c/(2h) the guide is multi-mode and the TL model stops being valid."""
    bandwidth = 0.35 / T_RISE          # usual rise-time to -3 dB rule of thumb
    assert bandwidth < 0.5 * C0 / (2.0 * GEOM.h)


# --------------------------------------------------------------------------
# FDTD ground truth.
# --------------------------------------------------------------------------

def test_wave_travels_at_the_speed_of_light(fdtd):
    measured = fdtd.measured_velocity()
    assert measured == pytest.approx(C0, rel=0.02)


def test_characteristic_impedance_matches_eta0_h_over_w(fdtd):
    measured = fdtd.measured_impedance()
    assert measured == pytest.approx(GEOM.impedance, rel=0.05)


def test_wavefront_arrives_at_the_far_end_after_one_transit_time(fdtd):
    arrival = fdtd.arrival_time(fdtd.n_probes - 1, level=0.5 * AMPLITUDE)
    # A raised-cosine step crosses 50% half a rise time after the front lands.
    expected = GEOM.transit_time + 0.5 * T_RISE
    assert arrival == pytest.approx(expected, abs=0.05 * GEOM.transit_time)


def test_open_end_doubles_the_voltage(fdtd):
    far_end = fdtd.V_line[:, -1]
    assert far_end.max() == pytest.approx(2.0 * AMPLITUDE, rel=0.10)


def test_matched_load_kills_the_reflection(fdtd):
    """Z0 across the far end must flatten the echo the open end produces.

    A damped case is needed at all because an ideal source driving an open
    line is a lossless resonator: it rings at 4T for ever while a lumped L-C
    rings at 2*pi*T, so those two never settle into agreement however slowly
    the source is ramped.

    The match is good but not perfect, and deliberately so: the open-end
    fringing capacitance of test_field_derived_end_capacitance_repairs_the_far_end
    sits in parallel with the load, with time constant Z0*C_end = h/c = 0.1T.
    That is why the residual below is a few tenths rather than zero, and why
    it only bites for edges as fast as the transit time.
    """
    matched = TEMLineFDTD(GEOM, SOURCE, load=GEOM.impedance, snapshot_every=None).run(
        3.0 * GEOM.transit_time)
    mid = (matched.n_probes - 1) // 2
    late = matched.t > 1.6 * GEOM.transit_time

    settled = matched.V_line[late, mid]
    open_ended = fdtd.V_line[fdtd.t > 1.6 * GEOM.transit_time, mid]

    assert np.abs(settled - AMPLITUDE).max() < 0.3 * AMPLITUDE
    assert np.abs(open_ended - AMPLITUDE).max() > 0.8 * AMPLITUDE


def test_the_matched_value_is_the_best_value():
    """Pins down the conductivity of the terminating cells, not just its effect."""
    def residual(resistance):
        sim = TEMLineFDTD(GEOM, SOURCE, load=resistance, snapshot_every=None).run(3.0 * GEOM.transit_time)
        mid = (sim.n_probes - 1) // 2
        v = sim.V_line[sim.t > 1.6 * GEOM.transit_time, mid]
        return np.sqrt(np.mean((v - AMPLITUDE) ** 2))

    at_match = residual(GEOM.impedance)
    assert at_match < residual(0.5 * GEOM.impedance)
    assert at_match < residual(2.0 * GEOM.impedance)


def test_lumped_model_improves_as_the_edge_slows():
    """The validity criterion, as a number: the error vanishes with t_rise/T.

    Once the edge is slow compared with the transit time the lumped model
    converges quadratically.  Below that it saturates -- no amount of extra
    rise time in that regime helps, because the disagreement is a delay the
    model cannot represent at all.
    """
    T = GEOM.transit_time
    Z0 = GEOM.impedance
    errors = []
    for ratio in (1.0, 4.0, 16.0):
        source = step_source(AMPLITUDE, ratio * T)
        sim = TEMLineFDTD(GEOM, source, load=Z0,
                          snapshot_every=None).run((ratio + 4.0) * T)
        _, v = solve_ladder(1, GEOM.L_per_m, GEOM.C_per_m, GEOM.length,
                            source, sim.t, load=Z0)
        errors.append(np.sqrt(np.mean((v[-1] - sim.V_line[:, -1]) ** 2)))

    assert errors[0] > 0.10 * AMPLITUDE          # fast edge: plainly wrong
    assert errors[-1] < 0.01 * AMPLITUDE         # slow edge: plainly right
    assert all(b < a for a, b in zip(errors, errors[1:])), errors


# --------------------------------------------------------------------------
# Circuit models.
# --------------------------------------------------------------------------

def test_single_section_reproduces_the_analytic_lc_step_response():
    """With one section the ladder *is* a series L-C, so check it in closed form.

    Solving this with an accurate integrator matters: it makes the failure in
    ``test_lumped_model_disagrees_with_fdtd`` a modelling failure rather than
    an integration error a student could dismiss.
    """
    L = GEOM.L_per_m * GEOM.length
    C = GEOM.C_per_m * GEOM.length
    omega0 = 1.0 / np.sqrt(L * C)
    t = np.linspace(0.0, 4.0 / omega0, 400)

    _, v = solve_ladder(1, GEOM.L_per_m, GEOM.C_per_m, GEOM.length,
                        step_source(amplitude=AMPLITUDE, t_rise=0.0), t)

    analytic = AMPLITUDE * (1.0 - np.cos(omega0 * t))
    assert np.allclose(v[-1], analytic, atol=1e-4 * AMPLITUDE)


def test_lumped_model_disagrees_with_fdtd(fdtd):
    """The whole point of the demo: one section is qualitatively wrong."""
    assert _midline_rms_error(fdtd, n_sections=1) > 0.30 * AMPLITUDE


def test_many_sections_converge_onto_fdtd(fdtd):
    """Away from the ends, the ladder is not an approximation but a limit."""
    assert _midline_rms_error(fdtd, n_sections=200) < 1e-3 * AMPLITUDE


def test_error_decreases_as_sections_are_added(fdtd):
    errors = [_midline_rms_error(fdtd, n) for n in (1, 5, 20, 200)]
    assert all(b < a for a, b in zip(errors, errors[1:])), errors


def test_no_ladder_reproduces_the_open_end(fdtd):
    """Refining the ladder does *not* fix the far end, however large N gets.

    An open end fringes and radiates.  Both are properties of the geometry
    there, not of the line, so the far-end error floors out instead of
    converging -- unlike the mid-line error, which falls by three orders of
    magnitude over the same range of N.
    """
    coarse = _far_end_rms_error(fdtd, n_sections=5)
    fine = _far_end_rms_error(fdtd, n_sections=200)
    assert fine > 0.5 * coarse
    assert fine > 0.15 * AMPLITUDE


def test_field_derived_end_capacitance_repairs_the_far_end(fdtd):
    """The circuit model can be patched, but only with a number from the fields.

    The fringing capacitance is about C' * h, i.e. the open end behaves like
    roughly one conductor height of extra line.  There is no way to arrive at
    that from the ladder itself.
    """
    uncorrected = _far_end_rms_error(fdtd, n_sections=200)
    corrected = _far_end_rms_error(fdtd, n_sections=200,
                                   end_capacitance=GEOM.C_per_m * GEOM.h)
    assert corrected < 0.6 * uncorrected


def _ladder_node(n_sections, fraction):
    """Index of the ladder node nearest ``fraction`` of the way along the line.

    Node k sits at x = k*l/N, so N = 1 has a single node at the far end -- the
    lumped model genuinely has nowhere else to report.
    """
    return max(1, int(round(fraction * n_sections))) - 1


def _midline_rms_error(fdtd, n_sections):
    """RMS disagreement at mid-line, before the open-end reflection arrives.

    The window stops at 1.45T so that the comparison sees only the line
    itself; what the open end does is tested separately above.
    """
    window = fdtd.t <= 1.45 * GEOM.transit_time
    t = fdtd.t[window]
    _, v = solve_ladder(n_sections, GEOM.L_per_m, GEOM.C_per_m, GEOM.length,
                        SOURCE, t)
    probe = (fdtd.n_probes - 1) // 2
    node = _ladder_node(n_sections, 0.5)
    return float(np.sqrt(np.mean((v[node] - fdtd.V_line[window, probe]) ** 2)))


def _far_end_rms_error(fdtd, n_sections, end_capacitance=0.0):
    """RMS disagreement at the open end over the first two transits."""
    window = fdtd.t <= 2.0 * GEOM.transit_time
    t = fdtd.t[window]
    _, v = solve_ladder(n_sections, GEOM.L_per_m, GEOM.C_per_m, GEOM.length,
                        SOURCE, t, end_capacitance=end_capacitance)
    return float(np.sqrt(np.mean((v[-1] - fdtd.V_line[window, -1]) ** 2)))
