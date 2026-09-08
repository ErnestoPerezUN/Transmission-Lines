# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Connecting a voltage source to a horizontal conductor: fields vs circuits.

Run this to see, in four stages, why a transmission line is not a lumped circuit:

  Stage 1  The fields.  A step source is switched onto a 1 m conductor 10 cm
           above a ground plane.  Watch the wave leave the source, cross the
           line at c, double at the open end and come back.

  Stage 2  The lumped model fails.  Collapse the whole line into one L and one
           C -- the textbook lumped approximation -- and compare.  When the
           source rises slowly compared with the transit time the two agree.
           When it rises quickly they do not, and the lumped answer is not
           slightly off, it is a different waveform with no delay in it.  The
           third panel turns that into a criterion: the error falls off as
           (t_r/T)^-2 once the edge is slow, and saturates once it is not.
           This stage terminates the line in Z0, because an open line driven
           by an ideal source is a lossless resonator with no quiescent
           regime for the two models to agree in at all.

  Stage 3  Recovering the physics costs sections.  Chop the line into N L-C
           sections and sweep N.  Along the line the ladder converges onto the
           field solution beautifully -- but only in the limit, and the price
           of one sentence of field physics ("a wave travels at c") is a
           circuit that has to grow without bound.

  Stage 4  Where L' and C' came from.  The circuit model needs numbers that
           only a field calculation supplies -- and at the open end it needs
           one more, a fringing capacitance, that no amount of refinement can
           produce from within the ladder.

Usage
-----
    python examples/voltage_source_on_line.py              # interactive
    python examples/voltage_source_on_line.py --save out   # write PNGs + GIF
    python examples/voltage_source_on_line.py --no-anim    # skip stage 1
"""
import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "fdtd"))

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap

from ladder_circuit import lossless_line_step_response, solve_ladder
from tem_line_fdtd import C0, LineGeometry, TEMLineFDTD
from waveforms import step_source

# --------------------------------------------------------------------------
# Palette.  Validated against the data-viz colour checks: the ladder sweep is
# an ORDINAL ramp (N is an ordered quantity, so one hue light->dark), while the
# field solution is a separate categorical hue because it is a different kind
# of thing, not a bigger N.
# --------------------------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e6e5e1"

FDTD_C = "#eb6834"      # categorical slot 2 -- the field solution, throughout
LUMPED_C = "#2a78d6"    # categorical slot 1 -- the lumped circuit
LADDER_C = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # ordinal blue ramp
MIDLINE_C = "#2a78d6"
OPENEND_C = "#e34948"

# Diverging map for the signed field: two poles, neutral GRAY midpoint.
FIELD_CMAP = LinearSegmentedColormap.from_list(
    "field", ["#2a78d6", "#f0efec", "#e34948"])

N_SWEEP = (1, 5, 20, 200)
ERROR_STRIDE = 4        # subsample FDTD traces when computing RMS errors


def style():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": MUTED,
        "axes.titlecolor": INK,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "axes.grid": True,
        "axes.axisbelow": True,          # recessive grid, under the data
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "lines.linewidth": 2.0,          # thin marks
        "font.size": 9,
    })


# --------------------------------------------------------------------------
# Stage 1 -- the fields
# --------------------------------------------------------------------------

def stage1_animation(sim, geom, out_dir):
    """Animate Ey, with the line voltage V(x) tracking underneath it."""
    steps = np.round(sim.snapshot_times / sim.dt).astype(int)
    x0, x1 = sim.conductor_span
    # Span exactly the doubled field at the open end: the travelling TEM wave
    # then sits at half saturation and the reflection reaches full colour.
    scale = 2.0 * sim.source(sim.t[-1]) / geom.h

    fig, (ax, axv) = plt.subplots(
        2, 1, figsize=(9, 6), height_ratios=[2, 1], constrained_layout=True)

    im = ax.imshow(sim.snapshots[0].T, origin="lower", extent=sim.extent,
                   cmap=FIELD_CMAP, vmin=-scale, vmax=scale, aspect="equal",
                   interpolation="bilinear")
    ax.plot([x0, x1], [geom.h, geom.h], color=INK, lw=3, solid_capstyle="butt")
    ax.axhline(0.0, color=INK, lw=3)
    ax.plot([x0, x0], [0.0, geom.h], color=FDTD_C, lw=3)
    ax.annotate("V(t)", (x0, geom.h / 2), xytext=(-8, 0), textcoords="offset points",
                ha="right", va="center", color=FDTD_C, fontweight="bold")
    ax.annotate("open end", (x1, geom.h), xytext=(6, 8), textcoords="offset points",
                color=MUTED, fontsize=8)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Vertical electric field $E_y$")
    ax.grid(False)
    clock = ax.text(0.985, 0.93, "", transform=ax.transAxes, ha="right",
                    color=INK, fontsize=10, fontweight="bold")

    # Absolute x, on the same limits as the field map above, so the wavefront
    # lines up vertically between the two panels.
    (vline,) = axv.plot(sim.probe_x + x0, sim.V_line[0], color=FDTD_C)
    axv.set_xlim(*sim.extent[:2])
    axv.set_ylim(-0.4, 2.4)
    axv.set_xlabel("x [m]")
    axv.set_ylabel("V [V]")
    axv.set_title(r"Line voltage $V(x) = -\int_0^h E_y \, dy$")
    axv.axhline(1.0, color=MUTED, lw=1, ls=":")
    axv.annotate("source, 1 V", (sim.extent[1], 1.0), xytext=(-4, 4),
                 textcoords="offset points", ha="right", color=MUTED, fontsize=8)

    def update(k):
        im.set_array(sim.snapshots[k].T)
        vline.set_ydata(sim.V_line[steps[k]])
        clock.set_text(f"t = {sim.snapshot_times[k] / geom.transit_time:.2f} T")
        return im, vline, clock

    anim = FuncAnimation(fig, update, frames=len(sim.snapshots),
                         interval=40, blit=False)
    if out_dir:
        path = out_dir / "stage1_propagation.gif"
        anim.save(path, writer=PillowWriter(fps=20), dpi=80)
        print(f"  wrote {path}")
        plt.close(fig)
    return anim


# --------------------------------------------------------------------------
# Stage 2 -- the lumped model fails
# --------------------------------------------------------------------------

def stage2_lumped_fails(geom, out_dir):
    """Same circuit, different rise times: when is the lumped model allowed?

    The line is terminated in Z0 here, unlike everywhere else in this script.
    That matters: with an ideal source and an open end the structure is a
    lossless resonator, ringing at 4T while the lumped L-C rings at 2*pi*T, so
    the two can never settle into agreement however slowly the source is
    ramped.  Matching the line damps it, and only then does a regime exist in
    which the lumped model is genuinely correct.
    """
    T = geom.transit_time
    Z0 = geom.impedance
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.1), constrained_layout=True)

    for ax, (ratio, label) in zip(axes, [
        (8.0, "slow edge:  $t_r = 8T$"),
        (0.2, "fast edge:  $t_r = T/5$"),
    ]):
        source = step_source(1.0, ratio * T)
        span = (ratio + 4.0) * T
        sim = TEMLineFDTD(geom, source, load=Z0, snapshot_every=None).run(span)
        _, v = solve_ladder(1, geom.L_per_m, geom.C_per_m, geom.length,
                            source, sim.t, load=Z0)

        ax.plot(sim.t / T, sim.V_line[:, -1], color=FDTD_C, label="fields (FDTD)")
        ax.plot(sim.t / T, v[-1], color=LUMPED_C, ls="--",
                label="lumped $L$-$C$ (one section)")
        ax.set_title(label)
        ax.set_xlabel("time  [transit times $T$]")
        ax.set_ylabel("far-end voltage [V]")
        ax.set_xlim(0, span / T)
        ax.legend(loc="lower right")

    axes[0].annotate("agree to ~1% RMS; the small\novershoot is the lumped resonance",
                     xy=(9.0, 1.03), xytext=(0.3, 0.80), color=MUTED, fontsize=8,
                     arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
    axes[1].annotate("no delay: the lumped model\nresponds before the wave arrives",
                     xy=(0.85, 0.25), xytext=(0.04, 1.02), color=MUTED, fontsize=8,
                     arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))

    # The criterion itself, as a number.
    ax = axes[2]
    ratios = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
    errors = []
    for ratio in ratios:
        source = step_source(1.0, ratio * T)
        sim = TEMLineFDTD(geom, source, load=Z0,
                          snapshot_every=None).run((ratio + 4.0) * T)
        _, v = solve_ladder(1, geom.L_per_m, geom.C_per_m, geom.length,
                            source, sim.t, load=Z0)
        errors.append(np.sqrt(np.mean((v[-1] - sim.V_line[:, -1]) ** 2)))
    errors = np.array(errors)

    ax.loglog(ratios, errors, color=LUMPED_C, marker="o", ms=4)
    ax.axvline(1.0, color=MUTED, lw=1, ls=":")
    guide = errors[-3] * (ratios[-3] / ratios[ratios >= 4]) ** 2
    ax.loglog(ratios[ratios >= 4], guide, color=MUTED, lw=1, ls="--")
    ax.annotate("$\\propto (t_r/T)^{-2}$", (ratios[-1], guide[-1]),
                xytext=(-4, 10), textcoords="offset points", ha="right",
                color=MUTED, fontsize=8)
    ax.annotate("saturates: below $t_r \\approx T$ the\nerror is a missing delay,\n"
                "and no rise time fixes it",
                xy=(0.15, errors[1]), xytext=(0.12, 0.012), color=MUTED, fontsize=8,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
    ax.annotate("$t_r = T$", (1.0, errors.min()), xytext=(4, 0),
                textcoords="offset points", color=MUTED, fontsize=8)
    ax.set_title("The validity criterion, measured")
    ax.set_xlabel("rise time / transit time,  $t_r / T$")
    ax.set_ylabel("lumped-vs-fields RMS error [V]")

    fig.suptitle("Stage 2 - a lumped $L$-$C$ is right only when the edge is "
                 "slow compared with the transit time",
                 color=INK, fontweight="bold", fontsize=12)
    _finish(fig, out_dir, "stage2_lumped_fails.png")
    return ratios, errors


# --------------------------------------------------------------------------
# Stage 3 -- how many sections it takes
# --------------------------------------------------------------------------

def stage3_convergence(sim, geom, source, out_dir):
    T = geom.transit_time
    mid_probe = (sim.n_probes - 1) // 2
    window = sim.t <= 1.45 * T          # before the open-end echo returns
    t = sim.t[window]
    coarse = t[::ERROR_STRIDE]          # RMS needs far fewer points than the FDTD has

    fig, (ax, axe) = plt.subplots(1, 2, figsize=(11, 4.2),
                                  constrained_layout=True)

    ax.plot(sim.t / T, sim.V_line[:, mid_probe], color=FDTD_C, lw=2.6,
            label="fields (FDTD)", zorder=5)
    for n, colour in zip(N_SWEEP, LADDER_C):
        _, v = solve_ladder(n, geom.L_per_m, geom.C_per_m, geom.length, source, t)
        node = max(1, int(round(0.5 * n))) - 1
        ax.plot(t / T, v[node], color=colour, lw=1.6, label=f"ladder, N = {n}")
    ax.set_title("Mid-line voltage as sections are added")
    ax.set_xlabel("time  [transit times $T$]")
    ax.set_ylabel("voltage at $x = \\ell/2$ [V]")
    ax.set_xlim(0, 1.45)
    ax.legend(loc="upper left")

    # Error vs N, at mid-line and at the open end.
    n_grid = np.array([1, 2, 3, 5, 8, 12, 20, 35, 60, 100, 200])
    far_window = sim.t <= 2.0 * T
    t_far = sim.t[far_window][::ERROR_STRIDE]
    ref_mid = sim.V_line[window, mid_probe][::ERROR_STRIDE]
    ref_far = sim.V_line[far_window, -1][::ERROR_STRIDE]
    e_mid, e_far = [], []
    for n in n_grid:
        _, v = solve_ladder(n, geom.L_per_m, geom.C_per_m, geom.length,
                            source, coarse)
        node = max(1, int(round(0.5 * n))) - 1
        e_mid.append(np.sqrt(np.mean((v[node] - ref_mid) ** 2)))
        _, vf = solve_ladder(n, geom.L_per_m, geom.C_per_m, geom.length,
                             source, t_far)
        e_far.append(np.sqrt(np.mean((vf[-1] - ref_far) ** 2)))

    axe.loglog(n_grid, e_mid, color=MIDLINE_C, marker="o", ms=4,
               label="along the line ($x = \\ell/2$)")
    axe.loglog(n_grid, e_far, color=OPENEND_C, marker="s", ms=4,
               label="at the open end")
    axe.set_title("How the disagreement scales with N")
    axe.set_xlabel("number of $L$-$C$ sections, N")
    axe.set_ylabel("RMS error vs fields [V]")
    axe.legend(loc="lower left")
    # Labelled at each curve's own end rather than with arrows: the two
    # endpoints are three decades apart, so any arrow long enough to reach
    # crosses the whole panel.
    axe.annotate("converges: the ladder\nreally is the line",
                 xy=(n_grid[-1], e_mid[-1]), xytext=(-8, 12),
                 textcoords="offset points", ha="right",
                 color=MUTED, fontsize=8)
    axe.annotate("floors out: the open end\nis not a circuit element",
                 xy=(n_grid[-1], e_far[-1]), xytext=(-8, 14),
                 textcoords="offset points", ha="right",
                 color=MUTED, fontsize=8)

    fig.suptitle("Stage 3 - the circuit model has to grow without bound to say "
                 "what one field equation says",
                 color=INK, fontweight="bold", fontsize=12)
    _finish(fig, out_dir, "stage3_convergence.png")
    return n_grid, np.array(e_mid), np.array(e_far)


# --------------------------------------------------------------------------
# Stage 4 -- the parameters the circuit model cannot supply itself
# --------------------------------------------------------------------------

def fit_end_capacitance(sim, geom, source, n_sections=60):
    """Find the shunt capacitance at the far node that best matches the fields."""
    window = sim.t <= 2.5 * geom.transit_time
    t = sim.t[window][::ERROR_STRIDE]
    target = sim.V_line[window, -1][::ERROR_STRIDE]
    trials = np.linspace(0.0, 3.0 * geom.C_per_m * geom.h, 17)
    errors = []
    for c_end in trials:
        _, v = solve_ladder(n_sections, geom.L_per_m, geom.C_per_m, geom.length,
                            source, t, end_capacitance=c_end)
        errors.append(np.sqrt(np.mean((v[-1] - target) ** 2)))
    return trials[int(np.argmin(errors))]


def stage4_open_end(sim, geom, source, c_end, out_dir):
    T = geom.transit_time
    window = sim.t <= 2.5 * T
    t = sim.t[window]

    _, plain = solve_ladder(200, geom.L_per_m, geom.C_per_m, geom.length, source, t)
    _, fixed = solve_ladder(200, geom.L_per_m, geom.C_per_m, geom.length, source, t,
                            end_capacitance=c_end)

    fig, ax = plt.subplots(figsize=(8, 4.4), constrained_layout=True)
    ax.plot(t / T, lossless_line_step_response(t, T), color=MUTED, lw=1.2, ls=":",
            label="ideal lossless line (analytic)")
    ax.plot(t / T, sim.V_line[window, -1], color=FDTD_C, lw=2.6,
            label="fields (FDTD)", zorder=5)
    ax.plot(t / T, plain[-1], color=LADDER_C[1], lw=1.6,
            label="ladder, N = 200")
    ax.plot(t / T, fixed[-1], color=LADDER_C[3], lw=1.6, ls="--",
            label=f"ladder, N = 200 + $C_{{end}}$ = {c_end * 1e12:.1f} pF")
    ax.set_xlabel("time  [transit times $T$]")
    ax.set_ylabel("open-end voltage [V]")
    ax.set_title("Stage 4 - the open end needs a number the ladder cannot produce")
    # Upper left is the one empty corner here; lower right belongs to the note.
    ax.legend(loc="upper left")
    ax.annotate(
        f"$C_{{end}} / C' = {c_end / geom.C_per_m * 1e3:.0f}$ mm $\\approx h$:\n"
        "the end behaves like extra line,\nand only the fields say how much",
        xy=(1.55, 1.80), xytext=(2.45, 0.55), ha="right",
        color=MUTED, fontsize=8,
        arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
    _finish(fig, out_dir, "stage4_open_end.png")


# --------------------------------------------------------------------------

def _finish(fig, out_dir, name):
    if out_dir:
        path = out_dir / name
        fig.savefig(path, dpi=150)
        print(f"  wrote {path}")
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--save", metavar="DIR",
                        help="write figures to DIR instead of showing them")
    parser.add_argument("--no-anim", action="store_true",
                        help="skip the stage 1 animation")
    args = parser.parse_args()

    out_dir = None
    if args.save:
        out_dir = pathlib.Path(args.save)
        out_dir.mkdir(parents=True, exist_ok=True)
        mpl.use("Agg")
    style()

    geom = LineGeometry()
    T = geom.transit_time
    source = step_source(amplitude=1.0, t_rise=T / 5)

    sim = TEMLineFDTD(geom, source, snapshot_every=12).run(4.0 * T)

    print(sim.summary())
    print()
    print("measured from the simulated fields, not assumed:")
    print(f"  propagation velocity  {sim.measured_velocity() / C0:.4f} c")
    print(f"  characteristic Z0     {sim.measured_impedance():.2f} ohm "
          f"(exact: {geom.impedance:.2f})")

    c_end = fit_end_capacitance(sim, geom, source)
    print(f"  open-end fringing C   {c_end * 1e12:.1f} pF "
          f"= {c_end / geom.C_per_m * 1e3:.0f} mm of extra line "
          f"(h = {geom.h * 1e3:.0f} mm)")
    print()

    anim = None
    if not args.no_anim:
        anim = stage1_animation(sim, geom, out_dir)
    ratios, e_rise = stage2_lumped_fails(geom, out_dir)
    n_grid, e_mid, e_far = stage3_convergence(sim, geom, source, out_dir)
    stage4_open_end(sim, geom, source, c_end, out_dir)

    print("what the figures come to:")
    print(f"  lumped model, t_r = T/5      {e_rise[ratios == 0.2][0]:.3f} V RMS error")
    print(f"  lumped model, t_r = 16T      {e_rise[ratios == 16][0]:.4f} V RMS error")
    print(f"  sections for 1% along line   {n_grid[np.argmax(e_mid < 0.01)]}")
    print(f"  best any ladder does at the open end: {e_far.min():.3f} V RMS, "
          "and it does not improve with N")

    if not out_dir:
        plt.show()
    return anim


if __name__ == "__main__":
    main()
