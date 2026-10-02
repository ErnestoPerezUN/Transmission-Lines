# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
# Refactored from the original EF_line.py script with the help of Claude (Anthropic); see NOTICE.md.
"""Electric field of an overhead line by the 2D charge simulation method.

Vocabulary used in this module and in the notebooks
---------------------------------------------------
* **Charge simulation** is the numerical method implemented here. Everything this module returns
  (fields, potentials, surface gradients) comes from it. It uses no closed formula for the gradient.
* **Fictitious line charges** are the unknowns of the method: a few auxiliary line charges placed inside
  each conductor. They are a mathematical tool, not real charges of the line. The real charge per unit
  length of a conductor is the sum of its fictitious charges (``charge_per_group``).
* **Formulas** (Peek, Kuffel / CIGRE, the equations of the CIGRE book) live in ``corona.py`` and are
  evaluated directly. Results are always labelled with the one that produced them.

What it computes
----------------
Given the geometry of a line (phases, sub-conductor bundles and ground wires) above a perfectly
conducting ground plane, it computes the electric field anywhere and, above all, the **surface
gradient of every conductor**, which is the number compared with Peek's critical gradient
(see ``corona.py``).

The method in four steps
------------------------
1. Each real conductor is replaced by ``n_charges`` **fictitious line charges** placed *inside* it, on a
   circle of radius ``charge_frac * r`` (default ``r/2``). Outside the conductor the field is that of
   these fictitious charges.
2. The ground is handled with the **image method**: a charge ``+q`` at ``(x, y)`` has an image
   ``-q`` at ``(x, -y)``, so the potential of the plane ``y = 0`` is zero.
3. Control points are chosen on the surface of every conductor and the potential there is forced to
   equal the conductor potential. This gives a linear system ``P q = V`` with the potential
   coefficient matrix ``P_ij = ln(r2/r1) / (2 pi eps0)``, where ``r1`` is the distance from
   control point ``i`` to charge ``j`` and ``r2`` its distance to the image of ``j``.
4. ``q`` is solved and the field is evaluated with the 2D Coulomb law
   ``E = q / (2 pi eps0) * r_vec / |r|^2``.

How ground wires are set to zero potential
------------------------------------------
A ground (shield) wire bonded to the tower at every structure is at **0 V**, but it is not a
conductor to ignore: it collects induced charge, and that charge changes the field of the phases.
The strategy is to keep it in the same linear system as one more conductor whose potential is zero
(``group_potentials[g] = 0``). It is not removed; its boundary condition simply has a zero
right-hand side.

For comparison, an **insulated** wire is available (``None`` in ``group_potentials``): its
potential becomes one more unknown and the condition "net charge = 0" is added. An insulated wire
floats close to the potential of the neighbouring phases and shields less than a grounded one.

Conventions
-----------
* Coordinates in metres; ``x`` horizontal, ``y`` height above the ground plane (``y > 0``).
* Potentials are **complex peak phasors** in volts. A line-to-line rms voltage ``V_LL`` gives a
  phase-to-ground peak of ``sqrt(2/3) * V_LL`` (see ``three_phase_potentials``).
* The problem is linear, so the field is also a complex phasor. On a conductor surface the field is
  normal to it and its peak value is the modulus of the normal component.
* Fields are returned in V/m. Multiply by ``V_PER_M_TO_KV_PER_CM`` to get kV/cm.

Limitations
-----------
2D model (infinite line, no sag, no towers), perfectly conducting ground, no corona space charge.
The error of the charge simulation (how far the surface potential is from the conductor potential,
relative to it) is reported in ``Solution.bc_error``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

EPS0 = 8.8541878128e-12          # vacuum permittivity [F/m]
K_COULOMB = 1.0 / (2.0 * np.pi * EPS0)
V_PER_M_TO_KV_PER_CM = 1.0e-5


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

@dataclass
class Wire:
    """A cylindrical conductor (or one sub-conductor of a bundle).

    x, y   : centre position [m]; y is the height above ground.
    r      : radio [m].
    group  : integer identifying the electrical conductor it belongs to. All the
             sub-conductors of one phase share the group (same potential).
    label  : text for plots (for example "A", "B", "C", "G1").
    """
    x: float
    y: float
    r: float
    group: int
    label: str = ""


def add_bundle(wires: list, x: float, y: float, n: int, spacing: float, r: float,
               group: int, label: str = "") -> list:
    """Append a bundle of ``n`` sub-conductors to ``wires`` and return the list.

    The sub-conductors sit at the vertices of a regular polygon of side ``spacing`` [m]
    (distance between adjacent sub-conductors) centred at ``(x, y)``. The bundle radius is
    ``R = spacing / (2 sin(pi/n))``. For ``n = 1`` a single conductor is added at ``(x, y)``.

    Orientation: n = 2 is horizontal, n = 4 is a square with horizontal and vertical sides, and
    n = 3 has one vertex pointing up.
    """
    if n == 1:
        wires.append(Wire(x, y, r, group, label))
        return wires
    R = spacing / (2.0 * np.sin(np.pi / n))
    angle0 = np.pi / 2 if n % 2 else np.pi / 2 + np.pi / n
    for k in range(n):
        a = angle0 + 2.0 * np.pi * k / n
        wires.append(Wire(x + R * np.cos(a), y + R * np.sin(a), r, group, label))
    return wires


def three_phase_potentials(v_ll_rms_kv: float, sequence: str = "abc") -> np.ndarray:
    """Phase-to-ground peak potentials [V] (complex) of a balanced three-phase system.

    ``v_ll_rms_kv`` is the line-to-line rms voltage in kV. Phase a is the real reference; b and c
    lag and lead by 120 degrees for the positive sequence ``"abc"`` (reversed for ``"acb"``).
    Returns a three-element array ordered [a, b, c].
    """
    vp = np.sqrt(2.0 / 3.0) * v_ll_rms_kv * 1e3
    a = np.exp(2j * np.pi / 3)
    if sequence == "abc":
        return vp * np.array([1.0, a ** 2, a])
    if sequence == "acb":
        return vp * np.array([1.0, a, a ** 2])
    raise ValueError("sequence must be 'abc' o 'acb'")


# ---------------------------------------------------------------------------
# Solution
# ---------------------------------------------------------------------------

@dataclass
class Solution:
    """Result of ``solve``. Fields are evaluated with the methods of this object."""
    wires: list
    charges: np.ndarray               # (N_wires, n_charges) fictitious line charges [C/m], complex
    charge_pos: np.ndarray            # (N_wires, n_charges, 2) position of each fictitious charge [m]
    group_potentials: dict            # potential of each group [V], floating ones included
    bc_error: float                   # max relative potential error on the surface
    n_charges: int = 0
    floating: tuple = field(default_factory=tuple)

    # -- potential and field at arbitrary points ------------------------------------------------
    def potential(self, X, Y):
        """Complex potential [V] at points (X, Y); valid outside the conductors."""
        X = np.asarray(X, float)
        Y = np.asarray(Y, float)
        v = np.zeros(X.shape, complex)
        for (qx, qy), q in zip(self.charge_pos.reshape(-1, 2), self.charges.reshape(-1)):
            r1 = np.hypot(X - qx, Y - qy)
            r2 = np.hypot(X - qx, Y + qy)
            v += K_COULOMB * q * np.log(r2 / r1)
        return v

    def field(self, X, Y):
        """Complex components (Ex, Ey) [V/m] at points (X, Y).

        Charge +q at (qx, qy) and its image -q at (qx, -qy)."""
        X = np.asarray(X, float)
        Y = np.asarray(Y, float)
        ex = np.zeros(X.shape, complex)
        ey = np.zeros(X.shape, complex)
        for (qx, qy), q in zip(self.charge_pos.reshape(-1, 2), self.charges.reshape(-1)):
            dx, dy = X - qx, Y - qy
            r1s = dx * dx + dy * dy
            dyi = Y + qy
            r2s = dx * dx + dyi * dyi
            ex += K_COULOMB * q * (dx / r1s - dx / r2s)
            ey += K_COULOMB * q * (dy / r1s - dyi / r2s)
        return ex, ey

    # -- instantaneous (time-domain) values --------------------------------------------------------
    # The phasors are peak values. The real value at the electrical angle wt is Re(phasor * e^{j wt}),
    # with wt = 2 pi f t: wt = 0 is the instant when phase a is at its positive peak.
    def potential_at(self, X, Y, wt=0.0):
        """Real potential [V] at the electrical angle ``wt`` [rad]."""
        return np.real(self.potential(X, Y) * np.exp(1j * wt))

    def field_at(self, X, Y, wt=0.0):
        """Real components (Ex, Ey) [V/m] at the electrical angle ``wt`` [rad]."""
        ex, ey = self.field(X, Y)
        ph = np.exp(1j * wt)
        return np.real(ex * ph), np.real(ey * ph)

    def surface_gradient_at(self, wt, n_points=360):
        """Signed normal field [V/m] on every wire surface at the angle ``wt``; shape (N_wires, n_points).

        Positive points out of the conductor (positive charge)."""
        phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
        out = np.zeros((len(self.wires), n_points))
        for i, w in enumerate(self.wires):
            ex, ey = self.field_at(w.x + w.r * np.cos(phi), w.y + w.r * np.sin(phi), wt)
            out[i] = ex * np.cos(phi) + ey * np.sin(phi)
        return out

    def max_gradient_vs_time(self, n_times=72, n_points=180):
        """Highest |normal field| [V/m] of each group along one cycle.

        Returns ``(wt, {group: array(n_times)})``. Each phase peaks at its own instant, 120 degrees
        apart, and the peak of each curve equals ``max_gradient_by_group``."""
        wts = np.linspace(0.0, 2.0 * np.pi, n_times, endpoint=False)
        res = {w.group: np.zeros(n_times) for w in self.wires}
        for k, wt in enumerate(wts):
            e = np.abs(self.surface_gradient_at(wt, n_points))
            for w, row in zip(self.wires, e):
                res[w.group][k] = max(res[w.group][k], row.max())
        return wts, res

    # -- surface gradient --------------------------------------------------------------
    def surface_gradient(self, n_points: int = 360):
        """Peak surface gradient of every wire.

        Returns ``(phi, E)``: ``phi`` are the angles [rad], shape (n_points,), and ``E`` is an
        (N_wires, n_points) array with the peak modulus of the normal component [V/m]."""
        phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
        out = np.zeros((len(self.wires), n_points))
        for i, w in enumerate(self.wires):
            px, py = w.x + w.r * np.cos(phi), w.y + w.r * np.sin(phi)
            ex, ey = self.field(px, py)
            en = ex * np.cos(phi) + ey * np.sin(phi)
            out[i] = np.abs(en)
        return phi, out

    def max_gradient_by_group(self, n_points: int = 360) -> dict:
        """Maximum peak surface gradient [V/m] of each group: ``{group: value}``."""
        _, E = self.surface_gradient(n_points)
        res = {}
        for w, e in zip(self.wires, E):
            res[w.group] = max(res.get(w.group, 0.0), float(e.max()))
        return res


def solve(wires: list, group_potentials: dict, n_charges: int = 12, charge_frac: float = 0.5) -> Solution:
    """Charge simulation: find the fictitious line charges that make every conductor surface equipotential.

    Each conductor is replaced by ``n_charges`` fictitious line charges inside it (see the module
    docstring). Their values are the solution of a linear system that imposes, at control points on
    the surface, the potential of the conductor. The field of the line then follows from them.

    Parameters
    ----------
    wires : list of ``Wire``.
    group_potentials : ``{group: potential}`` in volts (complex peak).
        * A number (``0`` included) fixes the potential of the group. **A grounded shield wire
          is declared with ``0``.**
        * ``None`` leaves the group *floating*: its potential is computed and its net charge
          is zero.
    n_charges : fictitious charges per conductor. More improve accuracy when conductors are close
        together (bundles); with 12 the boundary error is usually below 1 %.
    charge_frac : radius of the circle of fictitious charges as a fraction of the conductor radius.
        Must be < 1 (the charges sit inside the conductor).
    """
    if not 0.0 < charge_frac < 1.0:
        raise ValueError("charge_frac must be between 0 and 1")
    groups = sorted({w.group for w in wires})
    missing = [g for g in groups if g not in group_potentials]
    if missing:
        raise ValueError(f"missing potentials for groups {missing}")
    floating = [g for g in groups if group_potentials[g] is None]

    nw = len(wires)
    th = 2.0 * np.pi * np.arange(n_charges) / n_charges
    # Charges and control points share the angles th. Interleaving them (control points at
    # th + pi/n) makes P singular for even n_charges (the highest harmonic samples to zero).
    # Check points sit half-way between control points, where the boundary error is largest.
    cpos = np.zeros((nw, n_charges, 2))
    bpos = np.zeros((nw, n_charges, 2))
    epos = np.zeros((nw, n_charges, 2))       # check points (not used to solve)
    for i, w in enumerate(wires):
        cpos[i, :, 0] = w.x + charge_frac * w.r * np.cos(th)
        cpos[i, :, 1] = w.y + charge_frac * w.r * np.sin(th)
        bpos[i, :, 0] = w.x + w.r * np.cos(th)
        bpos[i, :, 1] = w.y + w.r * np.sin(th)
        epos[i, :, 0] = w.x + w.r * np.cos(th + np.pi / n_charges)
        epos[i, :, 1] = w.y + w.r * np.sin(th + np.pi / n_charges)

    def pmatrix(pts, chg):
        px, py = pts[:, 0][:, None], pts[:, 1][:, None]
        qx, qy = chg[:, 0][None, :], chg[:, 1][None, :]
        r1 = np.hypot(px - qx, py - qy)
        r2 = np.hypot(px - qx, py + qy)
        return K_COULOMB * np.log(r2 / r1)

    C = cpos.reshape(-1, 2)
    P = pmatrix(bpos.reshape(-1, 2), C)               # (N, N)
    Pe = pmatrix(epos.reshape(-1, 2), C)
    N = nw * n_charges
    wire_group = np.repeat([w.group for w in wires], n_charges)

    nf = len(floating)
    A = np.zeros((N + nf, N + nf), complex)
    b = np.zeros(N + nf, complex)
    A[:N, :N] = P
    for i in range(N):
        g = wire_group[i]
        if g in floating:
            A[i, N + floating.index(g)] = -1.0        # P q - V_f = 0
        else:
            b[i] = group_potentials[g]
    for k, g in enumerate(floating):                  # net charge zero
        A[N + k, :N] = (wire_group == g).astype(float)
    x = np.linalg.solve(A, b)
    q = x[:N]

    gp = dict(group_potentials)
    for k, g in enumerate(floating):
        gp[g] = x[N + k]

    v_target = np.array([gp[g] for g in wire_group])
    err = np.abs(Pe @ q - v_target)
    scale = max(np.max(np.abs(v_target)), 1e-30)
    return Solution(wires, q.reshape(nw, n_charges), cpos, gp, float(err.max() / scale),
                    n_charges, tuple(floating))


def charge_per_group(sol: Solution) -> dict:
    """Net line charge [C/m] (complex) of each group. A grounded shield wire carries induced charge."""
    res = {}
    for w, q in zip(sol.wires, sol.charges):
        res[w.group] = res.get(w.group, 0.0) + q.sum()
    return res


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

PALETTE = ["#C0392B", "#1F77B4", "#2CA02C", "#7F7F7F", "#9467BD", "#8C564B"]


def _colors(sol_or_wires, guard_groups=()):
    wires = sol_or_wires.wires if isinstance(sol_or_wires, Solution) else sol_or_wires
    groups = sorted({w.group for w in wires})
    col = {}
    k = 0
    for g in groups:
        if g in guard_groups:
            col[g] = "#555555"
        else:
            col[g] = PALETTE[k % len(PALETTE)]
            k += 1
    return col


def plot_layout(wires: list, guard_groups=(), zoom_group=None, figsize=(11, 5.5), title=None):
    """Draw the conductor layout above the ground plane.

    Left panel: full view with heights and ground. Right panel (if ``zoom_group`` is not None):
    the bundle of that group with its true radius.

    ``guard_groups`` draws those groups (grounded shield wires) in grey. Returns the figure.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    col = _colors(wires, guard_groups)
    ncol = 2 if zoom_group is not None else 1
    fig, axes = plt.subplots(1, ncol, figsize=figsize,
                             gridspec_kw={"width_ratios": [2, 1]} if ncol == 2 else None)
    ax = axes[0] if ncol == 2 else axes

    xs = [w.x for w in wires]
    ys = [w.y for w in wires]
    pad = 0.15 * (max(xs) - min(xs) + 4)
    ax.axhline(0, color="#8B7355", lw=3)
    ax.fill_between([min(xs) - pad, max(xs) + pad], -0.6, 0, color="#D9CBB0", alpha=0.6)
    # centre of each group, used for labels
    seen = {}
    for w in wires:
        seen.setdefault(w.group, []).append(w)
    for g, ws in seen.items():
        cx = np.mean([w.x for w in ws])
        cy = np.mean([w.y for w in ws])
        # drawn radius is exaggerated to be visible; the zoom panel shows the true radius
        for w in ws:
            ax.add_patch(Circle((w.x, w.y), max(w.r, 0.12), color=col[g], zorder=3))
        ax.annotate(ws[0].label or f"g{g}", (cx, cy), textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=11, fontweight="bold", color=col[g])
        ax.plot([cx, cx], [0, cy], ls=":", color=col[g], lw=1)
        ax.annotate(f"{cy:.1f} m", (cx, cy / 2), textcoords="offset points", xytext=(4, 0),
                    fontsize=9, color=col[g])
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(-0.6, max(ys) + 0.15 * max(ys) + 1.5)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("Height y [m]")
    ax.set_title(title or "Conductor layout (drawn radius exaggerated)")
    ax.grid(alpha=0.25)

    if ncol == 2:
        az = axes[1]
        ws = seen[zoom_group]
        cx = np.mean([w.x for w in ws])
        cy = np.mean([w.y for w in ws])
        rmax = max(np.hypot(w.x - cx, w.y - cy) for w in ws) + 2.5 * max(w.r for w in ws)
        for w in ws:
            az.add_patch(Circle((w.x, w.y), w.r, color=col[zoom_group], zorder=3))
        az.set_xlim(cx - rmax, cx + rmax)
        az.set_ylim(cy - rmax, cy + rmax)
        az.set_aspect("equal")
        az.set_title(f"Bundle {ws[0].label} to scale" + chr(10) + f"{len(ws)} sub-conductors, "
                     f"r = {ws[0].r * 100:.2f} cm", fontsize=10)
        az.set_xlabel("x [m]")
        az.grid(alpha=0.25)
    fig.tight_layout()
    return fig


def plot_surface_gradient(sol: Solution, guard_groups=(), e_crit_kv_cm=None, figsize=(8, 4.5)):
    """Peak surface gradient [kV/cm] versus angle for every energised conductor.

    If ``e_crit_kv_cm`` (peak critical gradient) is given, it is drawn as a horizontal line together
    with the 95 % line. Grounded shield wires are skipped (their gradient is low and does not
    govern the design)."""
    import matplotlib.pyplot as plt

    phi, E = sol.surface_gradient(360)
    col = _colors(sol, guard_groups)
    fig, ax = plt.subplots(figsize=figsize)
    for w, e in zip(sol.wires, E):
        if w.group in guard_groups:
            continue
        ax.plot(np.degrees(phi), e * V_PER_M_TO_KV_PER_CM, color=col[w.group], lw=1.2)
    if e_crit_kv_cm is not None:
        ax.axhline(e_crit_kv_cm, color="k", ls="--", label=f"$E_c$ = {e_crit_kv_cm:.1f} kV/cm")
        ax.axhline(0.95 * e_crit_kv_cm, color="k", ls=":", label="0.95 $E_c$")
        ax.legend()
    ax.set_xlabel("Angle on the conductor [deg] (0 = +x, 90 = up)")
    ax.set_ylabel("Peak surface gradient [kV/cm]")
    ax.set_title("Surface gradient of every sub-conductor")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def plot_field_map(sol: Solution, xlim, ylim, guard_groups=(), n=220, figsize=(9, 6), vmax_kv_m=None):
    """Map of the peak |E| [kV/m] around the line, with the conductors drawn.

    The modulus is sqrt(|Ex|^2 + |Ey|^2) of the phasors, an upper bound of the instantaneous
    (elliptically polarised) field. Points inside a conductor are masked."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    X, Y = np.meshgrid(np.linspace(*xlim, n), np.linspace(max(ylim[0], 1e-3), ylim[1], n))
    ex, ey = sol.field(X, Y)
    E = np.sqrt(np.abs(ex) ** 2 + np.abs(ey) ** 2) * 1e-3
    for w in sol.wires:
        E[np.hypot(X - w.x, Y - w.y) < w.r * 1.05] = np.nan
    col = _colors(sol, guard_groups)
    fig, ax = plt.subplots(figsize=figsize)
    cs = ax.contourf(X, Y, E, levels=np.linspace(0, vmax_kv_m or np.nanpercentile(E, 97), 16),
                     cmap="viridis", extend="max")
    fig.colorbar(cs, label="peak |E| [kV/m]")
    for w in sol.wires:
        ax.add_patch(Circle((w.x, w.y), max(w.r, 0.08), color=col[w.group], ec="w", zorder=3))
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Electric field modulus (conductors drawn with exaggerated radius)")
    fig.tight_layout()
    return fig


def _masked(sol, X, Y, *arrays, margin=1.02):
    """NaN-out the points that fall inside a conductor (the charges are only valid outside)."""
    inside = np.zeros(X.shape, bool)
    for w in sol.wires:
        inside |= np.hypot(X - w.x, Y - w.y) < w.r * margin
    return [np.where(inside, np.nan, a) for a in arrays]


def draw_bundle_zoom(ax, sol: Solution, group: int, wt: float = 0.0, half_width: float = None, n: int = 181,
                     e_max_kv_cm: float = None, n_equipotentials: int = 21, streamlines: bool = True):
    """Draw the instantaneous field around one bundle on ``ax``.

    Colour: |E| [kV/cm] at the electrical angle ``wt``. Black lines: equipotentials (spaced evenly
    between the smallest and largest potential in the window). Blue lines: field lines. The
    conductors are drawn with their true radius.

    ``half_width`` [m] is the half size of the window; the default is 2.2 times the bundle radius plus
    six sub-conductor radii. ``e_max_kv_cm`` fixes the colour scale so the frames of an animation
    are comparable.
    """
    from matplotlib.patches import Circle

    ws = [w for w in sol.wires if w.group == group]
    cx = float(np.mean([w.x for w in ws]))
    cy = float(np.mean([w.y for w in ws]))
    rb = max(np.hypot(w.x - cx, w.y - cy) for w in ws)
    hw = half_width or (2.2 * rb + 6.0 * ws[0].r)
    X, Y = np.meshgrid(np.linspace(cx - hw, cx + hw, n), np.linspace(cy - hw, cy + hw, n))
    ex, ey = sol.field_at(X, Y, wt)
    v = sol.potential_at(X, Y, wt) * 1e-3
    e = np.hypot(ex, ey) * V_PER_M_TO_KV_PER_CM
    e, v, ex, ey = _masked(sol, X, Y, e, v, ex, ey)
    vmax = e_max_kv_cm or float(np.nanpercentile(e, 99))
    # levels are denser at low |E| so the structure far from the conductors stays visible
    cs = ax.contourf(X, Y, e, levels=vmax * np.linspace(0, 1, 24) ** 1.8, cmap="YlOrRd", extend="max")
    if streamlines:
        ax.streamplot(X[0], Y[:, 0], np.nan_to_num(ex), np.nan_to_num(ey), color="#1F4E9C", linewidth=0.7,
                      density=0.9, arrowsize=0.8)
    ax.contour(X, Y, v, levels=np.linspace(np.nanmin(v), np.nanmax(v), n_equipotentials), colors="k", linewidths=0.7)
    for w in ws:
        ax.add_patch(Circle((w.x, w.y), w.r, fc="#BBBBBB", ec="k", zorder=5))
    ax.set_xlim(cx - hw, cx + hw)
    ax.set_ylim(cy - hw, cy + hw)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    return cs


def draw_field_map_at(ax, sol: Solution, wt: float, xlim, ylim, vmax_kv_m: float, n: int = 160):
    """Instantaneous |E| [kV/m] map of the whole line at the electrical angle ``wt`` (for animations)."""
    from matplotlib.patches import Circle

    X, Y = np.meshgrid(np.linspace(*xlim, n), np.linspace(max(ylim[0], 1e-3), ylim[1], n))
    ex, ey = sol.field_at(X, Y, wt)
    (e,) = _masked(sol, X, Y, np.hypot(ex, ey) * 1e-3)
    cs = ax.contourf(X, Y, e, levels=np.linspace(0, vmax_kv_m, 17), cmap="viridis", extend="max")
    for w in sol.wires:
        ax.add_patch(Circle((w.x, w.y), max(w.r, 0.08), fc="#DDDDDD", ec="k", zorder=5))
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    return cs
