# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""2D FDTD (TE_z) model of a voltage source driving a horizontal conductor.

The conductor runs above a ground plane, so the pair forms a parallel-plate TEM
transmission line.  The structure is invariant along z; everything is computed
per unit width w in that direction.

    y
    ^   ............. Mur-1 absorbing boundary .............
    |
   h+   ==============================================   PEC conductor (Ex = 0)
    |   ||
    |  V(t)     TEM region: Ey downward, Hz out of page      [open end]
    |   ||
   0+   ==============================================   PEC ground  (Ex = 0)
    +------------------------------------------------------> x

Choosing a parallel plate rather than a round wire is what makes the comparison
against a circuit model honest: the per-unit-length parameters are *exact*,

    L' = mu0 h / w      C' = eps0 w / h
    v  = 1/sqrt(L'C') = c        Z0 = sqrt(L'/C') = eta0 h / w

with no thin-wire approximation and no numerical extraction step.  Any
disagreement with a circuit model is therefore a failure of the circuit model,
not an error in its parameters.

Fields are the TE_z set on a standard Yee grid:

    Ex[i, j]  at ((i+1/2) dx,  j      dy)    shape (Nx,   Ny+1)
    Ey[i, j]  at ( i      dx, (j+1/2) dy)    shape (Nx+1, Ny  )
    Hz[i, j]  at ((i+1/2) dx, (j+1/2) dy)    shape (Nx,   Ny  )

    dHz/dt = -(1/mu0)  (dEy/dx - dEx/dy)
    dEx/dt =  (1/eps0)  dHz/dy
    dEy/dt = -(1/eps0)  dHz/dx
"""
from dataclasses import dataclass

import numpy as np

C0 = 299792458.0
MU0 = 4.0e-7 * np.pi
EPS0 = 1.0 / (MU0 * C0 ** 2)
ETA0 = np.sqrt(MU0 / EPS0)


@dataclass(frozen=True)
class LineGeometry:
    """Physical layout and grid resolution.

    Defaults give a 1 m line 10 cm above ground: transit time 3.34 ns,
    Z0 = 37.7 ohm, and 20 cells across the gap.
    """
    h: float = 0.10             # conductor height above ground [m]
    length: float = 1.00        # conductor length [m]
    width: float = 1.00         # assumed width in z [m]
    dx: float = 5e-3            # cell size, both directions [m]
    left_margin: float = 0.15   # free space before the source [m]
    right_margin: float = 0.35  # free space past the open end [m]
    height_above: float = 0.30  # free space above the conductor [m]

    @property
    def dy(self):
        return self.dx

    @property
    def L_per_m(self):
        """Series inductance per unit length [H/m]."""
        return MU0 * self.h / self.width

    @property
    def C_per_m(self):
        """Shunt capacitance per unit length [F/m]."""
        return EPS0 * self.width / self.h

    @property
    def impedance(self):
        """Characteristic impedance sqrt(L'/C') [ohm]."""
        return np.sqrt(self.L_per_m / self.C_per_m)

    @property
    def velocity(self):
        """Phase velocity 1/sqrt(L'C') [m/s] -- exactly c for a parallel plate."""
        return 1.0 / np.sqrt(self.L_per_m * self.C_per_m)

    @property
    def transit_time(self):
        """One-way travel time along the line [s]."""
        return self.length / self.velocity

    @property
    def mode_cutoff(self):
        """Lowest non-TEM cutoff frequency [Hz]; above this the TL model dies."""
        return C0 / (2.0 * self.h)


class TEMLineFDTD:
    """Yee-grid TE_z solver for the geometry above.

    The excitation is an *ideal* voltage source: Ey is hard-forced across the
    gap at the left end of the conductor so that its line integral is exactly
    V(t).  That gives zero source impedance and a reflection coefficient of -1
    for returning waves -- the same thing an ideal source does in the circuit
    model, which is what keeps the two comparable.
    """

    def __init__(self, geometry=None, source=None, load=None, courant=0.99,
                 snapshot_every=8):
        """
        Parameters
        ----------
        load : float, optional
            Resistance placed across the far end [ohm], realised as a lossy
            column of cells spanning the gap.  Default is an open circuit.
            Setting it to the geometry's impedance matches the line, which
            removes the reflection and with it the resonance -- necessary if
            you want to see the regime where a lumped model is actually valid,
            since an undamped resonator has no quasi-static limit to agree in.
        snapshot_every : int or None
            Keep a copy of the whole Ey field every this many steps, for
            animation.  Pass None when only the probe traces are wanted: the
            snapshots dominate both memory and run time on long runs.
        """
        self.geom = geometry or LineGeometry()
        self.source = source
        self.load = load
        self.snapshot_every = snapshot_every

        g = self.geom
        self.dx = g.dx
        self.dy = g.dy
        self.Nx = int(round((g.left_margin + g.length + g.right_margin) / self.dx))
        self.Ny = int(round((g.h + g.height_above) / self.dy))
        self.jh = int(round(g.h / self.dy))                  # conductor row
        self.ia = int(round(g.left_margin / self.dx))        # source / near end
        self.ib = self.ia + int(round(g.length / self.dx))   # open / far end

        # 2D Courant limit.
        self.dt = courant / (C0 * np.sqrt(1.0 / self.dx ** 2 + 1.0 / self.dy ** 2))

        self.Ex = np.zeros((self.Nx, self.Ny + 1))
        self.Ey = np.zeros((self.Nx + 1, self.Ny))
        self.Hz = np.zeros((self.Nx, self.Ny))

        # Mur first-order coefficients.
        self._kx = (C0 * self.dt - self.dx) / (C0 * self.dt + self.dx)
        self._ky = (C0 * self.dt - self.dy) / (C0 * self.dt + self.dy)

        # Resistive termination: one column of lossy cells across the gap.  A
        # slab of conductivity sigma, thickness dx, height h and width w has
        # resistance h/(sigma dx w), so invert that for the sigma wanted.
        if load is None:
            self._ca_load = self._cb_load = None
        else:
            sigma = g.h / (load * self.dx * g.width)
            k = sigma * self.dt / (2.0 * EPS0)
            self._ca_load = (1.0 - k) / (1.0 + k)
            self._cb_load = (self.dt / EPS0) / (1.0 + k)

    # -- grid metadata -----------------------------------------------------

    @property
    def n_probes(self):
        """Number of voltage/current probes spread along the conductor."""
        return self.ib - self.ia + 1

    @property
    def probe_x(self):
        """Probe positions measured from the source [m]."""
        return np.arange(self.n_probes) * self.dx

    @property
    def extent(self):
        """(x0, x1, y0, y1) in metres, for imshow."""
        return (0.0, self.Nx * self.dx, 0.0, self.Ny * self.dy)

    @property
    def conductor_span(self):
        """(x_start, x_end) of the conductor in metres."""
        return (self.ia * self.dx, self.ib * self.dx)

    # -- time stepping -----------------------------------------------------

    def _step(self, t):
        """Advance the fields from time t to t + dt.

        Returns the (V, I) probe arrays sampled at time t.  Both have to be
        read mid-step: Ey is still at step n here, and Hz is time-centred on n
        by averaging the half-steps either side of it.  Sampling after the E
        update instead would offset V and I by half a time step and bias the
        extracted Z0 by roughly dt/(2 t_rise).
        """
        dx, dy, dt = self.dx, self.dy, self.dt
        jh, ia, ib = self.jh, self.ia, self.ib

        # Boundary values at step n, needed by the Mur update further down.
        ex_top_old = self.Ex[:, self.Ny].copy()
        ex_top_in_old = self.Ex[:, self.Ny - 1].copy()
        ey_left_old = self.Ey[0, :].copy()
        ey_left_in_old = self.Ey[1, :].copy()
        ey_right_old = self.Ey[self.Nx, :].copy()
        ey_right_in_old = self.Ey[self.Nx - 1, :].copy()

        hz_old = self.Hz.copy()
        ey_load_old = None if self.load is None else self.Ey[ib, :jh].copy()

        # H at n-1/2 -> n+1/2.
        self.Hz -= (dt / MU0) * (
            (self.Ey[1:, :] - self.Ey[:-1, :]) / dx
            - (self.Ex[:, 1:] - self.Ex[:, :-1]) / dy
        )

        hz_centred = 0.5 * (hz_old + self.Hz)          # time-centred on n
        V = -self.Ey[ia:ib + 1, :jh].sum(axis=1) * dy
        # Ampere loop around the top conductor.  The half-cell offset in x is
        # removed by averaging the Hz columns either side of each Ey column.
        I = -0.5 * (
            hz_centred[ia - 1:ib, :jh].mean(axis=1)
            + hz_centred[ia:ib + 1, :jh].mean(axis=1)
        ) * self.geom.width

        # E at n -> n+1.  The untouched rows/columns are boundaries: Ex[:, 0]
        # is the PEC ground and stays zero for the whole run.
        self.Ex[:, 1:-1] += (dt / EPS0) * (self.Hz[:, 1:] - self.Hz[:, :-1]) / dy
        self.Ey[1:-1, :] -= (dt / EPS0) * (self.Hz[1:, :] - self.Hz[:-1, :]) / dx

        # PEC conductor: tangential E vanishes on a zero-thickness sheet.
        self.Ex[ia:ib, jh] = 0.0

        # Redo the terminating column with the lossy update, which the plain
        # one above has just overwritten.
        if self.load is not None:
            self.Ey[ib, :jh] = (
                self._ca_load * ey_load_old
                - self._cb_load * (self.Hz[ib, :jh] - self.Hz[ib - 1, :jh]) / dx)

        # Mur-1 on the three open sides, so radiation leaves instead of
        # reflecting.  The bottom side is the physical ground plane.
        self.Ex[:, self.Ny] = ex_top_in_old + self._ky * (
            self.Ex[:, self.Ny - 1] - ex_top_old)
        self.Ey[0, :] = ey_left_in_old + self._kx * (self.Ey[1, :] - ey_left_old)
        self.Ey[self.Nx, :] = ey_right_in_old + self._kx * (
            self.Ey[self.Nx - 1, :] - ey_right_old)

        # Ideal voltage source across the gap: -integral(Ey dy) = V(t+dt).
        self.Ey[ia, :jh] = -self.source(t + dt) / self.geom.h

        return V, I

    def run(self, t_max):
        """Step to t_max and return self, carrying the recorded history."""
        n_steps = int(round(t_max / self.dt))
        self.t = np.arange(n_steps) * self.dt
        self.V_line = np.zeros((n_steps, self.n_probes))
        self.I_line = np.zeros((n_steps, self.n_probes))
        self.snapshots = []
        self.snapshot_times = []

        for n in range(n_steps):
            t = n * self.dt
            self.V_line[n], self.I_line[n] = self._step(t)
            if self.snapshot_every and n % self.snapshot_every == 0:
                self.snapshots.append(self.Ey.astype(np.float32))   # astype already copies
                self.snapshot_times.append(t)

        self.snapshot_times = np.asarray(self.snapshot_times)
        return self

    # -- measurements ------------------------------------------------------

    def arrival_time(self, probe, level):
        """Time at which a probe first rises through level, interpolated."""
        v = self.V_line[:, probe]
        above = np.flatnonzero(v >= level)
        if above.size == 0:
            raise ValueError(f"probe {probe} never reaches {level}")
        n = above[0]
        if n == 0:
            return self.t[0]
        v0, v1 = v[n - 1], v[n]
        return self.t[n - 1] + (level - v0) / (v1 - v0) * self.dt

    def measured_velocity(self, level=0.5):
        """Propagation velocity from the delay between two mid-line probes.

        Measured at the quarter and three-quarter points to stay clear of both
        the source discontinuity and the open end.
        """
        near = int(0.25 * (self.n_probes - 1))
        far = int(0.75 * (self.n_probes - 1))
        amplitude = self.source(self.t[-1])
        delay = (self.arrival_time(far, level * amplitude)
                 - self.arrival_time(near, level * amplitude))
        return (self.probe_x[far] - self.probe_x[near]) / delay

    def measured_impedance(self):
        """Z0 = V/I on the forward wave at mid-line.

        Sampled after the front has fully passed but before the reflection from
        the open end gets back, i.e. between T/2 and 3T/2.
        """
        mid = (self.n_probes - 1) // 2
        T = self.geom.transit_time
        window = (self.t > 0.75 * T) & (self.t < 1.4 * T)
        if not window.any():
            raise ValueError("run is too short to see a clean forward wave")
        return float(np.median(self.V_line[window, mid] / self.I_line[window, mid]))

    def summary(self):
        """Human-readable dump of the grid and the resulting line parameters."""
        g = self.geom
        return "\n".join([
            f"grid            {self.Nx} x {self.Ny} cells, dx = dy = {self.dx * 1e3:.1f} mm",
            f"time step       {self.dt * 1e12:.2f} ps "
            f"({g.transit_time / self.dt:.0f} steps per transit)",
            f"conductor       {g.length:.2f} m long, {g.h * 1e2:.0f} cm above ground",
            f"transit time T  {g.transit_time * 1e9:.3f} ns",
            f"L' (exact)      {g.L_per_m * 1e9:.1f} nH/m",
            f"C' (exact)      {g.C_per_m * 1e12:.1f} pF/m",
            f"Z0 (exact)      {g.impedance:.2f} ohm",
            f"non-TEM cutoff  {g.mode_cutoff / 1e9:.2f} GHz",
        ])
