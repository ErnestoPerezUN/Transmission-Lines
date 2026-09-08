# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Source waveforms shared by the FDTD and the lumped-circuit models.

Both models must be driven by the *identical* excitation, otherwise any
disagreement between them says more about the source than about the physics.
"""
import numpy as np


def raised_cosine_step(t, amplitude=1.0, t_rise=1e-9, t_delay=0.0):
    """A step of height ``amplitude`` smoothed over ``t_rise`` seconds.

    The transition is a raised cosine, so the waveform is C1-continuous: an
    ideal step would excite the grid all the way up to the Nyquist frequency
    and ring. ``t_rise <= 0`` gives a true Heaviside step, which is only useful
    for comparing against closed-form circuit solutions.
    """
    u = np.asarray(t, dtype=float) - t_delay
    if t_rise <= 0.0:
        ramp = np.where(u >= 0.0, 1.0, 0.0)
    else:
        # clip() handles both tails: u < 0 -> 0, u > t_rise -> 1.
        ramp = 0.5 * (1.0 - np.cos(np.pi * np.clip(u / t_rise, 0.0, 1.0)))
    out = amplitude * ramp
    return float(out) if out.ndim == 0 else out


def step_source(amplitude=1.0, t_rise=1e-9, t_delay=0.0):
    """Bind the parameters of :func:`raised_cosine_step` into a callable V(t)."""
    def source(t):
        return raised_cosine_step(t, amplitude, t_rise, t_delay)
    return source
