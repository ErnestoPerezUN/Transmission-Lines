# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Orchestrate the full profile: route -> sampled points -> DEM/canopy -> CSV.

Named elevation_profile.py, not profile.py: `profile` is a Python
standard-library module (the pure-Python profiler), and naming this file
that would shadow it for anything importing line_profile's modules -- a
real conflict found while writing this file's first test, not a style
preference.
"""
from pathlib import Path

import pandas as pd

from remote_dem import configure_remote_access, get_or_build_local_clip, sample_raster
from route import read_route, sample_route


def build_profile(route_path: str | Path, paso_m: float = 30.0,
                   dem_keys: tuple[str, ...] = ("cop30", "cop90", "srtm_gl1"),
                   reference_dem: str = "cop30",
                   cache_dir: str | Path = "data/line_profile") -> pd.DataFrame:
    """Build the full profile table for a route. R1-R9.

    reference_dem is the DEM used for suelo_estimado_m (R7) and pendiente_pct
    (R8); it must be one of dem_keys.
    """
    cache_dir = Path(cache_dir)
    configure_remote_access(cache_dir)

    info = read_route(route_path)
    df = sample_route(info.line, paso_m)

    bbox = (df["lon"].min(), df["lat"].min(), df["lon"].max(), df["lat"].max())
    points = list(zip(df["lon"], df["lat"]))

    for dem_key in dem_keys:
        clip_path = get_or_build_local_clip(dem_key, bbox, cache_dir)
        df[dem_key] = sample_raster(clip_path, points)

    canopy_clip = get_or_build_local_clip("glad_canopy", bbox, cache_dir)
    df["dosel_m"] = sample_raster(canopy_clip, points)

    df["suelo_estimado_m"] = df[reference_dem] - df["dosel_m"]

    df["pendiente_pct"] = 0.0
    dz = df[reference_dem].diff()
    dd = df["distancia_m"].diff()
    df.loc[df.index[1:], "pendiente_pct"] = (dz.iloc[1:] / dd.iloc[1:] * 100.0).to_numpy()

    return df


_CSV_HEADER_WARNING = (
    "# suelo_estimado_m es una ESTIMACION (elevacion del DEM menos dosel), no una medicion del "
    "terreno.\n"
    "# dosel_m proviene de GLAD Forest Height 2020: vegetacion no actualizada a la fecha del "
    "trazado.\n"
    "# Codigo academico, sin garantia; no usar para decisiones de diseno sin verificacion "
    "independiente. Ver DISCLAIMER.md.\n"
)


def write_profile_csv(df: pd.DataFrame, out_path: str | Path) -> None:
    """Write the profile as CSV: comma-decimal, no index. R9.

    Prefixed with '#' comment lines warning that suelo_estimado_m is an
    estimate and that the canopy is GLAD Forest Height 2020 (R7's risk).
    """
    out_path = Path(out_path)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        f.write(_CSV_HEADER_WARNING)
        df.to_csv(f, index=False, lineterminator="\n")
