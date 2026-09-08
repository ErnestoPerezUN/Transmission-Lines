# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Read and validate a transmission-line route from KML or GeoPackage.

R1-R4 of docs/specs/2026-09-07-perfil-topografico-linea-requirements.md.
"""
import math
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import pandas as pd
import shapely
from pyproj import CRS, Geod, Transformer
from shapely.geometry import LineString

#: Two segment endpoints closer than this (in EPSG:9377, metres) are treated
#: as the same point when concatenating tramos exported separately.
_JOIN_TOLERANCE_M = 1.0


@dataclass(frozen=True)
class RouteInfo:
    """A single 2D route, always in EPSG:4326.

    n_vertices: total vertex count of `line`.
    length_km: geodesic length on the WGS84 ellipsoid (not planar).
    """
    line: LineString
    n_vertices: int
    length_km: float


def read_route(path: str | Path) -> RouteInfo:
    """Read a route from a .kml or .gpkg file. R1.

    Reads all LineString/MultiLineString geometries in the file and
    concatenates them in order of connection (segment endpoints within
    _JOIN_TOLERANCE_M in EPSG:9377 are considered touching; R2), drops any Z
    coordinate (R4), and returns a single RouteInfo.
    """
    gdf = gpd.read_file(path)
    is_line = gdf.geom_type.isin(["LineString", "MultiLineString"])
    if not is_line.any():
        found = ", ".join(f"{kind} ({count})" for kind, count in gdf.geom_type.value_counts().items())
        raise ValueError(
            f"'{path}' no contiene ninguna geometría de línea (LineString/MultiLineString); "
            f"encontrado: {found or 'nada'}"
        )
    segments = [shapely.force_2d(geom) for geom in gdf.loc[is_line, "geometry"]]
    line = segments[0] if len(segments) == 1 else _concatenate_segments(segments)

    geod = Geod(ellps="WGS84")
    length_km = geod.geometry_length(line) / 1000.0

    return RouteInfo(line=line, n_vertices=len(line.coords), length_km=length_km)


def _concatenate_segments(segments: list[LineString], tol_m: float = _JOIN_TOLERANCE_M) -> LineString:
    """Chain segments end-to-end, in whatever order they connect.

    Raises ValueError naming the gap distance if the remaining segments never
    connect within tol_m of the growing chain's current end.
    """
    to_9377 = Transformer.from_crs(4326, 9377, always_xy=True)

    def endpoints(seg):
        start, end = seg.coords[0], seg.coords[-1]
        return to_9377.transform(*start), to_9377.transform(*end)

    remaining = list(segments)
    chain = list(remaining.pop(0).coords)
    while remaining:
        chain_end = to_9377.transform(*chain[-1])
        for i, seg in enumerate(remaining):
            start_pt, end_pt = endpoints(seg)
            if math.dist(chain_end, start_pt) <= tol_m:
                chain.extend(list(seg.coords)[1:])
                remaining.pop(i)
                break
            if math.dist(chain_end, end_pt) <= tol_m:
                chain.extend(list(seg.coords)[-2::-1])
                remaining.pop(i)
                break
        else:
            gaps = [min(math.dist(chain_end, p) for p in endpoints(seg)) for seg in remaining]
            nearest_gap = min(gaps)
            raise ValueError(
                f"la trayectoria tiene un hueco de {nearest_gap:.1f} m entre tramos "
                f"(tolerancia {tol_m} m); revisa que los segmentos conecten"
            )
    return LineString(chain)


def to_planar(line_4326: LineString, epsg: int = 9377) -> tuple[LineString, CRS]:
    """Reproject a route to a planar CRS (metres). R3."""
    transformer = Transformer.from_crs(4326, epsg, always_xy=True)
    coords = [transformer.transform(x, y) for x, y in line_4326.coords]
    return LineString(coords), CRS.from_epsg(epsg)


def sample_route(line_4326: LineString, paso_m: float, epsg: int = 9377) -> pd.DataFrame:
    """Sample a route every paso_m metres, measured on its planar projection. R3, R4.

    Columns: distancia_m [m, cumulative from the first point], lon, lat
    [degrees, EPSG:4326], x_9377, y_9377 [m], es_vertice [bool, True at the
    route's original vertices]. The first and last point are always included
    (as vertices); the original vertices are kept even where they fall off
    the regular paso_m grid, so a torre at a corredor's bend is never lost.
    """
    planar, _ = to_planar(line_4326, epsg)
    to_4326 = Transformer.from_crs(epsg, 4326, always_xy=True)
    length_m = planar.length

    coords = list(planar.coords)
    vertex_distances = [0.0]
    for a, b in zip(coords, coords[1:]):
        vertex_distances.append(vertex_distances[-1] + math.dist(a, b))

    grid = list(range(0, math.ceil(length_m), max(int(paso_m), 1))) if paso_m > 0 else []
    grid = [float(d) for d in grid if d < length_m]
    if not grid or not math.isclose(grid[-1], length_m):
        grid.append(length_m)

    vertex_set = {round(d, 6) for d in vertex_distances}
    all_distances = sorted({round(d, 6) for d in grid} | vertex_set)

    rows = []
    for d in all_distances:
        p = planar.interpolate(d)
        lon, lat = to_4326.transform(p.x, p.y)
        rows.append({
            "distancia_m": d, "lon": lon, "lat": lat, "x_9377": p.x, "y_9377": p.y,
            "es_vertice": d in vertex_set,
        })
    return pd.DataFrame(rows)
