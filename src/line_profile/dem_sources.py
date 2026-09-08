# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Registry of the open, keyless DEM and canopy sources used by this module.

All URLs and tile-naming conventions were verified by reading a real point
((-74.5, 6.5)) on 2026-09-07. Copernicus and SRTM GL1 do NOT share one
tile-name format, even though both are 1-degree tiles: Copernicus tiles
carry an explicit arc-minute suffix ("N06_00_W075_00") while SRTM GL1 does
not ("N06W075"); the functions
below use the format actually measured, not the simplified one.
"""
import math
from dataclasses import dataclass
from typing import Callable


def _ns_ew_deg(lon: float, lat: float) -> tuple[str, str, int, int]:
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return ns, ew, abs(math.floor(lat)), abs(math.floor(lon))


def _tile_srtm_style(lon: float, lat: float) -> str:
    """1x1 deg tile, e.g. N06W075 (OpenTopography's SRTM_GL1 mirror)."""
    ns, ew, lat_deg, lon_deg = _ns_ew_deg(lon, lat)
    return f"{ns}{lat_deg:02d}{ew}{lon_deg:03d}"


def _tile_copernicus_style(lon: float, lat: float) -> str:
    """1x1 deg tile with explicit arc-minutes, e.g. N06_00_W075_00."""
    ns, ew, lat_deg, lon_deg = _ns_ew_deg(lon, lat)
    return f"{ns}{lat_deg:02d}_00_{ew}{lon_deg:03d}_00"


def _tile_glad_style(lon: float, lat: float) -> str:
    """10x10 deg tile named by its NW corner, e.g. 2020_10N_080W.

    Only verified for the northern hemisphere (Colombia); GLAD's southern
    convention is not exercised by this module's tests.
    """
    lat_band = (math.floor(lat / 10) + 1) * 10 if lat >= 0 else 0
    lon_band = math.floor(lon / 10) * -10
    return f"2020_{lat_band}N_{lon_band:03d}W"


@dataclass(frozen=True)
class DemSource:
    key: str
    url_template: str  # contains "{tile}"
    tile_deg: float
    tile_fn: Callable[[float, float], str]
    nodata_as_nan: bool


DEM_SOURCES: dict[str, DemSource] = {
    "cop30": DemSource(
        key="cop30",
        url_template=(
            "https://copernicus-dem-30m.s3.amazonaws.com/"
            "Copernicus_DSM_COG_10_{tile}_DEM/Copernicus_DSM_COG_10_{tile}_DEM.tif"
        ),
        tile_deg=1.0, tile_fn=_tile_copernicus_style, nodata_as_nan=True,
    ),
    "cop90": DemSource(
        key="cop90",
        url_template=(
            "https://copernicus-dem-90m.s3.amazonaws.com/"
            "Copernicus_DSM_COG_30_{tile}_DEM/Copernicus_DSM_COG_30_{tile}_DEM.tif"
        ),
        tile_deg=1.0, tile_fn=_tile_copernicus_style, nodata_as_nan=True,
    ),
    "srtm_gl1": DemSource(
        key="srtm_gl1",
        url_template="https://opentopography.s3.sdsc.edu/raster/SRTM_GL1/SRTM_GL1_srtm/{tile}.tif",
        tile_deg=1.0, tile_fn=_tile_srtm_style, nodata_as_nan=True,
    ),
    "glad_canopy": DemSource(
        key="glad_canopy",
        url_template="https://glad.umd.edu/users/Potapov/GLCLUC2020/Forest_height_2020/{tile}.tif",
        tile_deg=10.0, tile_fn=_tile_glad_style, nodata_as_nan=False,
    ),
}


def tile_name(dem_key: str, lon: float, lat: float) -> str:
    """URL of the mosaic that contains (lon, lat) for the given DEM source. R5, R6."""
    source = DEM_SOURCES[dem_key]
    return source.url_template.format(tile=source.tile_fn(lon, lat))


def _grid_starts(low: float, high: float, step: float) -> list[float]:
    """Multiples of step, from the one covering `low` up to (excluding) `high`."""
    start = math.floor(low / step) * step
    starts, v = [], start
    while v < high:
        starts.append(v)
        v += step
    return starts


def tiles_for_bbox(dem_key: str, bbox_4326: tuple[float, float, float, float]) -> list[str]:
    """All mosaic URLs that intersect bbox_4326 = (min_lon, min_lat, max_lon, max_lat).

    Covers a route that crosses the border between two mosaics: one URL per
    distinct tile, no duplicates.
    """
    source = DEM_SOURCES[dem_key]
    min_lon, min_lat, max_lon, max_lat = bbox_4326
    seen: set[str] = set()
    urls: list[str] = []
    for lat0 in _grid_starts(min_lat, max_lat, source.tile_deg):
        for lon0 in _grid_starts(min_lon, max_lon, source.tile_deg):
            tile = source.tile_fn(lon0 + source.tile_deg / 2, lat0 + source.tile_deg / 2)
            if tile not in seen:
                seen.add(tile)
                urls.append(source.url_template.format(tile=tile))
    return urls
