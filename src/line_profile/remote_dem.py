# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Remote, keyless access to the DEM/canopy COGs in dem_sources.py, with a
local cache. R5, R6, R10, R11.
"""
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Sequence

import certifi
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.merge import merge as _rasterio_merge

from dem_sources import tiles_for_bbox


def _ascii_safe_path(path: Path) -> Path:
    """Return an ASCII-only path pointing at the same file, if one exists.

    On Windows, every long path has an 8.3 "short name" alias that is always
    ASCII (e.g. "Ernesto Pérez" -> "ERNEST~1"); GetShortPathNameW only works
    for a path that already exists on disk. Elsewhere (or if that call
    fails), a non-ASCII path has no such alias and the caller must not rely
    on one.
    """
    text = str(path)
    if text.isascii():
        return path
    if sys.platform == "win32":
        import ctypes
        buf = ctypes.create_unicode_buffer(260)
        if ctypes.windll.kernel32.GetShortPathNameW(text, buf, 260) and buf.value.isascii():
            return Path(buf.value)
    raise RuntimeError(
        f"'{path}' no tiene una ruta ASCII equivalente; certifi y GDAL fallan con rutas no ASCII"
    )


def configure_remote_access(cache_dir: str | Path) -> None:
    """Make /vsicurl/ reads work even when the cert bundle or cache_dir sit
    under a non-ASCII path. R11.

    Measured failure: on a machine whose certifi bundle lives under
    "C:\\Users\\Ernesto Pérez\\...", GDAL's cURL backend fails to open any
    remote COG with a UnicodeDecodeError. Copying the bundle into cache_dir
    and pointing CURL_CA_BUNDLE at its ASCII short-path alias fixes it --
    the copy alone is not enough if cache_dir itself sits under a non-ASCII
    path (e.g. the caller's own home directory). Also disables an extra
    directory listing GDAL otherwise issues against S3 buckets before every
    /vsicurl/ open.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cert_source = Path(certifi.where())
    dest = cache_dir / "cacert.pem"
    if str(cert_source).isascii() and str(dest).isascii():
        target = cert_source
    else:
        shutil.copy(cert_source, dest)
        target = _ascii_safe_path(dest)

    os.environ["CURL_CA_BUNDLE"] = str(target)
    os.environ["GDAL_DISABLE_READDIR_ON_OPEN"] = "EMPTY_DIR"


def sample_raster(local_clip_path: str | Path, points_4326: Sequence[tuple[float, float]]) -> np.ndarray:
    """Nearest-neighbour sample of a local (or /vsimem/) raster at each (lon, lat).

    Each value is a real reading of one DEM cell, never smoothed; nodata
    (e.g. -32768 in SRTM GL1) becomes NaN, never leaks out as a sentinel.

    points_4326 are always given in EPSG:4326 (every real DEM source in
    dem_sources.py already is). If the raster itself declares a different
    CRS -- as the synthetic plane used by the closing physical assertion
    does, in EPSG:9377, so its "30 m cell" is an exact metre, not an
    approximation at one latitude -- the points are reprojected to match it
    before sampling.
    """
    with rasterio.open(local_clip_path) as src:
        nodata = src.nodata
        if src.crs is not None and src.crs.to_epsg() != 4326:
            to_raster_crs = Transformer.from_crs(4326, src.crs, always_xy=True)
            points = [to_raster_crs.transform(lon, lat) for lon, lat in points_4326]
        else:
            points = points_4326
        values = np.array([v[0] for v in src.sample(points)], dtype="float64")
    if nodata is not None:
        values = np.where(values == nodata, np.nan, values)
    return values


def _open_remote(url: str):
    """Open a remote COG for reading. Isolated so tests can redirect it to a
    local file instead of a real network read."""
    return rasterio.open(f"/vsicurl/{url}")


def _cache_filename(dem_key: str, bbox_4326: tuple[float, float, float, float]) -> str:
    return f"{dem_key}_" + "_".join(f"{v:.4f}" for v in bbox_4326) + ".tif"


def _snap_bounds_to_pixel_grid(bounds, transform) -> tuple[float, float, float, float]:
    """Expand `bounds` outward to whole pixels of `transform`'s own grid.

    Measured bug this fixes: computing a crop's bounds as plain float
    subtraction (e.g. -74.51 - 0.02) can land a hair off the source
    raster's true grid line. For a query point that sits exactly on a pixel
    boundary (an unlucky but real case: (-74.5, 6.5), a whole 1-arcsecond
    tick, gave 199.15 m from a naively-cropped clip vs. the correct 198.12 m
    read straight from the source tile). Snapping in pixel space -- via the
    transform's own forward/inverse matrix, not degree arithmetic -- makes
    the crop's grid bit-identical to the source's, so no query point can
    fall on a different side of a boundary than it would on the source.
    """
    min_lon, min_lat, max_lon, max_lat = bounds
    inv = ~transform
    col_a, row_a = inv * (min_lon, max_lat)
    col_b, row_b = inv * (max_lon, min_lat)
    col0, col1 = math.floor(min(col_a, col_b)), math.ceil(max(col_a, col_b))
    row0, row1 = math.floor(min(row_a, row_b)), math.ceil(max(row_a, row_b))
    snapped_min_lon, snapped_max_lat = transform * (col0, row0)
    snapped_max_lon, snapped_min_lat = transform * (col1, row1)
    return (snapped_min_lon, snapped_min_lat, snapped_max_lon, snapped_max_lat)


def get_or_build_local_clip(dem_key: str, bbox_4326: tuple[float, float, float, float],
                             cache_dir: str | Path, margin_deg: float = 0.02) -> Path:
    """Return a local GeoTIFF covering bbox_4326 (+ margin_deg on every side). R5, R10.

    Built once from the remote tiles that intersect the (padded) bbox --
    merged with rasterio.merge if there is more than one, so a route that
    crosses a mosaic boundary shows no seam -- and reused on every later call
    with the same dem_key and bbox.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    min_lon, min_lat, max_lon, max_lat = bbox_4326
    padded = (min_lon - margin_deg, min_lat - margin_deg, max_lon + margin_deg, max_lat + margin_deg)
    cache_path = cache_dir / _cache_filename(dem_key, padded)
    if cache_path.exists():
        return cache_path

    urls = tiles_for_bbox(dem_key, padded)
    if not urls:
        raise ValueError(f"'{dem_key}': ningún mosaico cubre el bbox {padded}")

    datasets = []
    try:
        for url in urls:
            try:
                datasets.append(_open_remote(url))
            except Exception as exc:
                raise RuntimeError(
                    f"'{dem_key}': no se pudo leer '{url}' y no hay un recorte previo en "
                    f"'{cache_path}' (sin red y sin caché)"
                ) from exc
        snapped = _snap_bounds_to_pixel_grid(padded, datasets[0].transform)
        mosaic, transform = _rasterio_merge(datasets, bounds=snapped)
        profile = datasets[0].profile.copy()
        profile.update(driver="GTiff", height=mosaic.shape[1], width=mosaic.shape[2], transform=transform)
        with rasterio.open(cache_path, "w", **profile) as dst:
            dst.write(mosaic)
    finally:
        for ds in datasets:
            ds.close()
    return cache_path
