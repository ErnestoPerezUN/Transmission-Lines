# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Tests for src/line_profile: reading a route, sampling it and reading DEMs.

Tests marked ``requiere_red`` hit real network services and are excluded
from the fast lane: ``python -m pytest tests -q -m "not requiere_red"``.
"""
import math
import os
import pathlib

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import LineString, Polygon

from types import SimpleNamespace

from dem_sources import DEM_SOURCES, tile_name, tiles_for_bbox
from elevation_profile import build_profile, write_profile_csv
from remote_dem import configure_remote_access, get_or_build_local_clip, sample_raster
from pyproj import Transformer

from route import read_route, sample_route, to_planar

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
REAL_KML = REPO_ROOT / "notebooks" / "line_profile" / "linea_la_virginia_medellin.kml"

# --------------------------------------------------------------------------
# T1: read_route() happy path (KML and GPKG).
# --------------------------------------------------------------------------

SYNTHETIC_LINE = LineString([(-75.6, 6.2), (-75.0, 6.4), (-74.5, 6.5)])


@pytest.fixture(params=["kml", "gpkg"])
def synthetic_route_file(tmp_path, request):
    ext, driver = request.param, {"kml": "KML", "gpkg": "GPKG"}[request.param]
    path = tmp_path / f"ruta.{ext}"
    gdf = gpd.GeoDataFrame({"name": ["tramo"]}, geometry=[SYNTHETIC_LINE], crs="EPSG:4326")
    gdf.to_file(path, driver=driver)
    return path


def test_read_route_kml_and_gpkg_roundtrip(synthetic_route_file):
    info = read_route(synthetic_route_file)
    assert info.n_vertices == 3
    assert info.line.has_z is False
    minx, miny, maxx, maxy = info.line.bounds
    assert -76.0 < minx < maxx < -74.0
    assert 6.0 < miny < maxy < 7.0
    assert info.length_km == pytest.approx(126.383, abs=0.5)


def test_read_route_real_kml():
    info = read_route(REAL_KML)
    assert info.n_vertices == 59
    assert info.length_km == pytest.approx(141.12, abs=0.1)


# --------------------------------------------------------------------------
# T2: reject a file that has no line geometry (R2). This is the exact defect
# of the old code, which derived a "centerline" from a substation polygon.
# --------------------------------------------------------------------------

def test_read_route_rejects_polygons_only(tmp_path):
    polygons = [
        Polygon([(0, 0), (0, 1), (1, 1), (1, 0)]),
        Polygon([(2, 2), (2, 3), (3, 3), (3, 2)]),
    ]
    path = tmp_path / "subestaciones.gpkg"
    gpd.GeoDataFrame({"name": ["a", "b"]}, geometry=polygons, crs="EPSG:4326").to_file(path, driver="GPKG")

    with pytest.raises(ValueError, match=r"Polygon.*2"):
        read_route(path)


# --------------------------------------------------------------------------
# T3: several disjoint segments (a route exported tramo by tramo). Tolerance
# is 1 m in EPSG:9377, per route.read_route()'s docstring.
# --------------------------------------------------------------------------

def _write_segments(tmp_path, segments):
    path = tmp_path / "tramos.gpkg"
    gpd.GeoDataFrame({"name": [f"t{i}" for i in range(len(segments))]},
                      geometry=segments, crs="EPSG:4326").to_file(path, driver="GPKG")
    return path


def test_read_route_concatenates_touching_segments(tmp_path):
    # segundo tramo empieza ~0.1 m del final del primero (bien dentro de la
    # tolerancia de 1 m), simulando dos exportaciones que casi coinciden.
    seg1 = LineString([(-75.6, 6.2), (-75.55, 6.25)])
    seg2 = LineString([(-75.55 + 1e-6, 6.25 + 1e-6), (-75.5, 6.30)])
    path = _write_segments(tmp_path, [seg1, seg2])

    info = read_route(path)

    # 3, no 4: el punto de arranque de seg2 queda a <1 m del final de seg1 y
    # se descarta como duplicado del empalme, en vez de dejar un segmento
    # fantasma de una fracción de metro entre dos puntos casi iguales.
    assert info.n_vertices == 3
    assert list(info.line.coords)[0] == pytest.approx((-75.6, 6.2))
    assert list(info.line.coords)[-1] == pytest.approx((-75.5, 6.30))


def test_read_route_rejects_disconnected_segments(tmp_path):
    seg1 = LineString([(-75.6, 6.2), (-75.55, 6.25)])
    seg2 = LineString([(-75.0, 6.9), (-74.9, 7.0)])  # a decenas de km, no a 1 m
    path = _write_segments(tmp_path, [seg1, seg2])

    with pytest.raises(ValueError, match=r"hueco"):
        read_route(path)


# --------------------------------------------------------------------------
# T4: a KML with altitude in its vertices (Google Earth sometimes writes it).
# That altitude is not an independent measurement (it is the same kind of
# source this module evaluates) and must never leak into length_km.
# --------------------------------------------------------------------------

def test_read_route_drops_z_from_kml(tmp_path):
    line_3d = LineString([(-75.6, 6.2, 1500), (-75.55, 6.25, 1800), (-75.5, 6.30, 1200)])
    path = tmp_path / "con_altitud.kml"
    gpd.GeoDataFrame({"name": ["tramo"]}, geometry=[line_3d], crs="EPSG:4326").to_file(path, driver="KML")

    info = read_route(path)

    assert info.line.has_z is False
    assert info.n_vertices == 3
    # sin la Z, la longitud geodésica es la de la proyección 2D, no la 3D
    # (que sería más larga por el desnivel de 1500->1800->1200 m).
    assert info.length_km < 15.7


# --------------------------------------------------------------------------
# T5: R3, distance accumulated over EPSG:9377 vs. WGS84 geodesic length.
# Fixtures fixed in requirements.md, R3.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("line, geodesic_km", [
    pytest.param(SYNTHETIC_LINE, 126.383, id="fixture-sintetica-3-vertices"),
])
def test_sample_route_planar_distance_matches_geodesic_within_tolerance(line, geodesic_km):
    df = sample_route(line, paso_m=1000.0)
    planar_km = df["distancia_m"].iloc[-1] / 1000.0
    assert planar_km == pytest.approx(geodesic_km, rel=0.001)  # tolerancia de R3: 0.1 %


def test_sample_route_planar_distance_matches_geodesic_for_real_route():
    info = read_route(REAL_KML)
    df = sample_route(info.line, paso_m=1000.0)
    planar_km = df["distancia_m"].iloc[-1] / 1000.0
    assert planar_km == pytest.approx(141.12, rel=0.001)


# --------------------------------------------------------------------------
# T6: configurable step + es_vertice flag (R4). No se fijan cifras geodésicas
# a mano: se cruzan contra to_planar(), ya verificado en T5.
# --------------------------------------------------------------------------

def test_sample_route_step_and_vertex_flag():
    paso_m = 50_000.0
    df = sample_route(SYNTHETIC_LINE, paso_m=paso_m)
    planar, _ = to_planar(SYNTHETIC_LINE)

    coords = list(planar.coords)
    vertex_distances = [0.0]
    for a, b in zip(coords, coords[1:]):
        vertex_distances.append(vertex_distances[-1] + math.dist(a, b))

    assert df.iloc[0]["distancia_m"] == 0.0 and bool(df.iloc[0]["es_vertice"])
    assert df.iloc[-1]["distancia_m"] == pytest.approx(vertex_distances[-1])
    assert bool(df.iloc[-1]["es_vertice"])

    marked = sorted(df.loc[df["es_vertice"], "distancia_m"])
    assert marked == pytest.approx(vertex_distances)

    gaps = df["distancia_m"].diff().dropna()
    assert (gaps <= paso_m + 1e-6).all(), "ningún salto debe superar el paso pedido"

    non_vertex = df.loc[~df["es_vertice"], "distancia_m"]
    assert all(
        math.isclose(d % paso_m, 0.0, abs_tol=1e-6) or math.isclose(d % paso_m, paso_m, abs_tol=1e-6)
        for d in non_vertex
    ), "los puntos que no son vértice deben caer en la grilla regular de paso_m"


# --------------------------------------------------------------------------
# T7: paso_m mayor que la longitud del tramo (caso borde de R4).
# --------------------------------------------------------------------------

def test_sample_route_step_larger_than_line_returns_two_points():
    short_line = LineString([(-75.6, 6.2), (-75.599, 6.201)])  # unas decenas de metros
    df = sample_route(short_line, paso_m=50_000.0)

    assert len(df) == 2
    assert df["es_vertice"].all()
    assert df["distancia_m"].iloc[0] == 0.0
    assert df["distancia_m"].iloc[-1] > 0.0


# --------------------------------------------------------------------------
# T8: dem_sources.tile_name() y DEM_SOURCES, convención de 1° (Copernicus y
# SRTM). Las URLs y el punto de referencia (-74.5, 6.5) son los medidos sin
# clave el 2026-09-07 (requirements.md, R5). Copernicus y SRTM no comparten
# literalmente el mismo formato de nombre de mosaico -- design.md lo daba por
# sentado ("N06W075 para SRTM/Copernicus"), pero el nombre real de Copernicus
# lleva minutos de arco explícitos ("N06_00_W075_00"); se corrige aquí con el
# formato medido, no con el simplificado.
# --------------------------------------------------------------------------

def test_tile_name_one_degree_convention():
    lon, lat = -74.5, 6.5

    cop30_url = tile_name("cop30", lon, lat)
    assert "copernicus-dem-30m.s3.amazonaws.com" in cop30_url
    assert "N06_00_W075_00" in cop30_url

    cop90_url = tile_name("cop90", lon, lat)
    assert "copernicus-dem-90m.s3.amazonaws.com" in cop90_url
    assert "N06_00_W075_00" in cop90_url

    srtm_url = tile_name("srtm_gl1", lon, lat)
    assert "opentopography.s3.sdsc.edu/raster/SRTM_GL1/" in srtm_url
    assert "N06W075" in srtm_url

    for key in ("cop30", "cop90", "srtm_gl1"):
        assert DEM_SOURCES[key].tile_deg == 1.0


# --------------------------------------------------------------------------
# T9: convención de 10° de GLAD Forest Height, medida el 2026-09-07.
# --------------------------------------------------------------------------

def test_tile_name_glad_ten_degree_convention():
    url = tile_name("glad_canopy", -74.5, 6.5)
    assert "glad.umd.edu" in url
    assert "2020_10N_080W" in url
    assert DEM_SOURCES["glad_canopy"].tile_deg == 10.0


# --------------------------------------------------------------------------
# T10: tiles_for_bbox() -- un bbox dentro de un solo mosaico, y uno que cruza
# la frontera de dos (caso borde: trayectoria que cruza el borde de mosaico).
# --------------------------------------------------------------------------

def test_tiles_for_bbox_single_and_crossing():
    single = tiles_for_bbox("srtm_gl1", (-74.9, 6.1, -74.6, 6.4))
    assert len(single) == 1

    crossing = tiles_for_bbox("srtm_gl1", (-74.9, 6.9, -74.6, 7.1))  # cruza N06/N07
    assert len(crossing) == 2
    assert len(set(crossing)) == 2


# --------------------------------------------------------------------------
# T11: remote_dem.configure_remote_access() -- rutas con tilde (R11). Caso
# real: en el equipo de Ernesto, certifi.where() vive bajo
# "C:\Users\Ernesto Pérez\...", y GDAL fallaba al abrir cualquier COG remoto
# hasta apuntar CURL_CA_BUNDLE a una copia en una ruta ASCII.
# --------------------------------------------------------------------------

def test_configure_remote_access_handles_non_ascii_cert_path(tmp_path, monkeypatch):
    accented_dir = tmp_path / "Ernesto Pérez"
    accented_dir.mkdir()
    fake_cert = accented_dir / "cacert.pem"
    fake_cert.write_text("-- CERTIFICADO DE PRUEBA --", encoding="utf-8")
    monkeypatch.setattr("remote_dem.certifi.where", lambda: str(fake_cert))

    cache_dir = tmp_path / "cache"
    saved_env = {k: os.environ.get(k) for k in ("CURL_CA_BUNDLE", "GDAL_DISABLE_READDIR_ON_OPEN")}
    try:
        configure_remote_access(cache_dir)

        copied = cache_dir / "cacert.pem"
        assert copied.exists()
        assert copied.read_text(encoding="utf-8") == fake_cert.read_text(encoding="utf-8")

        bundle = os.environ["CURL_CA_BUNDLE"]
        # el propio directorio de caché de la prueba hereda el nombre de
        # usuario del sistema (con tilde) de tmp_path: exactamente el caso
        # real que R11 tiene que sobrevivir, no solo el de certifi.where().
        assert bundle.isascii(), "CURL_CA_BUNDLE debe quedar en una ruta ASCII incluso si cache_dir no lo es"
        assert pathlib.Path(bundle).read_text(encoding="utf-8") == fake_cert.read_text(encoding="utf-8")
        assert os.environ["GDAL_DISABLE_READDIR_ON_OPEN"] == "EMPTY_DIR"
    finally:
        for key, value in saved_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _synthetic_geotiff(data, transform, nodata=None):
    """A GTiff written into an in-memory /vsimem/ path, given a 2D array."""
    memfile = rasterio.io.MemoryFile()
    with memfile.open(
        driver="GTiff", height=data.shape[0], width=data.shape[1], count=1,
        dtype=data.dtype, crs="EPSG:4326", transform=transform, nodata=nodata,
    ) as dst:
        dst.write(data, 1)
    return memfile  # el llamador debe mantenerlo vivo (memfile.name) y cerrarlo


# --------------------------------------------------------------------------
# T12: remote_dem.sample_raster() -- vecino más cercano, caso feliz.
# --------------------------------------------------------------------------

def test_sample_raster_nearest_neighbor():
    data = np.array([[10.0, 20.0], [30.0, 40.0]], dtype="float32")
    transform = from_origin(0, 2, 1, 1)  # celdas de 1x1 grado, origen (0, 2)
    with _synthetic_geotiff(data, transform) as memfile:
        values = sample_raster(memfile.name, [(0.5, 1.5), (1.5, 1.5), (0.5, 0.5), (1.5, 0.5)])

    assert list(values) == pytest.approx([10.0, 20.0, 30.0, 40.0])


# --------------------------------------------------------------------------
# T13: remote_dem.sample_raster() -- nodata se convierte en NaN.
# --------------------------------------------------------------------------

def test_sample_raster_nodata_becomes_nan():
    data = np.array([[10.0, -32768.0]], dtype="float32")
    transform = from_origin(0, 1, 1, 1)
    with _synthetic_geotiff(data, transform, nodata=-32768.0) as memfile:
        values = sample_raster(memfile.name, [(0.5, 0.5), (1.5, 0.5)])

    assert values[0] == pytest.approx(10.0)
    assert math.isnan(values[1])


def _write_synthetic_plane_tile(path, lon0, lat0, cell=0.1, n=10, mx=100.0, my=50.0, z0=0.0):
    """A GeoTIFF covering [lon0, lon0+n*cell] x [lat0, lat0+n*cell] with
    z(lon, lat) = z0 + mx*lon + my*lat -- a synthetic plane, not physical
    units. Used to check merges/caching preserve the field, not its scale.
    """
    arr = np.zeros((n, n), dtype="float32")
    transform = from_origin(lon0, lat0 + n * cell, cell, cell)
    for row in range(n):
        for col in range(n):
            lon = lon0 + (col + 0.5) * cell
            lat = (lat0 + n * cell) - (row + 0.5) * cell
            arr[row, col] = z0 + mx * lon + my * lat
    with rasterio.open(
        path, "w", driver="GTiff", height=n, width=n, count=1, dtype="float32",
        crs="EPSG:4326", transform=transform,
    ) as dst:
        dst.write(arr, 1)


# --------------------------------------------------------------------------
# T14: get_or_build_local_clip() combina mosaicos sin costura. La apertura
# remota (remote_dem._open_remote) y tiles_for_bbox se parchan para leer de
# disco local en vez de red -- design.md llama a esto "la ruta remota
# parcheada".
# --------------------------------------------------------------------------

def test_get_or_build_local_clip_merges_without_seam(tmp_path, monkeypatch):
    cell = 0.1
    tile_w, tile_e = tmp_path / "w.tif", tmp_path / "e.tif"
    _write_synthetic_plane_tile(tile_w, lon0=-75.0, lat0=6.0, cell=cell)
    _write_synthetic_plane_tile(tile_e, lon0=-74.0, lat0=6.0, cell=cell)

    monkeypatch.setattr("remote_dem.tiles_for_bbox", lambda dem_key, bbox: [str(tile_w), str(tile_e)])
    monkeypatch.setattr("remote_dem._open_remote", lambda url: rasterio.open(url))

    clip_path = get_or_build_local_clip("srtm_gl1", (-75.0, 6.0, -73.0, 7.0), tmp_path / "cache", margin_deg=0.0)

    points = [(-74.001, 6.5), (-73.999, 6.5)]  # a 2 mm de la costura, uno en cada mosaico
    values = sample_raster(clip_path, points)
    analytic = [100.0 * lon + 50.0 * lat for lon, lat in points]
    # cota del vecino más cercano: en x y en y el punto puede estar hasta
    # 0.5*celda del centro de su celda, cada eje por su lado (no combinados
    # por Pitágoras -- el desplazamiento máximo es una caja, no un círculo).
    bound = 0.5 * cell * (abs(100.0) + abs(50.0)) + 1e-3

    # cada lado de la costura, muestreado de un mosaico distinto, cae dentro
    # de la misma cota de vecino más cercano que cualquier punto interior de
    # un solo mosaico -- si el merge desalineara los mosaicos, insertara
    # nodata en la unión o leyera el mosaico equivocado, el error superaría
    # esta cota (una diferencia de una celda completa entre los dos lados no
    # prueba nada por sí sola: eso pasa también dentro de un único mosaico).
    assert not np.isnan(values).any()
    for measured, expected in zip(values, analytic):
        assert abs(measured - expected) <= bound


# --------------------------------------------------------------------------
# T15: sin red y sin caché -> error explícito (no una traza cruda).
# --------------------------------------------------------------------------

def test_get_or_build_local_clip_no_network_no_cache_raises(tmp_path, monkeypatch):
    monkeypatch.setattr("remote_dem.tiles_for_bbox", lambda dem_key, bbox: ["https://fake/tile.tif"])

    def _boom(url):
        raise OSError("simulated network failure")

    monkeypatch.setattr("remote_dem._open_remote", _boom)

    with pytest.raises(RuntimeError, match=r"srtm_gl1.*https://fake/tile\.tif"):
        get_or_build_local_clip("srtm_gl1", (-75.0, 6.0, -74.9, 6.1), tmp_path / "cache")


# --------------------------------------------------------------------------
# T16: caché local reutilizable (R10).
# --------------------------------------------------------------------------

def test_get_or_build_local_clip_reuses_cache(tmp_path, monkeypatch):
    tile = tmp_path / "single.tif"
    _write_synthetic_plane_tile(tile, lon0=-75.0, lat0=6.0)

    calls = {"n": 0}

    def counting_open(url):
        calls["n"] += 1
        return rasterio.open(url)

    monkeypatch.setattr("remote_dem.tiles_for_bbox", lambda dem_key, bbox: [str(tile)])
    monkeypatch.setattr("remote_dem._open_remote", counting_open)

    cache_dir = tmp_path / "cache"
    bbox = (-75.0, 6.0, -74.9, 6.1)
    first = get_or_build_local_clip("srtm_gl1", bbox, cache_dir, margin_deg=0.0)
    assert calls["n"] == 1

    second = get_or_build_local_clip("srtm_gl1", bbox, cache_dir, margin_deg=0.0)
    assert calls["n"] == 1, "la segunda llamada no debe volver a abrir la fuente remota"
    assert first == second

    point = [(-74.95, 6.05)]
    assert sample_raster(first, point) == pytest.approx(sample_raster(second, point))


# --------------------------------------------------------------------------
# T17: medición con red del tamaño del recorte de dosel GLAD para el
# sub-tramo real km 70-100 (design.md, "Riesgos y mitigaciones"). El tile
# completo pesa 487 MB (medido en el brainstorming); 50 MB es un límite de
# sanidad generoso, no la cifra que importa reportar -- esa es la medida.
# --------------------------------------------------------------------------

@pytest.mark.requiere_red
def test_glad_clip_size_is_reasonable_for_real_segment(tmp_path):
    cache_dir = tmp_path / "cache"
    configure_remote_access(cache_dir)  # R11: sin esto, /vsicurl/ falla en equipos con tilde en la ruta

    info = read_route(REAL_KML)
    df = sample_route(info.line, paso_m=500.0)
    sub = df[(df["distancia_m"] >= 70_000) & (df["distancia_m"] <= 100_000)]
    bbox = (sub["lon"].min(), sub["lat"].min(), sub["lon"].max(), sub["lat"].max())

    clip_path = get_or_build_local_clip("glad_canopy", bbox, cache_dir)
    size_mb = clip_path.stat().st_size / (1024 * 1024)
    print(f"\n[T17] recorte de dosel GLAD para km 70-100: {size_mb:.3f} MB (bbox={tuple(round(v, 4) for v in bbox)})")

    assert size_mb < 50.0


# --------------------------------------------------------------------------
# T19: elevation_profile.build_profile() -- suelo_estimado_m = elevación de
# reference_dem menos dosel_m (R7). read_route/sample_route/get_or_build_
# local_clip/sample_raster se parchan: esta tarea prueba solo la columna
# derivada, no la lectura de ruta ni el acceso remoto (ya probados).
# --------------------------------------------------------------------------

def test_build_profile_suelo_estimado_is_elevation_minus_canopy(monkeypatch, tmp_path):
    fake_route = SimpleNamespace(line=None, n_vertices=3, length_km=1.0)
    fake_sampled = pd.DataFrame({
        "distancia_m": [0.0, 30.0, 60.0],
        "lon": [-75.0, -75.0, -75.0], "lat": [6.0, 6.0, 6.0],
        "x_9377": [0.0, 30.0, 60.0], "y_9377": [0.0, 0.0, 0.0],
        "es_vertice": [True, False, True],
    })
    elevations = {"cop30": [100.0, 110.0, 120.0], "cop90": [95.0, 108.0, 118.0], "srtm_gl1": [98.0, 109.0, 119.0]}
    canopy = [10.0, 5.0, 0.0]

    monkeypatch.setattr("elevation_profile.read_route", lambda path: fake_route)
    monkeypatch.setattr("elevation_profile.sample_route", lambda line, paso_m: fake_sampled.copy())
    monkeypatch.setattr("elevation_profile.configure_remote_access", lambda cache_dir: None)
    monkeypatch.setattr("elevation_profile.get_or_build_local_clip", lambda dem_key, bbox, cache_dir: dem_key)

    def fake_sample_raster(clip_token, points):
        if clip_token == "glad_canopy":
            return np.array(canopy)
        return np.array(elevations[clip_token])

    monkeypatch.setattr("elevation_profile.sample_raster", fake_sample_raster)

    df = build_profile("ruta-falsa.kml", paso_m=30.0, cache_dir=tmp_path)

    assert list(df["suelo_estimado_m"]) == pytest.approx(
        [c - d for c, d in zip(elevations["cop30"], canopy)]
    )


# --------------------------------------------------------------------------
# T20: elevation_profile.write_profile_csv() -- columnas, formato y avisos
# del encabezado (R9).
# --------------------------------------------------------------------------

def test_write_profile_csv_format_and_header_warnings(tmp_path):
    df = pd.DataFrame({
        "distancia_m": [0.0, 30.0],
        "lon": [-75.0, -75.0001], "lat": [6.0, 6.0002],
        "x_9377": [0.0, 30.0], "y_9377": [0.0, 5.0],
        "es_vertice": [True, False],
        "cop30": [100.0, 105.0], "cop90": [98.0, 103.0], "srtm_gl1": [99.0, 104.0],
        "dosel_m": [10.0, 8.0], "suelo_estimado_m": [90.0, 97.0],
        "pendiente_pct": [0.0, 16.67],
    })
    out_path = tmp_path / "perfil.csv"

    write_profile_csv(df, out_path)

    text = out_path.read_text(encoding="utf-8")
    lines = text.splitlines()
    comment_lines = [l for l in lines if l.startswith("#")]
    comment_text = "\n".join(comment_lines)
    assert comment_lines, "el CSV debe empezar con al menos una línea de comentario '#'"
    assert "suelo_estimado_m" in comment_text and "estimaci" in comment_text.lower()
    assert "GLAD" in comment_text and "2020" in comment_text

    header_line = next(l for l in lines if not l.startswith("#"))
    columns = header_line.split(",")
    required = {
        "distancia_m", "lon", "lat", "x_9377", "y_9377", "es_vertice",
        "cop30", "cop90", "srtm_gl1", "dosel_m", "suelo_estimado_m", "pendiente_pct",
    }
    assert required.issubset(set(columns))
    assert "." in [l for l in lines if l.startswith("0.0,")][0]  # punto decimal, no coma

    read_back = pd.read_csv(out_path, comment="#")
    assert len(read_back) == 2
    assert "Unnamed: 0" not in read_back.columns  # sin índice escrito


# --------------------------------------------------------------------------
# T21: aserción física de cierre -- plano inclinado sintético en EPSG:9377,
# celda de 30 m, línea recta de 500 m muestreada con paso_m = 10 m mediante
# el pipeline real (route.sample_route + elevation_profile.build_profile);
# solo la lectura del archivo de ruta y el acceso remoto se parchan.
#
# Verificado aquí con números reales antes de escribir la aserción: con
# paso_m=10 m y celda=30 m, la mayoría de los pasos caen dentro de
# la MISMA celda (pendiente_pct = 0 entre ellos) y el paso que cruza una
# celda muestra hasta 7 puntos porcentuales de desvío frente a la pendiente
# analítica de 5 %; una tolerancia de 0.1 pp por paso individual es
# matemáticamente imposible con paso_m menor que la celda -- es la propia
# lección del módulo (muestrear más fino que la celda no agrega información,
# y aquí además introduce ruido). Lo que sí se sostiene, y es lo que se
# verifica, es la pendiente PROMEDIO de extremo a extremo: su cota de error
# es la del vecino más cercano en los dos puntos extremos únicamente
# (cota = celda*(|mx|+|my|)/longitud*100 = 0.42 pp aquí; medido: 0.08 pp).
# --------------------------------------------------------------------------

def test_closing_physical_assertion_synthetic_inclined_plane(tmp_path, monkeypatch):
    cell = 30.0
    z0, mx, my = 1000.0, 0.03, 0.04  # pendiente analítica s = 5 % (triángulo 3-4-5)
    s = math.hypot(mx, my) * 100.0

    to_9377 = Transformer.from_crs(4326, 9377, always_xy=True)
    to_4326_ = Transformer.from_crs(9377, 4326, always_xy=True)
    start_lon, start_lat = -74.5, 6.5  # mismo punto de referencia que R5/R6
    sx, sy = to_9377.transform(start_lon, start_lat)
    ex, ey = sx + 300.0, sy + 400.0  # 500 m exactos, dirección 3-4-5
    end_lon, end_lat = to_4326_.transform(ex, ey)
    line_4326 = LineString([(start_lon, start_lat), (end_lon, end_lat)])

    margin = 5 * cell
    min_x, max_x = min(sx, ex) - margin, max(sx, ex) + margin
    min_y, max_y = min(sy, ey) - margin, max(sy, ey) + margin
    n_cols = math.ceil((max_x - min_x) / cell)
    n_rows = math.ceil((max_y - min_y) / cell)
    origin_x, origin_y = min_x, min_y + n_rows * cell  # borde norte, para from_origin

    arr = np.zeros((n_rows, n_cols), dtype="float64")
    for row in range(n_rows):
        for col in range(n_cols):
            X = origin_x + (col + 0.5) * cell
            Y = origin_y - (row + 0.5) * cell
            arr[row, col] = z0 + mx * X + my * Y
    transform = from_origin(origin_x, origin_y, cell, cell)
    raster_path = tmp_path / "plano_inclinado.tif"
    with rasterio.open(
        raster_path, "w", driver="GTiff", height=n_rows, width=n_cols, count=1,
        dtype="float64", crs="EPSG:9377", transform=transform,
    ) as dst:
        dst.write(arr, 1)

    monkeypatch.setattr(
        "elevation_profile.read_route",
        lambda path: SimpleNamespace(line=line_4326, n_vertices=2, length_km=0.5),
    )
    monkeypatch.setattr("elevation_profile.configure_remote_access", lambda cache_dir: None)
    monkeypatch.setattr(
        "elevation_profile.get_or_build_local_clip",
        lambda dem_key, bbox, cache_dir: raster_path,
    )

    df = build_profile(
        "ruta-sintetica.kml", paso_m=10.0, dem_keys=("cop30",),
        reference_dem="cop30", cache_dir=tmp_path,
    )

    # 1. sample_raster() no se desvía del plano analítico en más de la cota
    #    del vecino más cercano (por eje, no por hipotenusa: ver T14).
    analytic = z0 + mx * df["x_9377"] + my * df["y_9377"]
    bound_point = 0.5 * cell * (abs(mx) + abs(my)) + 1e-6
    assert (df["cop30"] - analytic).abs().max() <= bound_point

    # 2. pendiente_pct: primera fila en 0 (caso borde de R8); la pendiente
    #    promedio de extremo a extremo reproduce s dentro de la cota agregada.
    assert df["pendiente_pct"].iloc[0] == 0.0
    avg_slope = (df["cop30"].iloc[-1] - df["cop30"].iloc[0]) / df["distancia_m"].iloc[-1] * 100.0
    bound_avg = cell * (abs(mx) + abs(my)) / df["distancia_m"].iloc[-1] * 100.0 + 1e-6
    assert abs(avg_slope - s) <= bound_avg

    # 3. distancia_m acumulada coincide con la longitud euclidiana (500 m)
    #    dentro de 1 mm.
    assert df["distancia_m"].iloc[-1] == pytest.approx(500.0, abs=1e-3)


# --------------------------------------------------------------------------
# T22: integración con red real, sin mocks, contra los valores medidos el
# 2026-09-07 en (-74.5, 6.5) (requirements.md, R5 y R6).
# --------------------------------------------------------------------------

@pytest.mark.requiere_red
def test_real_dems_match_measured_values(tmp_path):
    cache_dir = tmp_path / "cache"
    configure_remote_access(cache_dir)
    # (-74.5, 6.5) exacto cae justo en el empate de la grilla de las cuatro
    # fuentes (1 arco-segundo para los tres DEM, 0.00025° para GLAD): a esa
    # precisión, "vecino más cercano" está mal definido y una resta en coma
    # flotante al construir el recorte puede desplazar la lectura a la celda
    # vecina (medido: cop30 pasó de 198.1 a 199.1 m antes de este ajuste; ver
    # el log de esta tarea). Se consulta 1e-5° al lado (~1 m), bien dentro de
    # la misma celda que el punto de referencia -- ninguna trayectoria real
    # cae jamás exactamente sobre una marca de arco-segundo entera.
    eps = 1e-5
    point = [(-74.5 + eps, 6.5 - eps)]
    bbox = (point[0][0] - 0.01, point[0][1] - 0.01, point[0][0] + 0.01, point[0][1] + 0.01)

    expected = {"cop30": 198.1, "cop90": 200.5, "srtm_gl1": 195.0}
    for dem_key, expected_value in expected.items():
        clip_path = get_or_build_local_clip(dem_key, bbox, cache_dir)
        measured = sample_raster(clip_path, point)[0]
        assert measured == pytest.approx(expected_value, abs=0.5), dem_key

    canopy_clip = get_or_build_local_clip("glad_canopy", bbox, cache_dir)
    canopy_measured = sample_raster(canopy_clip, point)[0]
    assert canopy_measured == pytest.approx(14.0, abs=0.5)
