import math

def srtm_tile_names(bounds):
    """
    bounds = (minx, miny, maxx, maxy) in lon/lat
    """
    min_lon, min_lat, max_lon, max_lat = bounds

    tiles = []

    for lat in range(math.floor(min_lat), math.ceil(max_lat) + 1):
        for lon in range(math.floor(min_lon), math.ceil(max_lon) + 1):
            ns = "N" if lat >= 0 else "S"
            ew = "E" if lon >= 0 else "W"
            tile = f"{ns}{abs(lat):02d}{ew}{abs(lon):03d}.hgt"
            tiles.append(tile)

    return tiles


import requests
from pathlib import Path
from tqdm import tqdm
import warnings

def download_srtm_tiles(tile_names, output_dir):
    aws_base = "https://elevation-tiles-prod.s3.amazonaws.com/skadi"
    cgiar_base = "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF"

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    downloaded = []

    for tile in tqdm(tile_names, desc="Downloading SRTM tiles"):
        lat_band = tile[0:3]   # N08
        lon_band = tile[3:7]   # W074

        out_file = output_dir / f"{tile}.hgt.gz"

        if out_file.exists():
            downloaded.append(out_file)
            continue

        # --- Try AWS SRTM first ---
        aws_url = f"{aws_base}/{lat_band}/{lon_band}/{tile}.hgt.gz"
        r = requests.get(aws_url, stream=True)

        if r.status_code == 200:
            with open(out_file, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            downloaded.append(out_file)
            continue

        # --- Fallback: CGIAR SRTM (GeoTIFF, 5x5 deg) ---
        print(f"⚠ AWS missing {tile}, falling back to CGIAR SRTM")

        # CGIAR tiles are large, so we download by region
        # Colombia is covered by tiles like: srtm_12_05.zip
        warnings.warn(
                f"SRTM1 tile {tile} not available from AWS. Skipping."
                )
        continue
    


    return downloaded

def download_copernicus_dem(bounds, out_tif):
    """
    bounds = (minx, miny, maxx, maxy) in lon/lat
    """
    url = "https://portal.opentopography.org/API/globaldem"
    params = {
        "demtype": "COP30",
        "south": bounds[1],
        "north": bounds[3],
        "west": bounds[0],
        "east": bounds[2],
        "outputFormat": "GTiff",
    }

    r = requests.get(url, params=params, stream=True)
    if r.status_code != 200:
        raise RuntimeError("Failed to download Copernicus DEM")

    with open(out_tif, "wb") as f:
        f.write(r.content)



# Unzip and convert .hgt → GeoTIFF
import gzip
import shutil
import rasterio
from rasterio.transform import from_origin

def hgt_to_geotiff(hgt_gz_files, output_dir):
    output_dir = Path(output_dir)
    tifs = []

    for gz in hgt_gz_files:
        hgt_path = output_dir / gz.stem
        tif_path = output_dir / f"{gz.stem}.tif"

        if tif_path.exists():
            tifs.append(tif_path)
            continue

        # unzip
        with gzip.open(gz, "rb") as f_in, open(hgt_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        # read raw HGT
        size = 3601
        data = np.fromfile(hgt_path, dtype=">i2").reshape((size, size))

        # parse tile location
        lat = int(gz.stem[1:3]) * (1 if gz.stem[0] == "N" else -1)
        lon = int(gz.stem[4:7]) * (1 if gz.stem[3] == "E" else -1)

        transform = from_origin(lon, lat + 1, 1/3600, 1/3600)

        with rasterio.open(
            tif_path,
            "w",
            driver="GTiff",
            height=size,
            width=size,
            count=1,
            dtype=data.dtype,
            crs="EPSG:4326",
            transform=transform,
        ) as dst:
            dst.write(data, 1)

        tifs.append(tif_path)

    return tifs

#Mosaic tiles into one DEM

from rasterio.merge import merge
def mosaic_dem(tif_files, output_tif):
    srcs = [rasterio.open(t) for t in tif_files]
    mosaic, transform = merge(srcs)

    with rasterio.open(
        output_tif,
        "w",
        driver="GTiff",
        height=mosaic.shape[1],
        width=mosaic.shape[2],
        count=1,
        dtype=mosaic.dtype,
        crs="EPSG:4326",
        transform=transform,
    ) as dst:
        dst.write(mosaic[0], 1)

    for src in srcs:
        src.close()

#MAIN DRIVER (this is what you call)
import geopandas as gpd
import numpy as np
def build_srtm_for_line(geoshape_path, out_dir="data/srtm"):
    gdf = gpd.read_file(geoshape_path).to_crs(epsg=4326)
    geom = gdf.geometry.iloc[0]
    bounds = geom.bounds

    tiles = srtm_tile_names(bounds)
    gz_files = download_srtm_tiles(tiles, out_dir)
    tif_tiles = hgt_to_geotiff(gz_files, out_dir)

    output_tif = Path(out_dir) / "SRTM_Colombia.tif"
    mosaic_dem(tif_tiles, output_tif)

    return output_tif

#bounds = gdf.geometry.iloc[0].buffer(0.01).bounds
#download_copernicus_dem(bounds, "DEM_Colombia.tif")

dem_path = build_srtm_for_line("src//TL_Colombia_500kV.geojson")
print("DEM ready:", dem_path)