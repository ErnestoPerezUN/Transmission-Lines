#Code developed by ChatGPT with the following prompt:
# make a code to get the profile for a transmission line in Colombia using SRTM, The pathe is provided by a geoshape

import os as os
import geopandas as gpd
from shapely.geometry import LineString
# Load line geometry (GeoJSON or Shapefile)
print(os.getcwd())
gdf = gpd.read_file("src//TL_Colombia_500kV.geojson")
# Ensure WGS84

gdf = gdf.to_crs(epsg=4326)

line = gdf.geometry.iloc[0]

#Reproject to UTM (automatic for Colombia)
from pyproj import CRS
from pyproj.database import query_utm_crs_info
from pyproj.aoi import AreaOfInterest


lon = line.centroid.x
lat = line.centroid.y

utm_crs_list = query_utm_crs_info(
    datum_name="WGS 84",
    area_of_interest=AreaOfInterest(
        west_lon_degree=lon,
        south_lat_degree=lat,
        east_lon_degree=lon,
        north_lat_degree=lat,
    ),
)

utm_crs = CRS.from_epsg(utm_crs_list[0].code)

#utm_crs = CRS.from_user_input(
#    CRS.from_epsg(4326).utm_zone(line.centroid.x, line.centroid.y)
#)

gdf_utm = gdf.to_crs(utm_crs)
line_utm = gdf_utm.geometry.iloc[0]

print(type(line_utm))
print(line_utm.geom_type)


# Convert to a LineString
from centerline.geometry import Centerline
from shapely.geometry import Polygon, MultiPolygon

def polygon_to_centerline(poly):
    if isinstance(poly, MultiPolygon):
        poly = max(poly.geoms, key=lambda g: g.area)

    if not isinstance(poly, Polygon):
        raise TypeError("Input must be Polygon or MultiPolygon")

    cl = Centerline(poly)
    return cl.geometry

line_utm = polygon_to_centerline(line_utm)

print(type(line_utm))
print(line_utm.geom_type)
#Interpolate points alog line
import numpy as np

def sample_line(line, step_m=30):
    distances = np.arange(0, line.length, step_m)
    points = [line.interpolate(d) for d in distances]
    return distances, points

dist_m, points_utm = sample_line(line_utm, step_m=30)

#convert sample points back to lat/lon
from shapely.ops import transform
from pyproj import Transformer

transformer = Transformer.from_crs(utm_crs, 4326, always_xy=True)

points_ll = [
    transform(transformer.transform, p) for p in points_utm
]

import requests
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



bounds = gdf.geometry.iloc[0].buffer(0.01).bounds
download_copernicus_dem(bounds, "SRTM_Colombia.tif")


#sample elevation
import rasterio

def sample_elevation(raster_path, points):
    elevations = []
    with rasterio.open(raster_path) as src:
        for p in points:
            row, col = src.index(p.x, p.y)
            elev = src.read(1)[row, col]
            elevations.append(elev)
    return np.array(elevations)

srtm_file = "SRTM_Colombia.tif"
elevation_m = sample_elevation(srtm_file, points_ll)

distance_km = dist_m / 1000.0

#plot profile

import matplotlib.pyplot as plt

plt.figure(figsize=(12,4))
plt.plot(distance_km, elevation_m, linewidth=1.5)
plt.xlabel("Distance along line [km]")
plt.ylabel("Elevation [m]")
plt.title("Transmission Line Elevation Profile (SRTM 30 m)")
plt.grid(True)
plt.tight_layout()
plt.show()

# Export
import pandas as pd

profile = pd.DataFrame({
    "distance_m": dist_m,
    "elevation_m": elevation_m,
    "lon": [p.x for p in points_ll],
    "lat": [p.y for p in points_ll]
})

profile.to_csv("elevation_profile.csv", index=False)