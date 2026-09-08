# Avisos de terceros, datos y origen del contenido

Este archivo registra todo lo que **no** es obra original de Ernesto Pérez, o que tiene condiciones de uso propias.
Cualquier archivo nuevo de terceros debe agregar una fila aquí.

## Datos

| Archivo | Origen | Condiciones |
|---|---|---|
| `src/TL_Colombia_500kV.geojson` | Exportado de OpenStreetMap con overpass-turbo (2026-02-17) | © OpenStreetMap contributors, licencia [ODbL 1.0](https://opendatacommons.org/licenses/odbl/). Cualquier redistribución o derivado de estos datos debe conservar esta atribución. |

## Imágenes

| Archivo | Origen | Condiciones |
|---|---|---|
| `src/tower_distances/img/100kV.png` | Fotografía de Ernesto Pérez (origen por confirmar) | Misma licencia que el repositorio |
| `src/tower_distances/img/Tower_500kV.png` | Origen por confirmar; si resulta ser de Google Street View, se retira | — |

Otras fotografías del mismo directorio, con marca de agua de Google Street View, no se incluyen en este
repositorio porque las condiciones de Google Maps no permiten redistribuirlas.

## Código de terceros

| Archivo | Origen | Condiciones |
|---|---|---|
| (ninguno) | | |

Un script para medir distancias entre objetos en una fotografía, adaptado de un artículo de PyImageSearch
(Adrian Rosebrock, 2016), no se incluye en este repositorio porque ese código no tiene licencia de
redistribución. `src/tower_distances/Distance_measure.py` cubre el mismo fin con código propio.

## Código generado con asistentes de IA

Sin reclamo de terceros; se licencia como el resto del repositorio. Cada archivo lo declara en su encabezado:

- `notebooks/parameters/Calc_inductancia_ChatGPT.py`, `notebooks/Emax_TL/Emax_bundle.py` (ChatGPT).
- `src/tower_distances/Distance_measure.py` (ChatGPT).
- `examples/FDTD1D_propagating.py`, `examples/FDTD2D_propagating.py` (ChatGPT, bocetos sin verificar).
- `src/fdtd/`, `examples/voltage_source_on_line.py`, `tests/test_tem_line.py`, `src/line_profile/` (escritos con
  asistencia de Claude Code, verificados con tests de física).

## Ecuaciones y normas citadas

Se implementan ecuaciones; no se reproduce texto, tablas ni figuras de las fuentes.

- IEEE Std 738-2012, *Calculating the Current-Temperature Relationship of Bare Overhead Conductors*
  (`notebooks/current_capacity/`).
- J. A. Martínez-Velasco (ed.), *Power System Transients: Parameter Determination*, CRC Press, 2009
  (`src/parameters/`).
- E. Kuffel, W. S. Zaengl, J. Kuffel, *High Voltage Engineering Fundamentals*, 2.ª ed., ecuación 2.9
  (`notebooks/Emax_TL/`).
