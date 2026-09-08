# Transmission Lines — código de los cursos de Modelado y Diseño de Líneas de Transmisión

Código académico de Ernesto Pérez para dos cursos de la Universidad Nacional de Colombia:

- **Modelado de Líneas de Transmisión**: parámetros con tierra finita, campos vs. circuitos (FDTD), campo eléctrico.
- **Diseño de Líneas de Transmisión**: capacidad de corriente (IEEE 738), gradiente superficial, perfil topográfico,
  distancias desde fotografías.

> **Sin garantía.** Este código existe para enseñar. No ha sido validado contra mediciones de campo y **no debe
> usarse para decisiones de diseño, operación o seguridad de instalaciones reales** sin verificación independiente.
> Lee [`DISCLAIMER.md`](DISCLAIMER.md) (español e inglés) antes de usarlo.
>
> **Licencia**: por definir; la intención es una licencia de código abierto permisiva una vez se verifique el
> estatuto de propiedad intelectual de la UNAL. Mientras tanto, todos los derechos reservados. Los datos de
> OpenStreetMap incluidos van bajo ODbL (ver [`NOTICE.md`](NOTICE.md)).

*English: academic code, provided as is, without warranty of any kind; not for real-world design. See
[`DISCLAIMER.md`](DISCLAIMER.md).*

## Contenido

| Ruta | Curso | Qué hace | Tests |
|---|---|---|---|
| `src/parameters/` + `notebooks/parameters/` | Modelado | Parámetros Z y Y de una línea con tierra finita (profundidad compleja de Deri), clases `cable`, `conductor`, `Tower` | — |
| `src/fdtd/` + `examples/voltage_source_on_line.py` | Modelado | Propagación electromagnética en un conductor sobre tierra (FDTD 2D) comparada con los modelos de circuito concentrado y distribuido | 15 |
| `src/electric_field/` | Modelado | Campo eléctrico en la superficie de conductores por el método de simulación de cargas | — |
| `notebooks/current_capacity/` | Diseño | Corriente máxima de un conductor desnudo según IEEE 738-2012 (aproximación estática) | — |
| `notebooks/Emax_TL/` | Diseño | Gradiente superficial de un haz de conductores (Kuffel, ec. 2.9) | — |
| `src/line_profile/` | Diseño | Perfil de elevación de una línea a partir de SRTM / Copernicus DEM (requiere `requirements-geo.txt`) | — |
| `src/tower_distances/` | Diseño | Distancias entre conductores a partir de una fotografía de la torre (requiere `requirements-vision.txt`) | — |

Los archivos generados con asistentes de IA lo declaran en su encabezado; la lista completa está en `NOTICE.md`.

## Instalación

```bash
python -m venv .venv
.venv/Scripts/activate            # Windows; en Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt   # núcleo (numpy, scipy, matplotlib, pandas, pytest)
pip install -r requirements-geo.txt     # opcional: perfil topográfico
pip install -r requirements-vision.txt  # opcional: distancias desde fotografías
```

Probado con Python 3.12. Los notebooks se guardan sin salidas (`nbstripout`, configurado en `.gitattributes`):
ejecútalos para ver las figuras.

## Campos vs. circuitos: `examples/voltage_source_on_line.py`

Conecta una fuente de tensión ideal a un conductor de 1 m a 10 cm sobre un plano de tierra y muestra, en cuatro
etapas, por qué eso no es un circuito concentrado:

1. la onda sale de la fuente, cruza la línea a `c` y se duplica en el extremo abierto;
2. un solo `L`-`C` concentrado falla cuando el flanco de la fuente es rápido; el criterio de validez se mide, no se
   afirma (el error cae como `(t_r/T)^-2` y se satura por debajo de `t_r ~ T`);
3. una escalera de `N` secciones converge a la solución de campos a lo largo de la línea, pero solo en el límite;
4. los números que el modelo de circuito no puede producir por sí mismo.

La geometría es una placa paralela a propósito: `L' = mu0*h/w` y `C' = eps0*w/h` son exactos y la comparación no
tiene parámetros ajustados. Cualquier desacuerdo es una falla del modelo de circuito, no una mala estimación de
`L'`. La simulación recupera `v = 0.9999 c` y `Z0 = 37.67` ohm contra un valor exacto de 37.67, ambos medidos
desde los campos.

El extremo abierto es el caso interesante: ninguna escalera lo reproduce para ningún `N` (el error se estanca en vez
de converger) porque el extremo tiene campo de borde y radia. Ajustar la capacitancia faltante da unos 10 pF, que
es `C' * h`: el extremo se comporta como una altura de conductor adicional de línea, y solo una solución de campos
lo dice.

```bash
python examples/voltage_source_on_line.py              # interactivo
python examples/voltage_source_on_line.py --save out   # PNG + GIF animado
python -m pytest tests -q                              # aserciones de física
```

Los módulos del solver están en `src/fdtd/`. Los scripts `examples/FDTD1D_*.py` y `examples/FDTD2D_*.py` son
bocetos anteriores sin relación y se dejan como estaban.

## Para contribuir o mantener

Cada módulo con tests trae al menos una aserción de física contra un valor cerrado (ver `tests/`); las ecuaciones
citan su fuente (norma, libro o artículo) en el docstring de la función que las implementa; los notebooks se
guardan sin salidas. Antes de agregar una dependencia, revisa si ya está en `requirements.txt` o en uno de los
extras (`requirements-geo.txt`, `requirements-vision.txt`).
