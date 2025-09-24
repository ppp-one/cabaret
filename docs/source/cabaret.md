# API

## The Observatory
The `Observatory` class is the main entry point for simulating astronomical images. 

```{eval-rst}
.. autosummary::
   :toctree: generated
   :template: class.rst

   cabaret.Observatory
```

### The Observatory Devices
The `Observatory` class encapsulates the configuration of all key devices at an observatory, as well as its site.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :template: class.rst

   cabaret.Camera
   cabaret.Focuser
   cabaret.Telescope
   cabaret.Site
```

A special feature is the `Camera` class are pixel defects, which are covered [`here`](cabaret.camera.html).

## Generating images

The simulation of images uses two important concepts:

- **Filters**: Photometric bands (e.g., G, R, I) that determine the wavelength range in which the catalog fluxes are extracted for the simulation.
- **Sources**: Representations of stars with positions and fluxes, either queried from catalogs or provided directly.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :template: class.rst

   cabaret.Filters
   cabaret.Sources
```

```{toctree}
:maxdepth: 1
:hidden:

cabaret.camera
```
