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

A notable feature of the Camera class is its support for pixel defects, which are documented in the {doc}`pixel defects <cabaret.camera>` section.

## Generating images

The simulation of images uses two important concepts:

- **Filters**: Photometric bands (e.g., G, R, I) that determine the wavelength range in which the catalog fluxes are extracted for the simulation.
- **GaiaQuery**: Provides methods to query star catalogs such as Gaia and 2MASS, returning tables or ready-to-use Sources objects for simulations. The TAP service endpoint can be selected via `GaiaTAPSource` (or by passing ``"GAIA"`` / ``"VIZIER"`` as a string).
- **Sources**: Representations of stars with positions and fluxes, either queried from catalogs or provided directly.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :template: class.rst

   cabaret.Filters
   cabaret.GaiaQuery
   cabaret.GaiaTAPSource
   cabaret.Sources
```

```{toctree}
:maxdepth: 1
:hidden:

cabaret.camera
```
