# Pixel defects

Pixel defects are imperfections in camera sensors that can affect the quality of astronomical images. These defects may introduce noise, reduce sensitivity, or cause artifacts in the data. The classes below model different types of pixel defects, allowing you to simulate and study their effects in synthetic images generated with **cabaret**.

Use these classes to add realistic sensor behaviors to your simulations or to test how your analysis pipeline handles imperfect data.

You can create pixel defect objects directly, or more conveniently, by passing a dictionary of defect configurations to the `pixel_defects` argument when constructing a `Camera`. Each dictionary entry specifies the type and parameters of a defect, and the camera will automatically instantiate the appropriate defect classes.
See the `Camera` class for more details and examples on configuring pixel defects using dictionaries.

```{eval-rst}
.. currentmodule:: cabaret.camera

.. autosummary::
   :template: class.rst
   :toctree: generated

   ConstantPixelDefect
   ColumnPixelDefect
   RandomNoisePixelDefect
   QuantumEfficiencyMapPixelDefect
   ReadoutSmearPixelDefect
```

```{toctree}
:maxdepth: 1
:hidden:


```
