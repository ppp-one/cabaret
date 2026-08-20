# Installation

To install **cabaret**, you need Python 3.11 or newer.

Cabaret can be installed using pip
```bash
pip install cabaret          # basic installation
pip install "cabaret[plot]"  # with matplotlib for plotting
```
or from sources
```bash
git clone https://github.com/ppp-one/cabaret.git
cd cabaret
pip install .
```

For development, we recommend [uv](https://docs.astral.sh/uv/)
```bash
git clone https://github.com/ppp-one/cabaret.git
cd cabaret
uv sync --dev          # installs the package + dev + test dependencies
uv sync --group docs   # additionally install documentation dependencies
```
