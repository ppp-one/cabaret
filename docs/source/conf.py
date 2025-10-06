import importlib.metadata
import os
import sys

# Add the project source directory to the path so that autodoc can find the modules
sys.path.insert(0, os.path.abspath("../../src"))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "cabaret"
copyright = "2025, Peter Pedersen, Lionel Garcia, David Degen"
author = "Peter Pihlmann, Lionel Garcia, David Degen"

version = importlib.metadata.version("cabaret")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_nb",
    # "numpydoc",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx_design",
    "IPython.sphinxext.ipython_console_highlighting",
]

nb_render_image_options = {"align": "center"}

# Add mappings for intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_css_files = ["style.css"]
html_short_title = "cabaret"
html_title = f"{html_short_title}"

html_context = {
    "default_mode": "light",
}

html_theme_options = {
    "path_to_docs": "docs",
    "repository_url": "https://github.com/ppp-one/cabaret",
    "repository_branch": "main",
    "use_repository_button": True,
    "use_fullscreen_button": False,
    "use_issues_button": True,
    "use_download_button": False,
    "show_navbar_depth": 1,
    "navbar_end": ["navbar-icon-links"],
    "collapse_navigation": True,
}


# Auto-generate API documentation
autodoc_member_order = "bysource"
autodoc_typehints = "signature"
autoclass_content = "class"
autodoc_default_options = {"no-value": True}
autodoc_preserve_defaults = False

myst_enable_extensions = ["dollarmath", "colon_fence"]

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_use_ivar = True
napoleon_attr_annotations = True
