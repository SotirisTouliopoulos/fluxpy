
# -*- coding: utf-8 -*-
#
# cobra documentation build configuration file, created by
# sphinx-quickstart on Wed Jun 13 19:17:34 2012.
#
# This file is execfile()d with the current directory set to its containing
# dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys
from os.path import dirname, join

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
SRC_PATH = join(dirname(dirname(__file__)))
sys.path.insert(0, SRC_PATH)

# -- General configuration ----------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "autoapi.extension",
    "nbsphinx",
]
# Document Python Code
autoapi_dirs = [join(SRC_PATH, "fluxpy")]
autoapi_add_toctree_entry = False


# Enable typehints
autodoc_typehints = "signature"

# Napoleon settings
napoleon_numpy_docstring = True

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "fluxpy"
copyright = "2024, Haris Zafeiropoulos"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
# This import has to be here.
from fluxpy import __version__ as release  # noqa: E402

version = ".".join(release.split(".")[:2])

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", ".ipynb_checkpoints"]

pygments_style = "sphinx"

# -- Options for HTML output --------------------------------------------------

mathjax_path = (
    "https://cdn.mathjax.org/mathjax/latest/"
    "MathJax.js?config=TeX-AMS-MML_HTMLorMML"
)
html_theme = "furo"

# -- Options for LaTeX output --------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    "papersize": "a4paper",
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    "preamble": r"\usepackage{amsmath,amssymb}",
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
    (
        "index",
        "fluxpy.tex",
        "fluxpy Documentation",
        "Haris Zafeiropoulos",
        "manual",
    )
]

# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ("index", "fluxpy", u"fluxpy Documentation", [u"Haris Zafeiropoulos"], 1)
]

# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        "index",
        "fluxpy",
        "fluxpy Documentation",
        "Haris Zafeiropoulos",
        "fluxpy",
        "A Python library for facilitating metabolic modeling analysis tasks",
        "Miscellaneous",
    )
]

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "http://docs.python.org/": None,
    "http://docs.scipy.org/doc/numpy/": None,
    "http://docs.scipy.org/doc/scipy/reference": None,
}
intersphinx_cache_limit = 10  # days to keep the cached inventories
