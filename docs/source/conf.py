# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# so that sphinx can locate the pipeline:
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'P4ward'
html_title = 'P4ward'
copyright = '2024, Paula Jofily'
author = 'Paula Jofily'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx_design'
]

templates_path = ['_templates']
exclude_patterns = []

# add the tutorial html files to be copied to the build folder as well
html_extra_path = ['../../tutorial/plots-mz1.html']

autodoc_mock_imports = [
    "rdkit",
    "numpy",
    "plotly",
    "pandas",
    "openmm",
    "openbabel",
    "multiprocessing",
    "func_timeout"
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_permalinks_icon = '<span>#</span>'
html_theme = 'sphinxawesome_theme'

# Select theme for both light and dark mode
pygments_style = "default"
# Select a different theme for dark mode
pygments_style_dark = "default"

# html_theme = 'furo'
html_static_path = ['_static']
