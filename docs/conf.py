# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'RIA'
copyright = '2026, Richard Howey'
author = 'Richard Howey'
release = '1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [    
    'sphinx_rtd_theme',
    'sphinxcontrib.bibtex',
]

# Fix for incompatible versions where things have changed location.
# https://stackoverflow.com/questions/59636631/importerror-cannot-import-name-mutablemapping-from-collections

import collections
collections.Iterable = collections.abc.Iterable
collections.Mapping = collections.abc.Mapping
collections.MutableSet = collections.abc.MutableSet
collections.MutableMapping = collections.abc.MutableMapping

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']
# html_js_files = ['custom.js']
# html_favicon = 'images/favicon.ico'

# templates_path = ['_templates']

html_theme_options = {  
    'logo_only': False,
     #'html_logo': ,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#FF9999',  # Navbar background color
}
