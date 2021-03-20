# -*- coding: utf-8 -*-

# pip3 install recommonmark
# pip3 install exhale
# pip3 install breathe
# pip3 install pydata-sphinx-theme
# pip3 install sphinx-markdown-parser
# pip3 install pymdown-extensions
# pip3 install m2r2
# pip3 install sphinxcontrib-email
# pip3 install furo
#
# https://pradyunsg.me/furo/

# The master toctree document.
master_doc = 'index'

# -- Project information -----------------------------------------------------

project   = 'Clothoids'
copyright = '2021, Enrico Bertolazzi and Marco Frego'
author    = ':email:`Enrico Bertolazzi <enrico.bertolazzi@unitn.it>`'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The `extensions` list should already be in here from `sphinx-quickstart`
extensions = [
  'breathe',
  'exhale',
  #'recommonmark', # non funziona
  'm2r2',          # funziona!
  'sphinx.ext.autodoc',
  'sphinx.ext.doctest',
  'sphinx.ext.coverage',
  'sphinx.ext.mathjax',
  'sphinx.ext.ifconfig',
  'sphinx.ext.githubpages',
  'sphinx.ext.intersphinx',
  'sphinxcontrib.email'
]

source_suffix = ['.rst', '.md']

# Setup the breathe extension
breathe_projects = {
  "Clothoids": "../xml"
}
breathe_default_project = "Clothoids"

# Setup the exhale extension
exhale_args = {
  # These arguments are required
  "containmentFolder":     "./api",
  "rootFileName":          "library_root.rst",
  "rootFileTitle":         "Clothoids API",
  "doxygenStripFromPath":  "..",
  # Suggested optional arguments
  "createTreeView":        True,
  # TIP: if using the sphinx-bootstrap-theme, you need
  # "treeViewIsBootstrap": True,
  "exhaleExecutesDoxygen": True,
  #"exhaleDoxygenStdin":    "INPUT = ../../src"
  "exhaleDoxygenStdin":
'''
        EXTRACT_ALL       = YES
        SOURCE_BROWSER    = YES
        EXTRACT_STATIC    = YES
        HIDE_SCOPE_NAMES  = YES
        QUIET             = YES
        INPUT             = ../../src ../../src/Clothoids
        EXAMPLE_RECURSIVE = YES
        GENERATE_TREEVIEW = YES
'''
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

pygments_style      = "sphinx"
pygments_dark_style = "monokai"

#html_theme = 'pydata_sphinx_theme'
html_theme = 'furo'
html_logo  = '../logo.png'

email_automode = True
