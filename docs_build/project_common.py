import os

project   = 'Clothoids'
copyright = '2021, Enrico Bertolazzi and Marco Frego'
author    = ':email:`Enrico Bertolazzi <enrico.bertolazzi@unitn.it>` and :email:`Marco Frego <marco.frego@unibz.it>`'
version   = os.popen('git describe --tags --abbrev=0').read()
##release   = '1.0'

# The master toctree document.
master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The `extensions` list should already be in here from `sphinx-quickstart`
extensions = [
  #'breathe',
  #'exhale',
  'm2r2',          # funziona!
  # standard sphinx extensions
  'sphinx.ext.autodoc',
  'sphinx.ext.todo',

  # 3rd party extensions
  'sphinxcontrib.fulltoc',
  #'sphinx.ext.viewcode',  # mainly to make sure layout works properly

  # cloud's extensions
  'cloud_sptheme',
  'cloud_sptheme.ext.autodoc_sections',
  'cloud_sptheme.ext.relbar_links',
  'cloud_sptheme.ext.escaped_samp_literals',
  'cloud_sptheme.ext.issue_tracker',
  'cloud_sptheme.ext.table_styling',
  #'cloud_sptheme.ext.role_index',  # NOTE: used only to provide example role index

  #'sphinx.ext.doctest',
  #'sphinx.ext.coverage',
  'sphinx.ext.mathjax',
  #'sphinx.ext.ifconfig',
  #'sphinx.ext.githubpages',
  #'sphinx.ext.intersphinx',
  'sphinx.ext.graphviz',
  'sphinx.ext.inheritance_diagram',
  #'guzzle_sphinx_theme',
  #'sphinx_typo3_theme',
  'sphinxcontrib.email'
]

source_suffix = ['.rst', '.md']

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

html_theme = 'cloud'
html_logo  = '../logo.png'

email_automode       = True
autodoc_member_order = 'bysource'

html_theme_options = {
  "lighter_header_decor" : False,
  "borderless_decor"     : False,
  "bodyfont"             : "Arial, sans-serif",
  "headfont"             : "Arial, sans-serif",

  #styling for document body
  "bgcolor"         : "#f8f8f8",
  "linkcolor"       : "#006906",

  #styling for document headers
  "headlinkcolor"   : "#327438",

  #styling for section headers
  "sectiontextcolor"  : "inherit",
  "sectionbgcolor"    : "#75c47c",
  "sectiontrimcolor"  : "rgba(0,0,0,.1)",
  "rubricbgcolor"     : "#d2e7d0",
  "rubric_trim_color" : "rgba(0,0,0,0.05)",

  "object_default_color"   : "#e4e4e4",
  "object_function_color"  : "#eefbff",
  "object_class_color"     : "#fff3df",
  "object_attribute_color" : "#ffd5ff",
  "object_exception_color" : "#e9ffd0",

  #styling for footer / html background
  "footerbgcolor" : "#6f6700",

  #styling for sidebar
  "sidebarbgcolor"   : "#ededed",
  "sidebarlinkcolor" : "#006906",
  "sidebarhighcolor" : "#FFF5DD",
  "bodytrimcolor"    : "rgba(0,0,0,.15)",

  #styling for top & bottom relbars
  "relbarbgcolor" : "#57A75E",

  # code blocks
  "codebgcolor"   : "#e8ffe6",
  "codetrimcolor" : "#129100",

  # admonitions
  "admonition_note_color"       : "#D9E4F1",
  "admonition_warning_color"    : "#EBC5A7",
  "admonition_seealso_color"    : "#eeeeee",
  "admonition_deprecated_color" : "#ffebab",
  "admonition_todo_color"       : "#eeeeee",

  # inline literals
  "quotebgcolor"   : "rgba(0,0,0,.06)",
  "quotetrimcolor" : "transparent",

  # index page
  "index_category_color" : "#999999"
}

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
html_sidebars = {'**': ['searchbox.html', 'globaltoc.html']}

# If false, no module index is generated.
# html_domain_indices = False

# If false, no index is generated.
# html_use_index = False

# If true, the index is split into individual pages for each letter.
# html_split_index = True

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "%s v%s" % (project, version)

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "%s" % (project)
