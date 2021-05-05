# pip install sphinxcontrib-email
# pip install importlib_metadata
# pip install cloud_sptheme
# pip install exhale
# pip install sphinxcontrib-pyexec

import exhale
import os
import re
import glob
import os.path
from pprint import pprint

def do_filter(project, dir, filter = "/**/*.rst", regex=None):
  if regex is None:
    regex = re.compile(r"^(\s*\.\.\s*doxygen.*::.*)$", flags=re.IGNORECASE)
  for fl in glob.glob(dir+filter,recursive=True):
    #print(fl)
    with open(fl) as fp:
      file_data = fp.read()
    with open(fl, "w") as fp:
      for line in file_data.splitlines():
        match = regex.match(line)
        if match:
          fp.write(f"{match.group(1)}\n   :project: {project}\n")
        else:
          fp.write(line+'\n')

saved_exhale_environment_ready = exhale.environment_ready

def exhale_environment_ready(app):
  exhale_projects_args    = dict(app.config._raw_config['exhale_projects_args'])
  breathe_projects        = dict(app.config._raw_config['breathe_projects'])
  breathe_default_project = app.config._raw_config['breathe_default_project']

  for project in breathe_projects:
    app.config.breathe_default_project = project
    app.config._raw_config['rst_prolog'] = ".. |xml| replace:: %s\n" % (project)
    os.makedirs(breathe_projects[project], exist_ok=True)

    project_exhale_args    = exhale_projects_args.get(project, {})
    app.config.exhale_args = dict(project_exhale_args)
    #app.config.exhale_args["containmentFolder"] = os.path.realpath(app.config.exhale_args["containmentFolder"])
    print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
    print("="*75)
    print(project)
    print("-"*50)
    pprint(app.config.exhale_args)
    print("="*75)

    saved_exhale_environment_ready( app )

    dir = app.config.exhale_args["containmentFolder"];
    ##os.system("ruby ../filter_exhale_breathe.rb %s %s" % (dir, project));
    do_filter(project,dir)

    print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

  app.config.breathe_default_project = breathe_default_project

exhale.environment_ready = exhale_environment_ready

project   = 'Clothoids'
copyright = '2021, Enrico Bertolazzi and Marco Frego'
author    = ':email:`Enrico Bertolazzi <enrico.bertolazzi@unitn.it>` and :email:`Marco Frego <marco.frego@unibz.it>`'
version   = os.popen('git describe --tags --abbrev=0').read()
##release   = '1.0'

# The master toctree document.
master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = [ '../_templates' ]
html_static_path = ['../_static']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The `extensions` list should already be in here from `sphinx-quickstart`
extensions = [
  # standard sphinx extensions
  'sphinx.ext.autodoc',
  'sphinx.ext.todo',
  #'sphinx.ext.doctest',
  #'sphinx.ext.coverage',
  'sphinx.ext.mathjax',
  #'sphinx.ext.ifconfig',
  #'sphinx.ext.githubpages',
  #'sphinx.ext.intersphinx',
  #'sphinx.ext.graphviz',
  #'sphinx.ext.inheritance_diagram',

  # 3rd party extensions
  #'sphinxcontrib.fulltoc',
  'sphinxcontrib.email',

  # cloud's extensions
  'cloud_sptheme',
  'cloud_sptheme.ext.autodoc_sections',
  'cloud_sptheme.ext.relbar_links',
  'cloud_sptheme.ext.escaped_samp_literals',
  'cloud_sptheme.ext.issue_tracker',
  'cloud_sptheme.ext.table_styling',
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
}

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
##html_sidebars = {'**': ['searchbox.html', 'globaltoc.html']}
html_sidebars = { '**': ['searchbox.html', 'globaltoc.html'] }

# relbar_links = [("genindex2","MATLAB index")]

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

html_copy_source     = False
html_show_sourcelink = False
