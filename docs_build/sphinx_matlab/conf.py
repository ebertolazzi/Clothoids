# -*- coding: utf-8 -*-

# pip3 install exhale
# pip3 install breathe
# pip3 install m2r2
# pip3 install sphinxcontrib-email
# pip3 install cloud_sptheme

import os

# -- Project information -----------------------------------------------------
exec(open("../project_common.py").read())

#from project_common import *

# Setup the breathe extension
breathe_projects        = { project+"_matlab": "../xml-matlab" }
breathe_default_project = project+"_matlab"

extensions.append('exhale');
extensions.append('breathe');

# Setup the exhale extension
exhale_args = {
  # These arguments are required
  "containmentFolder":     "./api-matlab",
  "rootFileName":          "library_root.rst",
  "rootFileTitle":         "MATLAB API",
  "doxygenStripFromPath":  "..",
  # Suggested optional arguments
  "createTreeView":        True,
  # TIP: if using the sphinx-bootstrap-theme, you need
  "treeViewIsBootstrap":   False,
  "exhaleExecutesDoxygen": True,
  #"exhaleDoxygenStdin":    "INPUT = ../../src"
  "exhaleDoxygenStdin":
'''
        EXTRACT_ALL         = YES
        SOURCE_BROWSER      = YES
        EXTRACT_STATIC      = YES
        HIDE_SCOPE_NAMES    = NO
        CALLER_GRAPH        = YES
        GRAPHICAL_HIERARCHY = YES
        HAVE_DOT            = YES
        QUIET               = NO
        INPUT               = ../../toolbox/lib
        GENERATE_TREEVIEW   = YES
        XML_OUTPUT          = xml-matlab
        SHORT_NAMES         = YES
        IMAGE_PATH          = ../images

        XML_PROGRAMLISTING    = YES
        RECURSIVE             = YES
        FULL_PATH_NAMES       = YES
        ENABLE_PREPROCESSING  = YES
        MACRO_EXPANSION       = YES
        SKIP_FUNCTION_MACROS  = NO
        EXPAND_ONLY_PREDEF    = NO
        INHERIT_DOCS          = YES
        INLINE_INHERITED_MEMB = YES
        EXTRACT_PRIVATE       = YES
        PREDEFINED           += protected=private
        EXTENSION_MAPPING     = .m=C++
        FILE_PATTERNS         = *.m
        FILTER_PATTERNS       = *.m=./m2cpp.pl
        GENERATE_HTML         = NO
''',
  "lexerMapping": { r".*\.m": "MATLAB" }
}

html_theme_options['logotarget'] = "../index"
html_theme_options['roottarget'] = "../index"