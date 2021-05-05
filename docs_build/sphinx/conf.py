# -*- coding: utf-8 -*-
import os
from pathlib import Path

# -- Project information -----------------------------------------------------
exec(open("../project_common.py").read())


rst_prolog = ".. |xml| replace:: %s\n" % (project)

extensions.append('breathe');
extensions.append('exhale');

breathe_projects = {
  "doc_matlab": "_doxygen/"+"doc_matlab/xml-matlab",
  "doc_cpp":    "_doxygen/"+"doc_cpp/xml-cpp",
}

breathe_default_project = "doc_cpp"

dir_path_cpp = os.path.dirname(os.path.realpath(__file__))+"../../../src"
dir_path_cpp = Path(dir_path_cpp).resolve()

dir_path_matlab = os.path.dirname(os.path.realpath(__file__))+"../../../toolbox/lib"
dir_path_matlab = Path(dir_path_matlab).resolve()

doxygen_common_stdin = """
        EXTRACT_ALL         = YES
        SOURCE_BROWSER      = YES
        EXTRACT_STATIC      = YES
        HIDE_SCOPE_NAMES    = NO
        CALLER_GRAPH        = YES
        GRAPHICAL_HIERARCHY = YES
        HAVE_DOT            = YES
        QUIET               = NO
        GENERATE_TREEVIEW   = YES
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
        GENERATE_HTML         = NO
"""

exhale_projects_args = {
  # Third Party Project Includes
  "doc_matlab": {
    'verboseBuild':          True,
    "rootFileName":          "root.rst",
    "createTreeView":        True,
    "exhaleExecutesDoxygen": True,
    "doxygenStripFromPath":  str(dir_path_matlab),
    "exhaleDoxygenStdin":   '''
        INPUT               = ../../toolbox/lib
        PREDEFINED         += protected=private
        XML_OUTPUT          = xml-matlab
        EXTENSION_MAPPING   = .m=C++
        FILE_PATTERNS       = *.m
        FILTER_PATTERNS     = *.m=./m2cpp.pl
'''+doxygen_common_stdin,
    "containmentFolder":    os.path.realpath('./api-matlab'),
    "rootFileTitle":        "MATLAB API",
    "lexerMapping": { r".*\.m": "MATLAB" }
  },

  "doc_cpp": {
    'verboseBuild':          True,
    "rootFileName":          "root.rst",
    "createTreeView":        True,
    "exhaleExecutesDoxygen": True,
    "doxygenStripFromPath":  str(dir_path_cpp),
    "exhaleDoxygenStdin":   '''
        INPUT               = ../../src
        PREDEFINED         += protected=private
        XML_OUTPUT          = xml-cpp
'''+doxygen_common_stdin,
    "containmentFolder":    os.path.realpath('./api-cpp'),
    "rootFileTitle":        "C++ API",
  }
}

cpp_index_common_prefix = ['G2lib::']

# If false, no module index is generated.
html_domain_indices = True

# If false, no index is generated.
html_use_index = True
