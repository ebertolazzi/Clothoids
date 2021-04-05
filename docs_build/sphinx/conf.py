# -*- coding: utf-8 -*-
import os
from past.builtins import execfile

# -- Project information -----------------------------------------------------
execfile('../project_common.py')

intersphinx_mapping = {
  'api-c': ('../sphinx_c/_build/html/objects.inv', None),
  'api-cpp': ('../sphinx_cpp/_build/html/objects.inv', None),
  'api-json': ('../sphinx_json/_build/html/objects.inv', None),
}

# If false, no module index is generated.
html_domain_indices = False

# If false, no index is generated.
html_use_index = False
