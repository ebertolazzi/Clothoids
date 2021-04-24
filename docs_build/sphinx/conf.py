# -*- coding: utf-8 -*-
import os

# -- Project information -----------------------------------------------------
exec(open("../project_common.py").read())

intersphinx_mapping = {
  'api-c': ('../sphinx_c/_build/html/objects.inv', None),
  'api-cpp': ('../sphinx_cpp/_build/html/objects.inv', None),
  'api-json': ('../sphinx_json/_build/html/objects.inv', None),
}

# If false, no module index is generated.
html_domain_indices = False

# If false, no index is generated.
html_use_index = False
