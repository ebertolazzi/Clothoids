#!/bin/bash

# This script is used to build the Sphinx documentation.

# Welcome the user
echo -e "[32mWelcome to the DRIVEWISE documentation build script!(B[m\n"

# Activate the conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sphinx

# Build the documentation
# Delete the _doxygen folder if it exists
if [ -d "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/source/_doxygen" ]; then
   rm -rf "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/source/_doxygen"
fi
doxygen "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/source/Doxyfile" > /dev/null
# Delete the _breathe folder if it exists
if [ -d "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/source/_breathe" ]; then
   rm -rf "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/source/_breathe"
fi
breathe-apidoc -o "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/source/_breathe" -m -p "CLOTHOIDS" -q "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/source/_doxygen/xml" > /dev/null
sphinx-build -b html "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/source" "/Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/build/html" > /dev/null

# Deactivate the conda environment
conda deactivate

echo -e "\n[32mBuilt the HTML documentation in /Users/enrico/Ricerca/develop/PINS/pins-mechatronix/LibSources/submodules/Clothoids/doc_generate/docs/build/html.(B[m\n"
echo "[32mEnjoy!(B[m"
