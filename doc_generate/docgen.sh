#!/bin/bash

# This script is used to setup the Sphinx documentation environment following
# the DRIVEWISE documentation guidelines.

#1# Colours
GREEN=$(tput setaf 2)
YELLOW=$(tput setaf 3)
RED=$(tput setaf 1)
NORMAL=$(tput sgr0)

#1# Functions

# Add a spinning indicator to show the user this command is running
spinner() {
    local pid=$1
    local delay=0.1
    local spinstr='|/-\'
    tput civis # hide cursor
    printf "$2" # print the message
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
    printf " ${GREEN}Done!${NORMAL}\n\n"
    tput cnorm # show cursor
}

#1# Welcome message

echo -e "${GREEN}Welcome to the DRIVEWISE Documentation Generator!${NORMAL}\n"

echo -e "I will now ask you a few questions to setup the environment.\n"

# Check that Python is installed
if ! command -v python &> /dev/null
then
    echo -e "${RED}Python is not installed, please install it before continuing.${NORMAL}"
    exit
fi

#1# Conda setup

#2# Ask the user if they want to use a conda environment

USE_CONDA=true
while true; do
    read -p "Do you want to use a conda environment? ([y]/n) " yn
    case $yn in
        [Yy]* ) USE_CONDA=true; break;;
        [Nn]* ) USE_CONDA=false; break;;
        "" ) USE_CONDA=true; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

# If the user wants to use a conda environment, ask for its name. Then, check if
# it exists, otherwise create it (asking for the python version).
if [ "$USE_CONDA" = true ]; then
    # Check if conda is installed
    if ! command -v conda &> /dev/null
    then
        echo -e "${RED}Conda is not installed, please install it before continuing.${NORMAL}"
        exit
    fi

    # Source conda
    source $(conda info --base)/etc/profile.d/conda.sh

    # Ask for the conda environment name
    read -p "Enter the name of the conda environment: " CONDA_ENV

    # Check if the conda environment exists, otherwise create it
    if conda env list | grep -q "$CONDA_ENV"; then
        echo -e "${YELLOW}Conda environment $CONDA_ENV exists, using it.${NORMAL}"
    else
        echo -e "${YELLOW}Conda environment $CONDA_ENV does not exist, creating it.${NORMAL}"
        read -p "Enter the python version (e.g. 3.11): " python_version
        conda create -n $CONDA_ENV python=$python_version -y
    fi

    # Activate the conda environment
    conda deactivate
    conda activate $CONDA_ENV
fi

#1# Required dependencies

# Make sure that the required dependencies are installed
pip install sphinx sphinx-material > /dev/null & spinner $! "${YELLOW}Making sure that the required Python dependencies are installed...${NORMAL}"

#1# Documentation setup
echo -e "${GREEN}Everything is ready to start creating the documentation.${NORMAL}\n"

#2# Ask for the documentation path

default_docs_path="$(realpath .)/docs"
read -p "Enter the documentation absolute path (default: ./docs): " docs_path
if [ -z "$docs_path" ]; then
    docs_path="$default_docs_path"
fi

#2# Ask for the project name

default_project_name="DRIVEWISE"
read -p "Enter the project name (default: $default_project_name): " project_name
if [ -z "$project_name" ]; then
    project_name=$default_project_name
fi

#2# Ask for the author name(s)

default_author_name="DRIVEWISE"
read -p "Enter the author name(s) (default: $default_author_name): " author_name
if [ -z "$author_name" ]; then
    author_name=$default_author_name
fi

#2# Ask for the version

default_version="0.1"
read -p "Enter the version (default: $default_version): " version
if [ -z "$version" ]; then
    version=$default_version
fi

#2# Ask for the release

default_release="0.1"
read -p "Enter the release (default: $default_release): " release
if [ -z "$release" ]; then
    release=$default_release
fi

#2# Ask for the language

default_language="en"
read -p "Enter the language (default: $default_language): " language
if [ -z "$language" ]; then
    language=$default_language
fi

#2# Ask if there will be a Python package to document

DOC_PYTHON_PACKAGE=false
while true; do
    read -p "Will there be a Python package to document? ([y]/n) " yn
    case $yn in
        [Yy]* ) DOC_PYTHON_PACKAGE=true; break;;
        [Nn]* ) DOC_PYTHON_PACKAGE=false; break;;
        "" ) DOC_PYTHON_PACKAGE=true; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

# If there will be a Python package to document, ask for its absolute path
if [ "$DOC_PYTHON_PACKAGE" = true ]; then
    default_python_package_path="$(realpath .)/src"
    read -p "Enter the Python package absolute path (default: ./src): " python_package_path
    if [ -z "$python_package_path" ]; then
        python_package_path="$default_python_package_path"
    fi
fi

#2# Ask if there will be a C++ package to document

DOC_CPP_PACKAGE=false
while true; do
    read -p "Will there be a C++ package to document? ([y]/n) " yn
    case $yn in
        [Yy]* ) DOC_CPP_PACKAGE=true; break;;
        [Nn]* ) DOC_CPP_PACKAGE=false; break;;
        "" ) DOC_CPP_PACKAGE=true; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

# If there will be a C++ package to document, ask for its absolute path
if [ "$DOC_CPP_PACKAGE" = true ]; then
    # Check if doxygen is installed
    if ! command -v doxygen &> /dev/null
    then
        echo -e "${RED}Doxygen is not installed, please install it before continuing.${NORMAL}"
        exit
    fi

    # Make sure breathe is installed
    pip install breathe > /dev/null & spinner $! "${YELLOW}Making sure that the required Python dependencies are installed...${NORMAL}"

    # Ask for the doxygen path
    default_cpp_package_path="$(realpath .)/include"
    read -p "Enter the C++ package absolute path (default: ./include): " cpp_package_path
    if [ -z "$cpp_package_path" ]; then
        cpp_package_path="$default_cpp_package_path"
    fi
    doxygen_path="$docs_path"/source/_doxygen
fi

#1# Create the documentation

# Run sphinx-quickstart in quiet mode
sphinx-quickstart -q -p "$project_name" -a "$author_name" -v $version -r $release -l $language --sep --extensions=sphinx.ext.autodoc,breathe "$docs_path" > /dev/null

#2# Override the conf.py file

# Get the absolute path of the conf.py file
conf_file="$docs_path"/source/conf.py

# Use the Python file override_conf_py.py to override the conf.py file, whilst
# checking if it succeeds
python override_files/override_conf_py.py "$conf_file" "$DOC_PYTHON_PACKAGE" "$python_package_path" "$DOC_CPP_PACKAGE" "$doxygen_path" "$project_name"
if [ $? -eq 1 ]; then
    echo -e "${RED}Failed to override the conf.py file.${NORMAL}"
    exit
fi

#2# Override the index.rst file (only if there is a Python package to document)

# Get the absolute path of the index.rst file
index_file="$docs_path"/source/index.rst

# Use the Python file override_index_rst.py to override the index.rst file,
# whilst checking if it succeeds
python override_files/override_index_rst.py "$index_file" "$DOC_PYTHON_PACKAGE" "$DOC_CPP_PACKAGE"
if [ $? -eq 1 ]; then
    echo -e "${RED}Failed to override the index.rst file.${NORMAL}"
    exit
fi

#2# Setup doxygen (only if there is a C++ package to document)
if [ "$DOC_CPP_PACKAGE" = true ]; then
    # Create the Doxyfile
    doxygen -g "$docs_path"/source/Doxyfile > /dev/null

    # Use the Python file override_doxyfile.py to override the Doxyfile, whilst
    # checking if it succeeds
    python override_files/override_doxyfile.py "$docs_path"/source/Doxyfile $DOC_CPP_PACKAGE "$project_name" "$version" "$doxygen_path" "$cpp_package_path"
    if [ $? -eq 1 ]; then
        echo -e "${RED}Failed to override the Doxyfile.${NORMAL}"
        exit
    fi
fi

echo -e "\n${GREEN}Created the documentation.${NORMAL}"

#1# Create custom building scripts

# Build HTML script
build_script="$docs_path"/build_html.sh
echo "#!/bin/bash" >> "$build_script"
echo "" >> "$build_script"
echo "# This script is used to build the Sphinx documentation." >> "$build_script"
echo "" >> "$build_script"
echo "# Welcome the user" >> "$build_script"
echo "echo -e \"${GREEN}Welcome to the DRIVEWISE documentation build script!${NORMAL}\n\"" >> "$build_script"
echo "" >> "$build_script"
if [ "$USE_CONDA" = true ]; then
    echo "# Activate the conda environment" >> "$build_script"
    echo "source \$(conda info --base)/etc/profile.d/conda.sh" >> "$build_script"
    echo "conda activate $CONDA_ENV" >> "$build_script"
    echo "" >> "$build_script"
fi
echo "# Build the documentation" >> "$build_script"
source_abs_path="$(realpath "$docs_path"/source)"
build_abs_path="$(realpath "$docs_path"/build)"
if [ "$DOC_PYTHON_PACKAGE" = true ]; then
    package_abs_path="$(realpath "$python_package_path")"
    echo "# Delete the _autodoc folder if it exists" >> "$build_script"
    echo "if [ -d \"$source_abs_path/_autodoc\" ]; then" >> "$build_script"
    echo "   rm -rf \"$source_abs_path/_autodoc\"" >> "$build_script"
    echo "fi" >> "$build_script"
    echo "sphinx-apidoc -q -o \"$source_abs_path/_autodoc\" \"$package_abs_path\" > /dev/null"  >> "$build_script"
fi
if [ "$DOC_CPP_PACKAGE" = true ]; then
    echo "# Delete the _doxygen folder if it exists" >> "$build_script"
    echo "if [ -d \"$source_abs_path/_doxygen\" ]; then" >> "$build_script"
    echo "   rm -rf \"$source_abs_path/_doxygen\"" >> "$build_script"
    echo "fi" >> "$build_script"
    echo "doxygen \"$source_abs_path/Doxyfile\" > /dev/null" >> "$build_script"
    echo "# Delete the _breathe folder if it exists" >> "$build_script"
    echo "if [ -d \"$source_abs_path/_breathe\" ]; then" >> "$build_script"
    echo "   rm -rf \"$source_abs_path/_breathe\"" >> "$build_script"
    echo "fi" >> "$build_script"
    echo "breathe-apidoc -o \"$source_abs_path/_breathe\" -m -p \"$project_name\" -q \"$doxygen_path/xml\" > /dev/null" >> "$build_script"
fi
echo "sphinx-build -b html \"$source_abs_path\" \"$build_abs_path/html\" > /dev/null" >> "$build_script"
echo "" >> "$build_script"
if [ "$USE_CONDA" = true ]; then
    echo "# Deactivate the conda environment" >> "$build_script"
    echo "conda deactivate" >> "$build_script"
    echo "" >> "$build_script"
fi
echo "echo -e \"\n${GREEN}Built the HTML documentation in "$docs_path"/build/html.${NORMAL}\n\"" >> "$build_script"
echo "echo \"${GREEN}Enjoy!${NORMAL}\"" >> "$build_script"
