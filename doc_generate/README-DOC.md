# README

## Setup the code and environment

Install miniconda_ or a similar tool (e.g.  Anaconda_) and create a conda
environment for the book::

```{bash}
conda env create -f doxygen-doc-env.yml
```

- [miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [Anaconda](https://www.anaconda.com/products/individual)

Activate the conda environment::

```{bash}
conda activate doxygen-doc
```

```
doxysphinx build . out html
sphinx-build -b html . ../out
```


## Build and view the website

To build the website run::

```{bash}
make html
```

When complete, the website is then viewable in your browser::

```
<yourbrowser> _build/html/index.html
```

You can also run sphinx-autobuild (updates while while you edit) with::

```
make autobuild
```

## remove conda env

```{bash}
conda deactivate
conda remove --name doxygen-doc --all
```
