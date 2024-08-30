# You can set these variables from the command line.
SPHINX_OPTS     =
BUILD_DIR       = build
DOCS_SOURCE_DIR = .
DOXYFILE        = Doxyfile_custom
DOXYFILE2       = Doxyfile_custom2

.PHONY: help clean html doxygen doxysphinx sphinx

help:
	@echo "Please use make <target> where <target> is one of"
	@echo "  sphinx      to run sphinx."
	@echo "  doxygen     to run doxygen."
	@echo "  doxysphinx  to run doxysphinx."
	@echo "  html        to run everything in correct order to generate the documentation."
	@echo "  clean       to clean up everything"

clean:
	doxysphinx --verbosity=DEBUG clean $(DOCS_SOURCE_DIR) $(BUILD_DIR)/html $(DOXYFILE)
	doxysphinx --verbosity=DEBUG clean $(DOCS_SOURCE_DIR) $(BUILD_DIR)/html $(DOXYFILE2)
# cannot use rm here because it's not portable (win/linux)
	@rm -rf doxygen/*
	@rm -rf doxygen2/*
	@rm -rf build/*

# sphinx
sphinx:
	sphinx-build -M html "$(DOCS_SOURCE_DIR)" "$(BUILD_DIR)" '--keep-going' '-j' 'auto' $(SPHINX_OPTS)

# some aliases
doxygen:
	@doxygen $(DOXYFILE)
	@doxygen $(DOXYFILE2)
	@ruby filter.rb

doxysphinx:
	doxysphinx --verbosity=DEBUG build $(DOCS_SOURCE_DIR) $(BUILD_DIR)/html $(DOXYFILE)
	doxysphinx --verbosity=DEBUG build $(DOCS_SOURCE_DIR) $(BUILD_DIR)/html $(DOXYFILE2)

html: doxygen doxysphinx sphinx
	@echo ""
	@echo "Build finished. The Documentation is in $(BUILD_DIR)/html."
