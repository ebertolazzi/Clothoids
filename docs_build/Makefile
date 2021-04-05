INDEX = ./sphinx/_build/html/genindex.html

build:
	make clean
	@cd sphinx_cpp;    make html; cd ..
	@cd sphinx_matlab; make html; cd ..
	@cd sphinx;        make html; cd ..
	@cp -rf sphinx_cpp/_build/html    sphinx/_build/html/api-cpp
	@cd sphinx_matlab; ruby filter.rb; cd ..
	@cp -rf sphinx_matlab/_build/html sphinx/_build/html/api-matlab
	@make install

install:
	@cd sphinx; rm -rf ../../docs/*; cp -a _build/html/. ../../docs; cd ..

clean:
	@rm -rf xml
	@rm -rf xml-cpp
	@rm -rf xml-matlab
	@rm -rf sphinx*/api-*
	@rm -rf sphinx*/_build
	@cd sphinx_cpp;    make clean; cd ..
	@cd sphinx_matlab; make clean; cd ..
	@cd sphinx;        make clean; cd ..
