FROMCPP    = ./sphinx_cpp/api-cpp/
FROMMATLAB = ./sphinx_matlab/api-matlab/
TOCPP      = ./sphinx/api-cpp/
TOMATLAB   = ./sphinx/api-matlab/

build:
	make clean
	@cd sphinx; make html; cd ..
	@cd sphinx_matlab; make html; ruby filter.rb; cd ..
	@make install

install:
	@mv sphinx/_build/html/api-matlab/root.html sphinx/root.html
	@rm -rf sphinx/_build/html/api-matlab
	@cp -rf sphinx_matlab/_build/html/api-matlab sphinx/_build/html/api-matlab
	@cp -f sphinx/root.html sphinx/_build/html/api-matlab/root.html
	@cp -rf sphinx_matlab/_build/html/genindex.html sphinx/_build/html/genindex2.html
	@rm -rf ../docs/*; cp -a sphinx/_build/html/. ../docs

clean:
	@rm -rf xml
	@rm -rf xml-cpp
	@rm -rf xml-matlab
	@rm -rf sphinx*/api-*
	@rm -rf sphinx*/_build
	@cd sphinx_matlab; make clean; cd ..
	@cd sphinx;        make clean; cd ..
