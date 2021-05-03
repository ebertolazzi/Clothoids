build:
	make clean
	@cd sphinx; make html; ruby filter.rb; cd ..
	@make install

install:
	@rm -rf ../docs/*; cp -a sphinx/_build/html/. ../docs

clean:
	@rm -rf xml
	@rm -rf xml-cpp
	@rm -rf xml-matlab
	@rm -rf sphinx*/api-*
	@rm -rf sphinx*/_build
	@cd sphinx; make clean; cd ..
