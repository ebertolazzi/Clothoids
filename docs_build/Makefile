github:
	@rm -rf xml
	@rm -rf sphinx/_build
	@rm -rf sphinx/api
	@cd sphinx; make clean; make html; cd ..
	@make install

install:
	@cd sphinx; rm -rf ../../docs/*; cp -a _build/html/. ../../docs; cd ..
