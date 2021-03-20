github:
	@cd sphinx; make clean; make html; cd ..
	@make install

install:
	@cd sphinx; rm -rf ../../docs/*; cp -a _build/html/. ../../docs; cd ..
