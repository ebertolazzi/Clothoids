FIGSASY = $(shell ls *.asy)

all:
	make LC_NUMERIC=C -f Makefile_asy

clean:
	rm -f $(patsubst %.asy, %.pdf, $(FIGSASY))
	rm -f $(patsubst %.asy, %.eps, $(FIGSASY))