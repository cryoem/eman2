all: pdf

pdf: sparx.pdf
ps: sparx.ps
dvi: sparx.dvi

sparx.pdf: sparx.ps
	ps2pdf sparx.ps sparx.pdf

sparx.ps: sparx.dvi
	dvips -o sparx.ps sparx.dvi

sparx.dvi: *.tex chapters/*
	latex sparx.tex
	makeindex sparx.idx
	latex sparx.tex
	makeindex sparx.idx
	latex sparx.tex

force_dvi:
	latex sparx.tex
	makeindex sparx.idx
	latex sparx.tex
	makeindex sparx.idx
	latex sparx.tex

clean:
	rm -f *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.lot
dist: clean
	rm -f *~