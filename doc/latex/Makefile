all: pdf
	cp EMAN2.pdf ../
	cp EMAN2.ps ../

html: index.html

ps: EMAN2.ps

pdf: EMAN2.pdf

dvi: EMAN2.dvi

index.html: EMAN2.tex *.tex
	latex2html EMAN2.tex

EMAN2.ps: EMAN2.dvi
	dvips -o EMAN2.ps EMAN2.dvi

EMAN2.pdf: EMAN2.ps
	ps2pdf EMAN2.ps EMAN2.pdf

EMAN2.dvi: EMAN2.tex *.tex
	latex EMAN2
	makeindex EMAN2.idx
	latex EMAN2

clean:
	rm -f EMAN2.ps EMAN2.pdf *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out
min: clean
	rm -f *~