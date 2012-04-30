all: pse.pdf

pse.pdf: pse.tex chalom.bib Makefile capitulo1.tex
	-pdflatex pse
	-bibtex pse
	-pdflatex pse
	pdflatex pse

pse.tex: pse.Rnw
	R CMD Sweave pse.Rnw

capitulo1.tex: capitulo1.Rnw
	R CMD Sweave capitulo1.Rnw

clean:
	rm -rf *.dvi *.bbl *.blg *.log *.aux *~ *.pdf

