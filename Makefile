all: pse.pdf

pse.pdf: pse.tex chalom.bib Makefile
	-pdflatex pse
	-bibtex pse
	pdflatex pse

clean:
	rm -rf *.dvi *.bbl *.blg *.log *.aux *~ *.pdf

