PROJ=pse

all: ${PROJ}.pdf

${PROJ}.pdf: ${PROJ}.tex chalom.bib Makefile capitulo1.tex
	-pdflatex ${PROJ}
	-bibtex ${PROJ}
	-pdflatex ${PROJ}
	pdflatex ${PROJ}

%.tex:
	R CMD Sweave $*.Rnw

clean:
	rm -rf *.dvi *.bbl *.blg *.log *.aux *~ *.pdf capitulo1.tex ${PROJ}.tex

