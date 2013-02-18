PROJ=tutorial

all: ${PROJ}.pdf

${PROJ}.pdf: tutorial.tex Makefile 
	-pdflatex ${PROJ}
	-bibtex ${PROJ}
	-pdflatex ${PROJ}
	-pdflatex ${PROJ} | grep --color='auto' 'Warning\|Error'


tutorial.tex: tutorial.Rnw
	R CMD Sweave tutorial.Rnw

clean:
	rm -rf *.bbl *.blg *.log *.aux *~ *-*.pdf *.toc *.tex
