PROJ=tutorial

all: ${PROJ}.pdf

${PROJ}.pdf: tutorial.tex Makefile 
	-pdflatex ${PROJ}
	-bibtex ${PROJ}
	-pdflatex ${PROJ}
	-pdflatex ${PROJ} | grep --color='auto' 'Warning\|Error'

zip: tutorial.pdf pse.R corcorr.c
	zip -r tutorial.zip tutorial.pdf pse.R corcorr.c

tutorial.tex: tutorial.Rnw
	R CMD Sweave tutorial.Rnw

clean:
	rm -rf *.bbl *.blg *.log *.aux *~ *-*.pdf *.toc *.tex *.o *.so *.zip
