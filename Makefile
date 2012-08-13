PROJ=pse

all: ${PROJ}.pdf

${PROJ}.pdf: ${PROJ}.tex chalom.bib Makefile capitulo1.tex revlit.tex
	-pdflatex ${PROJ}
	-bibtex ${PROJ}
	-pdflatex ${PROJ}
	pdflatex ${PROJ} | grep --color='auto' 'Warning\|Error'
	evince ${PROJ}.pdf > /dev/null &

capitulo1.tex: capitulo1.Rnw corcorr.so
	./mkcap.sh

pse.tex: pse.Rnw
	R CMD Sweave pse.Rnw

clean:
	rm -rf *.dvi *.bbl *.blg *.log *.aux *~ *.pdf capitulo1.tex ${PROJ}.tex *.toc
