PROJ=pse

all: ${PROJ}.pdf

${PROJ}.pdf: ${PROJ}.tex chalom.bib Makefile revlit.tex abstract.tex intro.tex sampling.tex quantanal.tex ack.tex leslie.tex
	-pdflatex ${PROJ}
	-bibtex ${PROJ}
	-pdflatex ${PROJ}
	-pdflatex ${PROJ} | grep --color='auto' 'Warning\|Error'

integr: integr.tex
	pdflatex integr

caveat:
	pdflatex new/caveat

leslie.tex: leslie.Rnw R/Independent.Rdata R/Dependent.Rdata 

clean:
	rm -rf *.bbl *.blg *.log *.aux *~ *-*.pdf *.toc *.tex

%.tex: %.Rnw
	R --vanilla <<< "library(utils); Sweave(\"$<\", encoding=\"utf8\")"
