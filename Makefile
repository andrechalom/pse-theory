PROJ=pse

all: ${PROJ}.pdf

${PROJ}.pdf: ${PROJ}.tex chalom.bib Makefile revlit.tex intro.tex sampling.tex quantanal.tex leslie.tex preamble.tex prob.tex caveat.tex  integr.tex mini.tex tribolium.tex pse_tutorial.tex multiple.tex twolife.tex
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
	rm -rf *.bbl *.blg *.log *.aux *~ *-*.pdf *.toc intro.tex leslie.tex quantanal.tex sampling.tex mini.tex tribolium.tex pse_tutorial.tex mydata.csv multiple.tex twolife.tex integr.tex

%.tex: %.Rnw
	R --vanilla <<< "library(utils); Sweave(\"$<\", encoding=\"utf8\")"
