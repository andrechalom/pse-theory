PROJ=pse

all: ${PROJ}.pdf

${PROJ}.pdf: ${PROJ}.tex chalom.bib Makefile revlit.tex abstract.tex intro.tex sampling.tex qualanal.tex quantanal.tex ack.tex leslie.tex appendix.tex
	-pdflatex ${PROJ}
	-bibtex ${PROJ}
	-pdflatex ${PROJ}
	-pdflatex ${PROJ} | grep --color='auto' 'Warning\|Error'
	evince ${PROJ}.pdf > /dev/null &

leslie.tex: leslie.Rnw R/Independent.Rdata R/Dependent.Rdata R/pse.R
	./mkcap.sh leslie

abstract.tex: abstract.Rnw
	./mkcap.sh abstract

ack.tex: ack.Rnw
	./mkcap.sh ack

revlit.tex: revlit.Rnw
	./mkcap.sh revlit

intro.tex: intro.Rnw
	./mkcap.sh intro

qualanal.tex: qualanal.Rnw
	./mkcap.sh qualanal

appendix.tex: appendix.Rnw
	./mkcap.sh appendix

quantanal.tex: quantanal.Rnw
	./mkcap.sh quantanal

sampling.tex: sampling.Rnw
	./mkcap.sh sampling

pse.tex: pse.Rnw
	R CMD Sweave pse.Rnw

clean:
	rm -rf *.dvi *.bbl *.blg *.log *.aux *~ *.pdf *.toc *.tex
