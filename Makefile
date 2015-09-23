PROJ=pse

all: pse.pdf arxiv1.pdf arxiv2.pdf

pse.pdf: pse.tex chalom.bib Makefile revlit.tex intro.tex sampling.tex quantanal.tex leslie.tex preamble.tex prob.tex integr.tex mini.tex tribolium.tex pse_tutorial.tex multiple.tex plue.tex plue_tutorial.tex

arxiv1.pdf: intro.tex abstract.tex sampling.tex quantanal.tex leslie.tex tribolium.tex revlit.tex

leslie.tex: leslie.Rnw R/leslie.Rdata R/leslie.csv

tribolium.tex: tribolium.Rnw R/Triboliumpeq.Rdata R/Triboliumlge.Rdata

%.pdf: %.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*
	pdflatex $*

clean:
	rm -rf *.bbl *.blg *.log *.aux *~ *-*.pdf *.toc intro.tex leslie.tex quantanal.tex sampling.tex mini.tex tribolium.tex pse_tutorial.tex mydata.csv multiple.tex integr.tex pse.tex Rplots.pdf plue.tex plue_tutorial.tex pse.pdf arxiv2.tex

%.tex: %.Rnw
	R --vanilla <<< "library(utils); Sweave(\"$<\", encoding=\"utf8\")"
