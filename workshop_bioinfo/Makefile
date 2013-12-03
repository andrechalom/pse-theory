all: poster.pdf

poster.pdf: poster.tex chalom.bib Makefile beamerthemeDreuw.sty
	pdflatex poster
	bibtex poster
	pdflatex poster
	pdflatex poster

%.tex: %.Rnw
	R --vanilla <<< "library(utils); Sweave(\"$<\", encoding=\"utf8\")"

clean:
	rm poster.blg poster.out poster*.pdf poster.log poster.nav poster.toc poster.aux poster.snm poster.bbl poster.tex

