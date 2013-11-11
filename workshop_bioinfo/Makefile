all: poster.pdf

poster.pdf: poster.tex chalom.bib Makefile beamerthemeDreuw.sty
	pdflatex poster
	bibtex poster
	pdflatex poster
	pdflatex poster

poster.tex: poster.Rnw
	./mkcap.sh poster

clean:
	rm poster.blg poster.out poster.pdf poster.log poster.nav poster.toc poster.aux poster.snm poster.bbl poster.tex
