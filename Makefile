all: pse.pdf


pse.pdf: pse.tex 


hjkhjdkfsldchalom.bib Makefile
	-pdflatex pse
	-bibtex pse
	pdflatex pse

clean:
	rm -rf *.dvi *.bbl *.blg *.log *.aux *~ 

