all:
	pdflatex prob.tex
	bibtex prob.aux
	pdflatex prob.tex
	pdflatex prob.tex
	apvlv prob.pdf

clean:
	rm -rf *.bbl *.blg *.log *.aux *~ *-*.pdf *.toc 
