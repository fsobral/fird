all: doc

doc: doc.tex
	export BIBINPUTS=~/Documents/
	pdflatex doc.tex
	bibtex doc
	pdflatex doc.tex
	pdflatex doc.tex

clean:
	rm -f *~ *.aux *.log *.pdf *.toc *.bbl *.blg *.out

.PHONY: all clean doc
