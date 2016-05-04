all : ppgtk.pdf ppgtk-notes.pdf

ppgtk.pdf : ppgtk.tex #ppgtk.bib
	pdflatex ppgtk
	pdflatex ppgtk
	#bibtex ppgtk
	#pdflatex ppgtk
	@cp ppgtk.pdf ../manual.pdf

ppgtk-notes.pdf : ppgtk-notes.tex #ppgtk-notes.bibtex
	pdflatex ppgtk-notes
	pdflatex ppgtk-notes
	#bibtex ppgtk-notes
	#pdflatex ppgtk-notes
	@cp ppgtk-notes.pdf ../
