all : ppgtk.pdf ppgtk-notes.pdf extending-ppgtk.pdf

.PHONY: copy

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

extending-ppgtk.pdf : extending-ppgtk.tex
	pdflatex extending-ppgtk
	pdflatex extending-ppgtk
	#bibtex extending-ppgtk
	#pdflatex extending-ppgtk
	@cp extending-ppgtk.pdf ../

copy :
	@cp ppgtk.pdf ../manual.pdf
	@cp ppgtk-notes.pdf ../
	@cp extending-ppgtk.pdf ../
