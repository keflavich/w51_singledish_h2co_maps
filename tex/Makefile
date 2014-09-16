# latex Makefile
texpath=/usr/texbin/
PDFLATEX=${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
LATEX=${PDFLATEX}
BIBTEX=bibtex
DVIPS=dvips
PS2PDF=ps2pdf

all: aa_normal aa_referee w51_maps 

w51_maps: aa_normal latex

latex:
	cat preface.tex >> preface_aa.tex
	${LATEX} w51_maps.tex
	${BIBTEX} w51_maps
	${LATEX} w51_maps.tex
	${BIBTEX} w51_maps
	${LATEX} w51_maps.tex
	gs -dSAFER -dBATCH -dNOPAUSE  -sDEVICE=pdfwrite -sOutputFile=w51_maps_compressed.pdf w51_maps.pdf

latex_referee: aa_referee latex

aa_normal:
	@rm -f w51_maps*.aux w51_maps*.bbl w51_maps*.blg w51_maps*.dvi w51_maps*.log w51_maps*.lot w51_maps*.lof
	#document="\documentclass[referee]{aa} \input{w51_maps.tex}"
	echo "\documentclass{aa}" > preface_aa.tex

aa_referee: 
	@rm -f w51_maps*.aux w51_maps*.bbl w51_maps*.blg w51_maps*.dvi w51_maps*.log w51_maps*.lot w51_maps*.lof
	#document="\documentclass[referee]{aa} \input{w51_maps.tex}"
	echo "\documentclass[referee]{aa}" > preface_aa.tex
