# latex Makefile
texpath=/usr/texbin/
PDFLATEX=${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
LATEX=${PDFLATEX}
BIBTEX=bibtex
DVIPS=dvips
PS2PDF=ps2pdf

all: w51_maps

w51_maps: 
	@rm -f w51_maps*.aux w51_maps*.bbl w51_maps*.blg w51_maps*.dvi w51_maps*.log w51_maps*.lot w51_maps*.lof
	${LATEX} w51_maps.tex
	${BIBTEX} w51_maps
	${LATEX} w51_maps.tex
	${BIBTEX} w51_maps
	${LATEX} w51_maps.tex
