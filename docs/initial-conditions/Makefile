PYTHON ?= ../../simulationCases/initialConditions/.venv/bin/python
TEX = young-laplace-initial-conditions

.PHONY: all figures pdf clean

all: pdf

figures:
	$(PYTHON) figures/make_figures.py

pdf: figures
	pdflatex -interaction=nonstopmode $(TEX).tex
	bibtex $(TEX)
	pdflatex -interaction=nonstopmode $(TEX).tex
	pdflatex -interaction=nonstopmode $(TEX).tex

clean:
	rm -f $(TEX).aux $(TEX).log $(TEX).out $(TEX).toc $(TEX).bbl $(TEX).blg
