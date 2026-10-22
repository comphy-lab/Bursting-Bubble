.PHONY: figures figures-main figures-supp clean

PYTHON ?= python3
NO_TEX ?= 0
FIGURE_FLAGS :=
ifeq ($(NO_TEX),1)
FIGURE_FLAGS := --no-tex
endif

figures: figures-main figures-supp

figures-main:
	$(PYTHON) figure-scripts/make_figures_main.py $(FIGURE_FLAGS)

figures-supp:
	$(PYTHON) figure-scripts/make_figures_supp.py

clean:
	rm -rf output
