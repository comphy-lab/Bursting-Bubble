# Young–Laplace initial conditions

Public report on the Stage-1 cavity generators in
`simulationCases/initialConditions/`.

- Source: [`young-laplace-initial-conditions.tex`](young-laplace-initial-conditions.tex)
- PDF: [`young-laplace-initial-conditions.pdf`](young-laplace-initial-conditions.pdf)

```bash
cd docs/initial-conditions
python3 figures/make_figures.py
pdflatex young-laplace-initial-conditions.tex
bibtex young-laplace-initial-conditions
pdflatex young-laplace-initial-conditions.tex
pdflatex young-laplace-initial-conditions.tex
```

The opening-angle figure solves a Bond continuation from $10^{-4}$ to
$900$ and is the slow step. Use the repository `initialConditions`
virtualenv if the system Python has no SciPy.
