# Experimental comparison notebook

`comparison_experiments.ipynb` contains Javier Rodríguez's comparison with
the digitised Zeff and Cattaneo experimental data in `data/`. Run the notebook
from this directory. The supplied notebook is preserved byte-for-byte under
`provenance/javi-comparison-2026-09-09/`, alongside input SHA-256 hashes.
The working notebook removes the macOS-only backend selection and clears
stored outputs; equations, parameters and data are unchanged.

For a non-interactive run using the existing capsule environment:

```bash
.venv/bin/python figure-scripts/reproduce_comparison_experiments.py \
  --output /path/to/review-output
```

The runner executes the code cells in order and saves three PDF/PNG plots,
an executed notebook and `light-check.json`. It compares printed numerical
results with the supplied notebook at its saved precision and independently
checks the three similarity roots with a bracketed solver. The original
manual prefactors and critical heights are retained; this is reproduction,
not refitting. It does not verify the source papers or introduce a viscous
regularisation at the critical height.
