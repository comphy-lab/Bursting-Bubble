# Singular Jets figure files

Run `make figures-main` for the main-paper figures or `make figures-supp` for
the three supplemental figures. The targets place PDFs in `output/`.

Figure 1 is supplied as the exact vector render; `schematic/make_figure1.py`
also writes a redraw. Figures 2–5 and the supplemental figures are supplied as
the exact vector renders used by the paper. Figure 2 therefore needs no simulation
fields or raw logs. The `comparisons/figure2c_javi_comparison.ipynb` notebook
is retained with its source data; the comparison code was written by Javier
Rodríguez-Rodríguez (Javi).

`make figures-main` writes exact `figure1.pdf`–`figure5.pdf` assets plus
regenerated `figure1-redrawn.pdf`, `figure4a.pdf`, and `figure4b.pdf`.
`make figures-supp` writes exact `figure-s1.pdf`–`figure-s3.pdf` assets.
`schematic/make_figure1.py` draws the schematic. The Figure 4(a) script solves
`P_nu(-cos(beta)) = 0`; Figure 4(b) is rendered from `data/figure4b.csv`.
Install Python with NumPy, SciPy, and Matplotlib. A LaTeX installation is
needed for publication typography. Use `make figures NO_TEX=1` when LaTeX is
unavailable.

The CSV and text inputs under `data/` are the tabular inputs used by the
included figure scripts and comparison notebook.
