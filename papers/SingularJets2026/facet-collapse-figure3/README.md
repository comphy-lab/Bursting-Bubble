# Figure 3 — self-similar interface collapse (Bo=0, Oh=0.03)

Data + script to reproduce the interface-profile collapse figure (Duchemin/
Cattaneo-style: pre/post-inception, raw + rescaled).

- `fig3_collapse.py` — self-contained (numpy + matplotlib). Reads `facets/`,
  chains the raw facet segments into continuous solid polylines, keeps the
  largest connected component per frame, fits t0, and writes
  `fig3_Oh0.03_collapse.{png,pdf}`.
- Run from the parent folder: `bash reproduce_facets.sh` (or directly
  `python3 fig3_collapse.py`).
- `extract_full.py` — the snellius-side extractor that produced `facets/`
  (runs `getFacet` on each raw dump, filters to r<0.6, rebuilds the index from
  the run log). Kept for provenance; not needed to reproduce the figure.
- `facets/` — FULL raw interface facets (`getFacet`) from case 5003 (L15,
  focus 14), 44 pre + 81 post frames:
  - `facetpremain_<t>.txt` — pre-inception (cavity collapse)
  - `facetmain_<t>.txt`    — post-inception (jet)
  - each file: `z r` interface segments (blank line between segments)
  - `index_pre.txt`, `index.txt` — per time t: `t  r_j  z_base` (from the run
    log: r_j = getBase jet radius, z_base = tag-based main-body base). `z_base`
    is the shift reference (a raw min(z) would be polluted by the entrained
    bubble).
  - The entrapped bubble and shed droplets are NOT in the drawing: they form
    their own small connected components, and `fig3_collapse.py` keeps only the
    largest component per frame. (Earlier versions used `getFacetMain`, whose
    main-body tag dropped ~75% of the neck facets right after inception and left
    a segmented neck — hence the switch to the full interface.)

## The rescaling

Lengths in R0, time in t_ic (raw sim units). tau = |t - t0|. With OUR cone
exponent **alpha = 0.629** (Oh=0.03, beta=38.4 deg; NOT the inertio-capillary
2/3):

    x = (r/R0) / (tau/t_ic)^alpha ,   y = (z - z_base)/R0 / (tau/t_ic)^alpha

The post-inception jet flanks collapse onto a single self-similar cone. To use
the L15/focus-15 case (5004) instead, drop its main-facets in `facets/` with the
same naming and re-run.
