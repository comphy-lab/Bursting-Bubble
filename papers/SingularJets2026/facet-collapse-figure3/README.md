# Figure 3 — self-similar interface collapse (Bo=0, Oh=0.03)

Data + script to reproduce the interface-profile collapse figure (Duchemin/
Cattaneo-style: pre/post-inception, raw + rescaled).

- `fig3_collapse.py` — self-contained (numpy + matplotlib). Reads `facets/`,
  chains the raw facet segments into continuous solid polylines, keeps the
  largest connected component per frame, fits the inception time t0, rescales
  about the virtual origins t0^- and t0^+, and writes
  `fig3_Oh0.03_collapse.{png,pdf}`.
- `fit_t0_collapse.py` — computes the virtual origins t0^- (before inception)
  and t0^+ (after) that best collapse the profiles for a given alpha, and writes
  `t0_collapse.json` (needs scipy and `compare_collapse_r1.py`).
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

Lengths in R0, time in t_ic (raw sim units). With OUR cone exponent
**alpha = 0.629** (Oh=0.03, beta=38.4 deg; NOT the inertio-capillary 2/3):

    x = (r/R0) / (tau/t_ic)^alpha ,   y = (z - z_base)/R0 / (tau/t_ic)^alpha

with tau = t0^- - t before inception and tau = t - t0^+ after inception.

- **Inception time t0.** From the post-inception jet radius, r_j^(1/alpha)
  linear in t over 0.01 < r_j < 0.06: t0 = 0.494430. It selects the frames and
  sets the colour scale |t - t0|.
- **Frames.** After inception 1.5e-4 <= t - t0 <= 1e-2 (79 frames; the data
  reach 6.3e-3). Before inception 1e-3 <= t0 - t <= 2e-2 (21 frames): closer to
  inception the neck pinches off the entrapped bubble faster than the cavity
  collapses, and earlier cavities are not yet self-similar.
- **Virtual origins.** `fit_t0_collapse.py` scans t0^± and minimises the pairwise
  symmetric Chamfer RMS between the rescaled profiles inside the plotted window
  (|x| <= 4, with -0.3 <= y <= 2.5 before and -0.5 <= y <= 8 after inception),
  measured in units of a 420 x 464 panel. For alpha = 0.629:
  t0^- = 0.494191 and t0^+ = 0.494355. For alpha = 2/3 the same procedure gives
  t0^- = 0.494391 and t0^+ = 0.494270, and the profiles collapse comparably, so
  the interface shape alone does not separate the two exponents.
- **Window.** The rescaled panels show |x| <= 4, the jet and the cavity bottom.
  Further out, the outer cavity flank after inception collapses with 2/3 but not
  with 0.629 for any t0^+.

The post-inception jet flanks collapse onto a single self-similar cone. To use
the L15/focus-15 case (5004) instead, drop its main-facets in `facets/` with the
same naming and re-run.
