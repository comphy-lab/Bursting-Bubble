# SingularJets2026 bursting-bubble data and reproduction capsule

This standalone capsule contains the solver logs, exact run provenance,
regular-grid fields, interface segments and locked Python environment needed to
reproduce the manuscript figures from a clean checkout.  The default workflow
is offline and does not require SSH, Basilisk, `qcc`, `getFacet` or `pdftoppm`.

## Layout — classified by Ohnesorge number

- `plotJetMetricsTheory.py` — the (self-contained) plotting tool, shared by all
  cases. Also on the CoMPhy Lab repo, merged into `main` (PR #6).
- **`figure-scripts/`** — scripts and small CSV inputs for manuscript figures
  that are not produced by the shared jet-metrics diagnostic script: the Fig. 1
  schematic (`make_sketch_paper_worthington.py`), the approved manuscript Fig. 2
  (`make_fig2_flux_scalings.py`), and the End Matter natural-numbered Fig. 4
  diagnostic (`make_nu_vs_beta.py`, `beta_alpha_Oh.csv`).
- **`data-fig2a/`** — the four committed regular-grid field arrays and interface
  segment sets used by Fig. 2(a), plus checksummed raw-snapshot metadata.  These
  are the default inputs for both the panel renderer and full-figure compositor.
- **`data-Oh-0.03/`** — the case we focus on (Bo=0, Oh=0.03). L13, L14, and
  three L15 runs (pre-inception focus 13, 14, and 15). `bash reproduce-Oh0.03.sh`
  → `data-Oh-0.03/gridconv3_Oh0.03_L13-L14-L15.{png,pdf}`.
- **`facet-collapse-figure3/`** — Figure 3: the self-similar interface-collapse figure (pre/post inception, raw + rescaled by |t-t0|^alpha). `bash reproduce_facets.sh` -> `facet-collapse-figure3/fig3_Oh0.03_collapse.{png,pdf}`. See that folder's README.
- **`data-Oh-0.02297/`** — Bo=0, Oh=0.02297. L13 + L14. `bash reproduce-Oh0.02297.sh` →
  `data-Oh-0.02297/gridconv_Oh0.02297_L13-L14.{png,pdf}`.
- **`metadata.json`** — canonical scientific choices, including the Fig. 2 fit
  window and regression prefactors.
- **`cases.csv`** and **`provenance/`** — the eight-case ledger, exact 5008
  restart/job evidence, solver Git bundle and build/runtime environment.
- **`pyproject.toml`** and **`uv.lock`** — the locked figure environment.

Each `data-Oh-*/` folder holds `<case>_L<level>[_focus<f>]_log.txt` (+
`case.params`), one `facet_inception.txt` (the interface `x y` cloud at
inception, for the cone fit), the generated figure (the `reproduce-Oh*.sh` runners live at the folder root).

## Manuscript figure map

- **Fig. 1:** `sketch_paper_worthington.{pdf,png}`; regenerate with
  `python3 figure-scripts/make_sketch_paper_worthington.py`.
- **Fig. 2:** `fig2_flux_scalings.{pdf,png}`; regenerate with
  `python3 figure-scripts/make_fig2_flux_scalings.py`.
  The current manuscript version is `fig2_v2.{pdf,png}`; regenerate with
  `bash reproduce-fig2-v2.sh`.  That runner first invokes the standalone
  `make_fig2a_streamlines.py` panel/subpanel renderer and then the separate
  `make_fig2_v2.py` full-figure compositor.  Both PDF and PNG outputs are
  written atomically from the same Matplotlib figures. Panel (c) samples an interpolated $Q_j(r_j)$
  branch before computing $We_j=Q_j^2/(\pi^2r_j^3)$ to avoid repeated $r_j$
  markers.
  The supporting raw grid-convergence diagnostic is
  `data-Oh-0.03/gridconv3_Oh0.03_L13-L14-L15.{pdf,png}`, regenerated with
  `bash reproduce-Oh0.03.sh`.
- **Fig. 3:** `facet-collapse-figure3/fig3_Oh0.03_collapse.{pdf,png}`;
  regenerate with `bash reproduce_facets.sh`.
- **Fig. 4 / End Matter diagnostic:** `nu_vs_beta.pdf`; regenerate with
  `python3 figure-scripts/make_nu_vs_beta.py`.

### Grid / focus notation

`L<n>` = MAXlevel n (the finest refinement, reached ~when v_jet is highest).
`focus <m>` = the pre-inception cap: the max resolution allowed *just before*
the jet forms (the cavity-focus collapse is capped at level m, then released to
the full MAXlevel once the jet is issued). At Oh=0.03 the converged data is
insensitive to focus 13 vs 14 in the asymptotic window.  Case 5004 was not
simply cancelled: its decisive rerun reached the focus-15 window and then
tripped the default kinetic-energy gate (`keStopMax=1e4`) at `t=0.49505452`.
Case 5008 restarted its snapshot `0.490625`, raised the gate to `1e6`, and ran
until Slurm job 24472158 reached its six-hour limit at `t=0.49776226`.
`cases.csv` records both the earlier manually cancelled 5004 attempt and the
decisive rerun.  The restart checksum is intentionally not recorded.

The 5001 focus13 log is a full-tail L15 run. The 5003 focus14 log included here
captures the asymptotic jet window and reached L15; use it as a focus-sensitivity
check, not as a full-tail video replacement.

## Log format

After two header lines, each row is whitespace-separated:

```
i  dt  t  ke  maxlevel  r_b  z_b  r_base  z_base  q_jet  q_l
```

The plots use, at the jet-base plane:
- `r_base` (col 8) = **r_j**, jet-base radius (robust getBase probe).
- `q_jet`  (col 10) = **∫ v_z r dr**; manuscript Fig. 2 uses
  **Q_j = 2π q_jet** and **q_j = Q_j/(π r_j)**.
- `q_l`    (col 11) = **∫ v_z dr**; retained as a diagnostic, but not used
  for the approved manuscript Fig. 2 line-flux panel.

`Oh` and `Bond` are in each `*_case.params` (Bond only picks the initial shape;
the solver has no gravity term).

## Notes

- `--rmin 0.003` shows the full resolved range; the jet base resolves to
  r_j ≈ 0.005, so the historical 0.04 clip is cosmetic, not a resolution limit.
- Theory lines: cone (this work, α=1/(2−ν(β)), β fit from the inception facet),
  inertio-capillary (α=2/3), PRF 2023 (α=1/2). Fit the cone and
  inertio-capillary prefactors over the r_j→0 window
  `--cone-fit-window 0.005 0.023952`.  This inclusive window is loaded from
  `metadata.json`; the tests pin its exact Level-15/focus-15 prefactors.
- Pass `--no-tex` (or set `SINGULARJETS_NO_TEX=1` for the shell runner) on
  systems without a LaTeX installation.

## Offline Fig. 2 workflow

With `uv` installed:

```bash
uv sync --frozen
SINGULARJETS_NO_TEX=1 bash reproduce-fig2-v2.sh
uv run --offline pytest
```

The committed `data-fig2a/fields/*.npz` files contain `z`, `r`, `f`, `uz`,
`ur` and `speed` on the exact 294 × 190 plotting grid.  The companion
`interfaces/*.npz` files contain every mirrored line segment drawn on the four
frames.  `data-fig2a/manifest.json` records their shapes and checksums.

The raw Basilisk snapshots are preserved separately under the requested
Dropbox folder
`2-Resource-Research/1-Github-Files/Bursting-Bubble/SingularJets2026-data/5003`.
Their byte sizes, SHA-256 values and Dropbox file IDs are in
`data-fig2a/raw-snapshots.json`.  Its public `dl=1` links were downloaded
without an authenticated browser and verified against the archived byte counts
and SHA-256 values.  They are used only by the explicit maintainer
`--refresh-data` path; ordinary reproduction remains offline.

## Provenance, checksums and reuse

- `cases.csv` is the authoritative case/parent/job/terminal-state ledger and
  explicitly encodes `5004 → 5008`.
- `provenance/environment.txt` pins the byte-matched simulation source commit,
  Basilisk tag/commit, modules, compile command and thread layout.
- `provenance/simulation-source.bundle` is a self-contained Git bundle for the
  solver commit; `provenance/raw/5008/` preserves the unmodified scheduler
  streams and restart submission.
- `SHA256SUMS` covers every committed capsule payload except itself.
- Data and figures are CC BY 4.0 under `LICENSE-DATA.md`; scripts retain the
  repository-level GPL-3.0 licence.  Citation metadata is in `CITATION.cff`.
