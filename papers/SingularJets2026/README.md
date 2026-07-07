# Bursting-bubble jet measurements — data + reproduction (for Javi)

Solver logs + a self-contained python script to reproduce the jet-observable
grid-convergence figures vs the self-similar cone theory. The full snapshot
dumps are large (~5–6 GB per Level-15 case, from the refinement near jet
inception), so only the logs + facet clouds are shipped here — everything the
plots need.

## Layout — classified by Ohnesorge number

- `plotJetMetricsTheory.py` — the (self-contained) plotting tool, shared by all
  cases. Also on the CoMPhy Lab repo, merged into `main` (PR #6).
- **`figure-scripts/`** — scripts and small CSV inputs for manuscript figures
  that are not produced from the jet-base logs: the Fig. 1 schematic
  (`make_sketch_paper_worthington.py`) and the End Matter natural-numbered
  Fig. 4 diagnostic (`make_nu_vs_beta.py`, `beta_alpha_Oh.csv`).
- **`data-Oh-0.03/`** — the case we focus on (Bo=0, Oh=0.03). L13, L14, and two
  L15 runs (pre-inception focus 13 and 14). `bash reproduce-Oh0.03.sh` →
  `data-Oh-0.03/gridconv3_Oh0.03_L13-L14-L15.{png,pdf}`.
- **`facet-collapse-figure3/`** — Figure 3: the self-similar interface-collapse figure (pre/post inception, raw + rescaled by |t-t0|^alpha). `bash reproduce_facets.sh` -> `facet-collapse-figure3/fig3_Oh0.03_collapse.{png,pdf}`. See that folder's README.
- **`data-Oh-0.02297/`** — Bo=0, Oh=0.02297. L13 + L14. `bash reproduce-Oh0.02297.sh` →
  `data-Oh-0.02297/gridconv_Oh0.02297_L13-L14.{png,pdf}`.

Each `data-Oh-*/` folder holds `<case>_L<level>[_focus<f>]_log.txt` (+
`case.params`), one `facet_inception.txt` (the interface `x y` cloud at
inception, for the cone fit), the generated figure (the `reproduce-Oh*.sh` runners live at the folder root).

## Manuscript figure map

- **Fig. 1:** `sketch_paper_worthington.{pdf,png}`; regenerate with
  `python3 figure-scripts/make_sketch_paper_worthington.py`.
- **Fig. 2:** `data-Oh-0.03/gridconv3_Oh0.03_L13-L14-L15.{pdf,png}`;
  regenerate with `bash reproduce-Oh0.03.sh`.
- **Fig. 3:** `facet-collapse-figure3/fig3_Oh0.03_collapse.{pdf,png}`;
  regenerate with `bash reproduce_facets.sh`.
- **Fig. 4 / End Matter diagnostic:** `nu_vs_beta.pdf`; regenerate with
  `python3 figure-scripts/make_nu_vs_beta.py`.

### Grid / focus notation

`L<n>` = MAXlevel n (the finest refinement, reached ~when v_jet is highest).
`focus <m>` = the pre-inception cap: the max resolution allowed *just before*
the jet forms (the cavity-focus collapse is capped at level m, then released to
the full MAXlevel once the jet is issued). At Oh=0.03 the converged data is
insensitive to focus 13 vs 14 in the asymptotic window. The planned Level-15,
focus-15 run (5004) was cancelled once the small-`r_j` discrepancy was traced
to the historical plotting clip (`rmin = 0.04`), not to missing resolution.

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
- `q_jet`  (col 10) = **q_j = ∫ v_z r dr**.
- `q_l`    (col 11) = **q_ℓ = ∫ v_z dr**.

`Oh` and `Bond` are in each `*_case.params` (Bond only picks the initial shape;
the solver has no gravity term).

## Notes

- `--rmin 0.003` shows the full resolved range; the jet base resolves to
  r_j ≈ 0.005, so the historical 0.04 clip is cosmetic, not a resolution limit.
- Theory lines: cone (this work, α=1/(2−ν(β)), β fit from the inception facet),
  inertio-capillary (α=2/3), PRF 2023 (α=1/2). Fit prefactors over the r_j→0
  window `--fit-window 0.008 0.025`.
- If matplotlib complains about LaTeX, set `text.usetex = False` at the top of
  `plotJetMetricsTheory.py`.
