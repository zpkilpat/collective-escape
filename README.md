# Social Discounting Enables Fast and Reliable Collective Escape

Code to reproduce every figure panel in the paper and its supplementary
materials.

Kilpatrick, Z. P. (2026). *Social Discounting Enables Fast and Reliable
Collective Escape.* Archived at Zenodo, DOI `PLACEHOLDER`.

---

## Quick start

```matlab
run_all
```

Panels are written to `output/` as PDF and PNG. Multi-panel figures are
assembled externally, so panel letters are not drawn here.

Runtime is dominated by one step, the `N = 100` fixed point, which takes
twenty-five damped iterations of two `O(nt^2)` Volterra sweeps each. Expect
tens of minutes for a full run. Setting `SKIP_SOLVER = true` in `run_all.m`
loads the solver struct that the first complete run saves to
`bayes_theta3568.mat`.

MATLAB only, no toolboxes:

```matlab
[~, p] = matlab.codetools.requiredFilesAndProducts('make_fig5.m');  p.Name
```

returns `MATLAB` alone. Local functions inside scripts require R2016b or
later. `caxis` is used rather than `clim` so the code runs before R2022a.

---

## Data

The field data are not redistributed. `run_all` skips the empirical figures
with a warning if they are absent.

| File | Contents |
|---|---|
| `sciadv_adt8600_data_s1.csv` | Per-event table, 177 kiskadee attacks and 81 flybys |
| `sciadv_adt8600_data_s2.xlsx` | Shoal structure, 76 sessions across three sites |

Both are Data S1 and S2 of Pacher et al., *Science Advances* **11**, eadt8600
(2025). Download them from that paper's archive into the working directory.

Two features of Data S1 that the loaders handle. The header sits on line 20
beneath nineteen lines of legend, the file is semicolon-delimited, and a
second stacked block below the real data must be dropped. And `t_first` is
recorded in video frames at 25 fps, not seconds.

---

## Files

### Entry point

`run_all.m` regenerates everything in dependency order.

### Figure scripts

One per figure, each writing its own panels. `make_fig1`, `make_fig2`,
`make_fig3` and `make_figS1` are scripts; the rest are functions returning a
struct of what they computed.

| Script | Figure | Needs |
|---|---|---|
| `make_fig1` | Fig. 1B, C | nothing |
| `make_fig2` | Fig. 2 | nothing |
| `make_fig3` | Fig. 3 | nothing |
| `make_fig4` | Fig. 4 | `bayes` from `solve_bayes` |
| `make_fig5` | Fig. 5 | Data S1 |
| `make_figS1` | Fig. S1 | nothing |
| `make_figS2` | Fig. S2 | `bayes` from `solve_bayes` |
| `make_figS3` | Fig. S3 | nothing |
| `make_figS4` | Fig. S4 | Data S1 |
| `make_figS5` | Fig. S5 | Data S1 |
| `make_figS6` | Fig. S6 | Data S1 and S2 |

### Shared primitives

Called by the figure scripts, not run directly.

| File | Does |
|---|---|
| `solve_bayes` | Bayesian survival-correction fixed point for an N-agent group |
| `solve_lambda_fixedpoint` | The same closure, as used by Fig. 3 |
| `fpt_moving` | Durbin-Buonocore moving-boundary first-passage solver |
| `survival_drift` | Survival correction on a given drift, over `fpt_moving` |
| `calib_pool` | Monte Carlo calibration of the pooling count, cached in a persistent |
| `read_pacher` | Data S1 and S2 loader, builds the recording-by-bout cluster key |
| `alpha_hat` | Eq. (13), discounting rate from two response rates at a pooling count |
| `linf_ceiling` | Closed-form saturation ceiling `L_inf(N)` |
| `fit_pool_mle` | Scale-free likelihood cross-check on the pooling count |
| `ddm_ig` | Inverse-Gaussian first-passage density, survival and hazard |
| `igrnd` | Inverse-Gaussian variates, Michael-Schucany-Haas |
| `wilson` | Wilson score interval |
| `paper_style`, `apply_axes`, `export_panel` | Plotting helpers for the script-style figures |

### Validation, run by hand

`validate_dyad_mc.m` checks the dyad observables against an independent
Euler-Maruyama simulation sharing no code with the Volterra quadrature.
Nothing calls it and it draws nothing by default, since the supplementary
figure it once produced was cut. It is retained because the data-availability
statement names the Monte Carlo validation and because the closed-form
agreement it reproduces is quoted in the supplement.

---

## Reproducing the numbers

`run_all` prints checks as it goes. The ones worth watching:

```
cluster key: 73 clusters, 47 attack, 26 flyby
TP = 0.7175, q = 0.3210, stat = 3.928 -> Mhat = 13.5, alpha = 0.9489
bootstrap: M in [8.82, 29.40], alpha in [0.906, 0.982]
           4000/4000 valid, P(alpha > L_inf) = 1.0000, min excess = 0.2043
```

The cluster count is a guard, not a diagnostic. Both loaders assert it,
because an earlier column-index error keyed clusters on the recording
*location* rather than the recording *file*, giving 31 clusters with six
carrying flybys and silently invalidating every clustered interval in the
paper. If it prints anything but 73/47/26, stop.

Three quantities in the paper are all called "the gap" and are not
interchangeable:

| Value | Meaning |
|---|---|
| 0.376 | `alpha_hat - L_inf(M)` at the point estimate |
| 0.260 | minimum over the admissible wedge `K <= M`, attained on `K = M` |
| 0.204 | minimum over the 4000 bootstrap replicates |

Bootstrap intervals are reproducible. All three resamples share
`seed = 20260811` and 4000 replicates, and the seed is applied immediately
before the resampling loop rather than at the top of the file, since
`calib_pool` advances the random stream enormously on a cold cache.

---

## Design notes

**One first-passage primitive.** Three implementations of the same Volterra
sweep coexisted during development, differing only in where the minus sign
sat in the kernel. Since that sign convention was once wrong in the
supplement, shipping three copies was the wrong call. `fpt_moving.m` accepts
handles or vectors and is the only one; `survival_drift` and both fixed-point
solvers route through it, so their self-tests exercise the same code every
figure uses.

**Shared primitives extracted.** `linf_ceiling`, `igrnd`, `wilson` and
`alpha_hat` were defined locally in up to six files each. One file apiece
now, called everywhere.

**Illustrative thresholds in the expository figures.** `make_fig1` and
`make_fig3` set their own thresholds rather than the fitted `theta_1`, since
both exist to show mechanism at a scale where individual crossings are
visible. The captions state the values used, and no number in the paper is
read off those panels. Only `make_figS1` needed the fitted anchor, because
its hump scaling runs through `theta_N`.

**Small helpers left duplicated on purpose.** `igd`/`igpdf`, `ncdf`,
`quantile_local`, `irls_logit`, the shape statistic, `assert_key`, and the
per-script export boilerplate. MATLAB's one-public-function-per-file rule
makes factoring these out a net loss in navigability. The two `assert_key`
copies differ only in naming the loader each one guards.

---

## Known inconsistencies

Documented rather than fixed, since the figures are verified as they stand
and refactoring before submission would risk more than it buys.

**Two style systems.** `make_fig1`, `make_fig2`, `make_fig3` and `make_figS1`
use `paper_style` with `apply_axes` and `export_panel`; the rest carry an
inline `P` struct with a local `export`. Both produce the same house style,
22 by 19 cm with LaTeX interpreters throughout.

**Two Data S1 loaders.** `make_fig5` parses the CSV itself; everything else
calls `read_pacher`. They disagreed during development, 148 latencies against
125, because one had not yet excluded the twenty-three answered flybys. Both
are correct now and both assert the cluster key independently.

**Two fixed-point iterations.** `solve_lambda_fixedpoint` serves Fig. 3;
`solve_bayes` serves Fig. 4 and Fig. S2. They differ in damping and in their
late-time survival floor, and agree on the saturation tail to about 0.3
percent.

**Discretization.** The time step is `Tmax/nt` with `nt = 4000` fixed while
`Tmax` grows with `N`, so the dyad has the coarsest grid. Doubling `nt` moves
the `N = 2` tail from 0.1636 to 0.1638 and leaves `N = 100` unchanged at
0.8077, which is why the paper quotes those recoveries to one significant
figure.

---

## License

See `LICENSE`. The Pacher et al. data are not covered by it and are not
redistributed here.
