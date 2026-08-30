# Social discounting enables fast and reliable collective escape

Code to reproduce every figure panel in the paper and its SI Appendix.

Kilpatrick, Z. P. (2026). *Social discounting enables fast and reliable
collective escape.*

---

## Quick start

```matlab
run_all
```

Panels are written to `output/` as PDF and PNG. Multi-panel figures are
assembled externally, so panel letters are not drawn here. `PANELMAP.md` says
which script writes which panel — worth reading first, because several scripts
do not emit the letter their name suggests.

Runtime is dominated by one step, the `N = 100` fixed point, which takes
twenty-five damped iterations of two `O(nt^2)` Volterra sweeps each. Expect
tens of minutes for a full run. Setting `SKIP_SOLVER = true` in `run_all.m`
loads the solver struct that the first complete run saves to
`bayes_theta3568.mat`.

MATLAB only, no toolboxes:

```matlab
[~, p] = matlab.codetools.requiredFilesAndProducts('make_fig4.m');  p.Name
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

Three features of Data S1 that the loaders handle. The header sits on line 20
beneath nineteen lines of legend and the file is semicolon-delimited. A second
stacked block below the real data must be dropped. And `t_first` is recorded
in video frames at 25 fps, not seconds.

---

## Files

### Entry point

`run_all.m` regenerates everything in dependency order.

### Figure scripts

One per figure, each writing its own panels. `make_fig1`, `make_fig2`,
`make_figS1` and `make_figS3` are scripts; the rest are functions returning a
struct of what they computed.

| Script | Writes | Needs |
|---|---|---|
| `make_fig1` | Fig. 1B, C, D | nothing |
| `make_fig2` | Fig. 2A–D **and Fig. S2C** | nothing |
| `make_fig3` | Fig. 3A–D | `bayes` from `solve_bayes` |
| `make_fig4` | Fig. 4A–D | Data S1 |
| `make_figS1` | Fig. S1A, B | nothing |
| `make_figS2` | Fig. S2A, B | `bayes` from `solve_bayes` |
| `make_figS3` | Fig. S3A, B | nothing |
| `make_figS4` | Fig. S4A–C | nothing (simulation study) |
| `make_figS5` | Fig. S5A–C | Data S1 |
| `make_figS6` | Fig. S6A–C | Data S1 |
| `make_figS7` | Fig. S7A–C | Data S1 and S2 |

### Solvers and model code

| File | Role |
|---|---|
| `fpt_moving.m` | Durbin–Buonocore Volterra sweep, moving boundary |
| `solve_bayes.m` | damped fixed point for `lambda^(N)`, N-agent |
| `solve_lambda_fixedpoint.m` | dyad fixed point used by `make_fig2` |
| `survival_drift.m` | survival under a time-varying drift |
| `ddm_ig.m` | inverse-Gaussian first-passage density, survival, hazard |
| `dyad_cascade_analytic.m` | closed-form cascade integrals G and C |
| `dyad_frontier_heuristic.m` | heuristic (α, θ) family for the frontier |
| `validate_dyad_mc.m` | Monte Carlo cross-check of the dyad observables |

### Empirical primitives

| File | Role |
|---|---|
| `read_pacher.m` | shared Data S1 / S2 loader, `file|bout` cluster key |
| `calib_pool.m` | self-consistent M calibration at `theta_1(M)` |
| `alpha_hat.m` | Eq. (3), α̂ from two pooled rates |
| `fit_pool_mle.m` | scale-free likelihood cross-check on M |
| `linf_ceiling.m` | closed-form `L_inf(N)` |
| `igrnd.m` | inverse-Gaussian sampler (Michael–Schucany–Haas) |
| `wilson.m` | Wilson score interval |

### Style

`paper_style.m`, `apply_axes.m`, `export_panel.m`.

---

## Expected console output

A clean run reproduces these. Any drift means something upstream moved.

```
solve_bayes tails      0.1636 0.3661 0.5005 0.6151 0.7098 0.8077
  against L_inf(N)     0.1716 0.3820 0.5195 0.6345 0.7269 0.8182
  iterations           12 10 10 13 17 25

Data S1 cluster key    73 clusters, 47 attack, 26 flyby
rates                  TP = 0.7175   q = 0.3210   n_timed = 125
latency                mean 4.9184 s   lambda_IG 19.3194   ratio 3.9280
identification         Mhat 13.5241   alpha-hat 0.9489
                       q_1 0.0282     theta_1 3.5678     K_max 35.93
bootstrap (4000)       Mci [8.825, 29.403]   alphaci [0.906, 0.982]
                       pOver 1.0000
```

**Three distinct gaps, all correct, not interchangeable.** `0.376` is
α̂ − L∞(M̂) at the point estimate; `0.260` is the minimum over the admissible
(M, K) wedge, attained on K = M; `0.204` is the minimum over the 4000
bootstrap replicates.

**Seed convention.** `rng(seed)` is applied immediately before the bootstrap
call, never at the top of a script. `calib_pool` draws 2e6 samples per grid
point on a cold cache and advances the stream enormously, so a seed set at the
top does not align across files.

---

## Known inconsistencies

Left in deliberately, documented rather than fixed, because changing them
risks moving panels that currently reproduce.

- **Two style systems.** `paper_style`/`apply_axes`/`export_panel` for the
  script-style figures; an inline `P` struct for the function-style ones.
- **Two fixed-point iterations.** `solve_lambda_fixedpoint` (dyad, used by
  `make_fig2`) and `solve_bayes` (N-agent). They agree to about 0.3%; the
  paper's numbers come from `solve_bayes`.
- **Illustrative thresholds.** `make_fig1` runs panels B and C at
  θ = 3.5678 and D at θ = 1, 2, 3; `make_fig2` panels A and B at θ = 1, 2, 3.
  These are stated in the captions and no reported number is read off them.
- **Discretization.** Panels are generated at `nt = 4000`. Halving the step
  moves the N = 2 tail from 0.1636 to 0.1638 and leaves N = 100 unchanged, so
  the shortfall against `L_inf` is quoted to one significant figure.
- `make_fig4` takes `(csv, opts)` where `opts` may be an options struct or an
  output-directory string; `make_figS4`–`S7` take `(csv, outdir)`.

---

## License

MIT. See `LICENSE`.
