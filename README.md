# Social Discounting Enables Fast and Reliable Collective Escape

MATLAB code for the paper. Reproduces every figure from scratch: the
drift-diffusion model, the moving-boundary first-passage solver, the Bayesian
fixed-point iteration, the Monte Carlo validation, and the empirical
identification from the sulphur molly field data.

Zachary P. Kilpatrick, Department of Applied Mathematics, University of
Colorado Boulder. `zpkilpat@colorado.edu`

## Requirements

MATLAB R2020b or later. Base MATLAB only for the empirical scripts, which
implement Wilson intervals, logistic regression, quantiles, and gamma
sampling locally rather than calling the Statistics Toolbox. `make_fig4.m`
calls `normcdf`, so that panel needs the toolbox unless you swap in the local
`ncdf` used elsewhere.

## Data

The field data are **not redistributed here**. Download Data S1 and Data S2
from the archive accompanying

> Pacher et al., *Better safe than sorry: the response of sulphur mollies to
> avian predators*, Science Advances (2025).

and place them in the working directory as

```
sciadv_adt8600_data_s1.csv     % per-event table, 258 scored events
sciadv_adt8600_data_s2.xlsx    % shoal structure, 76 measurement sessions
```

Data S1 is semicolon-delimited with a 19-line legend and the header on line
20. Two things about it are easy to get wrong and are handled in
`read_pacher.m`:

- Line 2 declares the sheet to be attacks only. It is not: the 81 flybys are
  in the same file with `risk = 'flyby'`. Below the real data sits a second
  stacked block (26 blanks, a repeated header, and 47 rows shifted one column
  left) which must be dropped. Filtering `risk` to `predator_attack` or
  `flyby` removes it and reproduces 177 attacks / 127 responses, 81 flybys /
  26 responses.
- `t_first` is in video frames at 25 fps, and 23 flybys carry a timing.
  Those are false-alarm response times, not detections. The pooling-count fit
  uses `D.tdet`, the 125 attack latencies only; using `D.tfirst` gives 148 and
  a different answer.

## Reproducing the figures

```matlab
mkdir output
run_all                      % everything, in dependency order
```

or individually:

```matlab
make_fig1                                     % cartoon panels B, C
make_fig2                                     % single-agent frontier
make_fig3                                     % dyad
bayes = run_bayes_nagent([2 5 10 20 40 100], 3.5678, 0.005);
out   = make_fig4(bayes, 'output');           % N-agent panels
out   = make_fig5('sciadv_adt8600_data_s1.csv');
make_figS1                                    % closed-form validation
out   = make_figS2;                           % MC validation, numbers only
out   = make_figS2(bayes, 'output');
out   = make_figS3('output');
out   = make_figS4('sciadv_adt8600_data_s1.csv', 'output');
out   = make_figS5('sciadv_adt8600_data_s1.csv', 'output');
out   = make_figS6('sciadv_adt8600_data_s1.csv', ...
                   'sciadv_adt8600_data_s2.xlsx', 'output');
```

Each script writes standalone panels (one PDF and one 400-dpi PNG per panel)
into `output/`; multi-panel figures are assembled externally. Panel letters
are deliberately not drawn in MATLAB.

### Order matters in two places

`run_bayes_nagent` must run **before** `make_fig4` and `make_figS2`, and at
the threshold `theta_1 = 3.5678`. Panel 4A reads the converged correction
from that struct; the other panels rebuild the threshold themselves.

`calib_pool` caches its calibration in a persistent variable. If you change
the response rates, call `calib_pool('reset')` or `clear functions` first, or
you will silently invert against the old curve.

## What each file does

**Solvers and shared machinery**

| file | role |
|---|---|
| `fpt_moving.m` | **the** first-passage primitive: Durbin--Buonocore Volterra sweep to a moving boundary |
| `survival_drift.m` | fixed threshold, time-varying drift; changes variables and calls `fpt_moving` |
| `solve_bayes.m` | damped fixed-point iteration for the survival correction |
| `solve_lambda_fixedpoint.m` | the dyad fixed point used by Fig. 3 |
| `run_bayes_nagent.m` | the production N-agent sweep |
| `calib_pool.m` | pooling-count calibration; the single place `M` is inverted |
| `read_pacher.m` | Data S1/S2 loader |
| `fit_pool_mle.m` | scale-free likelihood cross-check on `M` |
| `dyad_cascade_analytic.m` | closed-form naive-dyad cascade integrals (Eqs. S5--S9) |
| `dyad_frontier_heuristic.m` | dyad observables for the constant-alpha heuristic |
| `ddm_ig.m` | inverse-Gaussian first passage, single agent |
| `linf_ceiling.m` | closed-form saturation ceiling `L_inf(N)` |
| `alpha_hat.m` | Eq. (13), the rate inversion |
| `igrnd.m`, `wilson.m` | inverse-Gaussian draws; Wilson score interval |
| `paper_style.m`, `apply_axes.m`, `export_panel.m` | figure style for Figs. 1--3, S1 |

Everything that solves a first passage goes through `fpt_moving.m`. Three
separate implementations coexisted during development; see `MISSING.md`.

**Figure scripts**: `make_fig1`--`make_fig5`, `make_figS1`--`make_figS6`.

## Two conventions worth knowing

**The threshold is anchored at K = 1.** `theta_1 = -log(q_1)` is the
*solitary* threshold, so the solitary false alarm rate carries no `alpha`.
Holding the group rate there gives `theta_N = (theta_1 + log N)/(1+alpha)`.

**`q_1` is per unit, not the group rate.** The observed `q_N = 0.321` is
shoal-level, pooled over the `M` responders whose earliest crossing registers
as a group response. The solitary rate is `q_1 = 1-(1-q_N)^(1/M) = 0.0282` at
`Mhat = 13.52`, giving `theta_1 = 3.5678`. Earlier versions of `make_fig4`
used `q_N` directly and produced panels for a lone fish false-alarming a third
of the time.

## Headline numbers

Reproduced by `make_fig5` on Data S1:

```
TP = 0.7175   q = 0.3210   n_timed = 125   mean latency 4.918 s
lambda/Tbar = 3.928  ->  Mhat = 13.52
alpha-hat = 0.949    q_1 = 0.0282    theta_1 = 3.568
L_inf(Mhat) = 0.572  gap 0.377       K_max = 35.9
```

## Self-tests

Several scripts check themselves and will warn rather than fail silently.
`run_bayes_nagent` compares the solver against the exact inverse Gaussian at
zero coupling (max error ~8e-16). `make_figS2` checks the dyad cascade
probability against its analytic bounds `[q^2, 2q-q^2]` before reporting; an
earlier version of those kernels disagreed with Monte Carlo by 75% and passed
unnoticed because nothing checked it.

## License

MIT. See `LICENSE`.

## Citation

If you use this code, please cite the paper and this archive (Zenodo DOI in
the repository description).
