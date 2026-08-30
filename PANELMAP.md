# Panel map

Which script writes which panel. Several scripts do not emit the letter their
own name suggests, because the paper went from five main figures and six
supplementary to four and seven without renaming panels inside the code.

## Main text

| Paper | File | Written by | Content |
|---|---|---|---|
| Fig. 1A | — | Illustrator | attack / flyby cartoon |
| Fig. 1B | `fig1B` | `make_fig1` | DDM sample paths |
| Fig. 1C | `fig1C` | `make_fig1` | naive vs discounting on shared noise |
| Fig. 1D | `fig1D` | `make_fig1` | survival LLR, θ = 1,2,3 |
| Fig. 2A | `fig2A` | `make_fig2` (internal B) | cascade probability ratio vs kick |
| Fig. 2B | `fig2B` | `make_fig2` (internal C) | cascade delay vs kick |
| Fig. 2C | `fig2C` | `make_fig2` (internal E) | survival, naive / heuristic / Bayes |
| Fig. 2D | `fig2D` | `make_fig2` (internal F) | speed–accuracy frontier |
| Fig. 3A–D | `fig3A`–`fig3D` | `make_fig3` | N-agent theory |
| Fig. 4A | `fig4A` | `make_fig4` | response rates vs shoal area |
| Fig. 4B | `fig4B` | `make_fig4` | latency distribution + IG fit |
| Fig. 4C | `fig4C` | `make_fig4` | calibration of λ_IG/T̄ against M |
| Fig. 4D | `fig4D` | `make_fig4` | α̂(M) against the admissible optimum |

## Supplement

| Paper | File | Written by | Content |
|---|---|---|---|
| Fig. S1A, S1B | `figS1A/B` | `make_figS1` | single-agent density, frontier |
| Fig. S2A, S2B | `figS2A/B` | `make_figS2` | heuristic vs Bayesian fixed point |
| Fig. S2C | `figS2C` | **`make_fig2`** (internal D) | converged dyad correction λ(t) |
| Fig. S3A, S3B | `figS3A/B` | `make_figS3` | saturation ceiling, hump scaling |
| Fig. S4A–C | `figS4A–C` | `make_figS4` | pooling-count robustness |
| Fig. S5A–C | `figS5A–C` | `make_figS5` | over-discounting wedge |
| Fig. S6A–C | `figS6A–C` | `make_figS6` | rate estimation and intervals |
| Fig. S7A–C | `figS7A–C` | `make_figS7` | area, headcount, density |

Fig. S2C is emitted by `make_fig2` because it shares that script's dyad
fixed-point solve; putting it in `make_figS2` would cost a second Volterra
sweep for one curve.

## Where each came from

Against the Science Advances version: old Fig. 2B → Fig. 1D; old Fig. 2A, 2C →
Fig. S1; old Fig. 3B, 3C, 3E, 3F → Fig. 2A–D; old Fig. 3D → Fig. S2C; old
Fig. 3A dropped; old Fig. 4 → Fig. 3; old Fig. 5 → Fig. 4 with the calibration
inset promoted to a standalone panel C. Supplement: old S1 → S3, old S3 → S4,
old S4 → S6, old S5 → S5, old S6 → S7, S2 unchanged. (SI order is set by first
mention in the main text; S4 and S6 were swapped again after that audit.)
