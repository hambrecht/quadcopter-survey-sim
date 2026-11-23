# Quadcopter survey simulations

Code and outputs for the simulation study in "When is a quadcopter the right tool for wildlife
surveys? A simulation of survey performance across designs and population densities".

The study simulates distance-sampling surveys of wildlife in 2 × 2 km blocks with a quadcopter. It
compares three flight designs and 4–25 blocks, at densities of 50–900 individuals per 1,000 km², and
finds the least effort that gives a CV ≤ 0.15 with an absolute bias ≤ 10%.

All performance measures are computed in R, in the same environment as the simulations; Python
only draws the figures.

## How to run

Open the scripts in R (or RStudio) with the repository root as working directory and run them from
top to bottom, line by line if you like:

1. `renv::restore()` installs the pinned R packages (R 4.6.1), including the dsims fork.
2. `R/run_simulations.R` runs the 120 scenarios × 999 replicates and saves them in `output/sims/`
   (about 3 h on 28 cores). Set `reps` and `cores` at the top; the number of cores does not change
   the results. Skip this step to use the stored simulations.
3. `R/compute_metrics.R` computes all performance measures and writes `tables/` (about 30 s).
4. `python/make_figures.py` draws `figures/` from the tables (Python 3 with `requirements.txt`):
   `python3 python/make_figures.py`.

A re-run reproduces the simulated data (animals, detections, effort) and the nadir estimates exactly.
Gimbal estimates come from fitting a half-normal detection function and can differ from the stored
results at about the 6th significant digit (optimiser convergence); no reported number changes.

`R/validate_setup.R` (optional, about 2 min) checks the set-up; see [Checks](#checks).
`R/setup.R` holds the definitions shared by these scripts and is loaded by them.

## Repository layout

| Path | Content |
|---|---|
| `R/setup.R` | Definitions: camera geometry, designs, study area and blocks, fixed transects, density surface, detection functions, population generator |
| `R/run_simulations.R` | Runs the simulations → `output/`, `transects/` |
| `R/compute_metrics.R` | All performance measures, effort and decision tables → `tables/` |
| `R/validate_setup.R` | Checks of the set-up |
| `python/make_figures.py`, `python/style.py` | Figures from `tables/` → `figures/` (PNG, 300 dpi; for journal submission set `DPI = 600` and `TIFF = True` in `style.py`) |
| `python/palette_check.py` | Optional colour-vision check of the figure palette |
| `transects/QC_<design>_<blocks>.shp` | Fixed transects flown in every replicate |
| `output/sims/*.RData` | dsims Simulation objects (`QC_sim`) and scenario information (`meta`), one file per scenario, with every replicate result. The transect path stored in each object is relative to the repository root |
| `output/effort.csv` | Effort per design × number of blocks |
| `output/sessionInfo.txt` | R session of the run that produced `output/sims/` |
| `tables/`, `figures/` | Tables and figures of the manuscript |
| `renv.lock`, `requirements.txt` | R and Python package versions |

Scenario IDs are `QC_<design>_design_<blocks>-A<N>`, where N is the number of animals in the
100 km² study area (density = 10 × N per 1,000 km²).

## Simulation set-up

**Study area and blocks**

- The study area is 10 × 10 km. It has a smooth density surface: a ring of hotspots and one
  low-density area (`R/setup.R`).
- 4, 9, 16 or 25 blocks of 2 × 2 km are placed on a regular grid with spacing 10/√k km.
- The grid origin is drawn at random in every replicate, and the grid wraps at the study-area
  boundary. Every location therefore has the same probability of being surveyed.

**Population**

- In every replicate, N animals (fixed N) are drawn over the whole study area. Only the animals
  inside the blocks are surveyed.
- The true value is the density of the whole study area, so the CV includes the variation from
  surveying only part of it.
- Common random numbers: every design and block count sees the same animals and the same grid
  origin in replicate i.

**Designs**

| Design | Camera and detection | Truncation w | Spacing | Analysis |
|---|---|---|---|---|
| 0% overlap (nadir) | 640 × 512 thermal camera, HFOV 32.2°, 120 m AGL; certain detection within the strip | 34.7 m | 69.4 m | uniform key |
| 200 m spacing (nadir) | as above | 34.7 m | 200 m | uniform key |
| Gimbal | HFOV 25°, 55° lateral tilt; half-normal detection, σ = 160 m | 289.7 m | 579.4 m (2w) | half-normal key |

- The nadir designs with 10% overlap (62.5 m spacing) and 65% overlap (24.3 m) record the same
  animals as 0% overlap. They are not simulated separately; they use the 0% estimates with their
  own effort.
- Design, detection and analysis truncation are always identical. No adjustment terms are used.

**Transects and block edges**

- Fixed lines, no random start. In every block the same north–south lines are flown at the design
  spacing: as many as fit with every strip (2w wide) inside the block, centred so that the
  unsurveyed margin is split between the two edges.

  | Design | Spacing | Lines per block | Margin at each edge | Block area surveyed |
  |---|---|---|---|---|
  | 65% overlap | 24.3 m | 80 | 5 m | 99.5% |
  | 10% overlap | 62.5 m | 31 | 28 m | 97.2% |
  | 0% overlap | 69.4 m | 28 | 28 m | 97.2% |
  | 200 m spacing | 200 m | 10 | 65 m (to the strip) | 34.7% |
  | Gimbal | 579.4 m | 3 | 131 m | 86.9% |

- Nothing outside a block is viewed. No animal between blocks or in a neighbouring block can be
  counted, and none is counted twice. The surveyed area equals the effort used in the analysis
  (2wL), so there is no edge bias.
- The survey is still a random sample, because the block grid origin is random in every replicate.
- `transects/` holds the files for the simulated designs. Effort for the 10% and 65% designs is
  computed from the same rule (`line_positions()` in `R/setup.R`).

**Replicates without detections** give an estimate of 0 (dsims `zero.n = "zero"`). They are kept
in every estimate-based metric, not dropped as failures. There were no failed replicates.

## Software

- **R 4.6.1** with the package versions in `renv.lock`.
- **dsims** is the fork [hambrecht/dsims](https://github.com/hambrecht/dsims) at commit
  `281ed86` (version 1.2.0). It adds three things to dsims 1.1.0:
  - a `population.generator` hook (used for the whole-area population and the random grid
    origin);
  - `zero.n = "zero"` (estimate 0 for replicates without detections);
  - the uniform detection key in the analysis, with a fix for replicates with no animal near a
    transect.
- **dssd** is pinned to commit `e8100d1` (version 1.0.3).
- `renv.lock` pins both commits. `R/setup.R` reads the dsims commit from `renv.lock` and stops with an
  error if the installed dsims is a different commit (run `renv::restore()`).
- `output/sessionInfo.txt` records the session of the run. Its "installed from" field is empty,
  but the package came from the fork commit above.
- **Python 3** with numpy, pandas, matplotlib and pillow (`requirements.txt`).

## Outputs and definitions

Densities are in individuals per 1,000 km². The true density is N × 10: 50, 100, 200, …, 900.

**Metrics** (`R/compute_metrics.R`)

Where dsims provides a measure, it is taken from `summary(sim, description.summary = FALSE)` of the
saved Simulation object; everything else is computed in R over the same replicates dsims includes.
dsims densities (m⁻²) are multiplied by 10⁹. The dsims `Truth` is the per-replicate truth from the
population generator, i.e. the study-area density (the script checks this).

| Measure | Column | Source |
|---|---|---|
| True density | `true_density` | dsims `@individuals$D$Truth` |
| Mean estimate | `mean_estimate` | dsims `@individuals$D$mean.Estimate` |
| SD of the estimates | `sd_estimate` | dsims `@individuals$D$sd.of.means` |
| Bias (%) | `bias_pct` | dsims `@individuals$D$percent.bias` |
| RMSE | `rmse` | dsims `@individuals$D$RMSE` |
| Mean SE | `mean_se` | dsims `@individuals$D$mean.se` |
| 95% CI coverage | `ci_coverage` | dsims `@individuals$D$CI.coverage.prob` |
| Mean detections | `mean_detections` | dsims `@individuals$summary$mean.n` |
| Mean covered area | `mean_covered_area_km2` | dsims `@individuals$summary$mean.Cover.Area` |
| Replicates, failures, zero detections | `n_reps`, `n_failed`, `n_zero_detections` | dsims `@total.reps`, `@failures`, `@zero.detection.reps` |
| Empirical CV | `emp_cv` | R: `sd.of.means / mean.Estimate` |
| Relative RMSE | `rrmse` | R: `RMSE / Truth` |
| Median, 2.5% and 97.5% quantiles | `median_estimate`, `q025`, `q975` | R: `quantile(type = 7)` of the replicate estimates |
| Median bias (%) | `median_bias_pct` | R |
| Mean analytic CV | `mean_cv`, `n_cv` | R: mean of the replicate dsims CVs |
| Covered fraction, Poisson and binomial CV | `covered_fraction`, `cv_poisson`, `cv_binomial` | R: from `mean.Cover.Area` and `mean.n` |
| Flags, criteria and effort | | R |

- **Empirical CV** = SD / mean of the 999 replicate estimates.
- **Mean analytic CV** = mean of the dsims CVs. Replicates without detections have no CV.
- **95% CI coverage** = share of replicates with lcl < true < ucl, over replicates with a defined CI
  (dsims definition).
- **Bias** (`bias_pct`) = 100 × (mean − true) / true, signed. This is the bias used in the criteria.
- **Median bias** (`median_bias_pct`) = 100 × (median − true) / true, signed. Reported next to the mean
  bias because a few extreme gimbal estimates at low density inflate the mean.
- **Reliably estimable** (`reliably_estimable`, next to `mean_detections`) is FALSE for a **gimbal**
  scenario with fewer than 20 detections on average. A half-normal detection function cannot be fitted
  reliably from so few distances (Buckland et al. 2015 recommend 60–80). The nadir designs use a uniform
  key (count / covered area), fit no detection function and are never flagged. Flagged: gimbal with 4 and
  9 blocks (all densities), 16 blocks (50–500) and 25 blocks (50–300).
  - Figures 4, 5 and S2 draw flagged scenarios as open grey-edged diamonds (grey whiskers in Figure 4),
    not joined by the gimbal line.
  - Figures 5 and S2 cap the y-axis at 3 (above every non-flagged value); flagged values above the cap are
    drawn as a grey triangle at the cap with the value.
- **RMSE** = √mean((estimate − true)²). **Relative RMSE** = RMSE / true.

**Criteria** (Table 2 and Figure 6): empirical CV ≤ 0.15 and abs(mean bias) ≤ 10%, and the scenario
is reliably estimable: flagged gimbal scenarios never count as meeting the criteria.

- The Table 2 notes (`tables/Table_2_notes.txt`) give the exact CV of every scenario within 0.01
  of 0.15.
- Figure 6 plots the lowest simulated density from which the criteria hold at that density and at
  every higher simulated density. `tables/thresholds.csv` also gives the plain lowest density that
  meets them.
- Table 2 gives the smallest number of blocks that meets the criteria at both densities of a class.
- There is no interpolation, and nothing is extrapolated beyond 900.

**Effort** is reported as total transect length (km): lines per block × 2 km × number of blocks.
One block corresponds to one survey day. Covered area is the ground area inside the strips (2w wide;
overlapping strips count once); for the simulated designs it equals the dsims covered area. Flight time
and batteries are not outputs of the simulations and are not reported. `thresholds.csv` also keeps the
within-block flight path length (`flight_path_km`: lines plus the legs between them).

**Sampling-theory approximations of the CV** (nadir designs only; NA for the gimbal, whose estimate
also depends on the fitted detection function). With a fixed population size N and every animal inside
the surveyed strips with probability f, the nadir count is binomial:

- `covered_fraction` f = mean covered area / 100 km² (study area).
- `cv_binomial` = √((1 − f) / n), with n = mean number of detections.
- `cv_poisson` = 1 / √n, the limit of the binomial for small f.
- Figure 5 shows both for the 0% overlap design (grey: Poisson; black: binomial).
- `tables/approximation_check.csv`: over the 80 nadir scenarios the binomial approximation deviates
  from the empirical CV by at most 8.2% (median 1.6%); the Poisson approximation by a median of 18.6%
  and up to 515% (25 blocks, where f ≈ 0.97). The largest binomial deviations are at 4 blocks and high
  density, where all animals share one random block position over an uneven density surface.

**Low gimbal densities.** A few replicates with 1–3 detections give tiny P̂a values and very large
estimates. These drive the large positive mean bias, the CV and the 97.5% limits in Figures 4, 5 and
S2, while the median bias stays moderate. This is a property of fitting a detection function to very
few distances; these scenarios are flagged as not reliably estimable.

| Output | Content |
|---|---|
| `tables/scenario_metrics.csv` | All measures per design × blocks × density (input of the figures) |
| `tables/Table_S1_effort.csv` | Lines per block, transects, transect length, covered area, survey days |
| `tables/Table_2_decision.csv` | Smallest number of blocks meeting the criteria per density class |
| `tables/Table_2_notes.txt` | Notes to Table 2, with the exact CV of scenarios within 0.01 of 0.15 |
| `tables/Table_S2_scenario_metrics.csv` | Measures for the supplement |
| `tables/approximation_check.csv` | Deviation of the Poisson and binomial approximations from the empirical CV (nadir designs) |
| `tables/thresholds.csv` | Lowest density meeting the criteria per design × blocks, with transect length (Figure 6) |
| `figures/Figure_1` | Nadir camera geometry |
| `figures/Figure_2` | Detection functions (a) and gimbal geometry (b) |
| `figures/Figure_3` | Total transect length per design and number of blocks |
| `figures/Figure_4` | Estimated / true density: median and 2.5–97.5% range |
| `figures/Figure_5` | Empirical CV with the 0.15 line and the Poisson and binomial approximations (0% overlap) |
| `figures/Figure_6` | Lowest density meeting the criteria, against total transect length |
| `figures/Figure_S2` | Relative RMSE |

Figure numbers follow the manuscript. Table S1 replaces the former Supplementary Figure 2.

## Checks

`R/validate_setup.R` runs four checks (2,000 replicates at N = 90, no simulation output):

1. The mean density inside the blocks equals the study-area density for every number of blocks
   (equal inclusion probability).
2. Nadir and gimbal designs get identical animals in the same replicate (common random numbers).
3. The generated populations pass dsims validation and lie inside their block.
4. Every strip (2w) lies inside its own block, for every design and number of blocks.

## Licence

- **Code** (`R/`, `python/`): MIT, see `LICENSE`.
- **Outputs** (`output/`, `transects/`, `tables/`, `figures/`): CC BY 4.0, see `LICENSE-outputs`.
- **dsims and dssd** are GPL (≥ 2). They are installed as dependencies, not distributed here.

## Citation

See `CITATION.cff`.
