# Performance measures, survey effort and decision tables from the saved simulations.
# Run from the repository root, after R/run_simulations.R (about 30 s). Writes tables/:
#   scenario_metrics.csv          all measures per design x blocks x density (input of the figures)
#   Table_S1_effort.csv           effort: transects, transect length, covered area, survey days
#   Table_2_decision.csv          smallest number of blocks meeting the criteria per density class
#   Table_2_notes.txt             notes of Table 2, with the exact CV of scenarios near the 0.15 limit
#   Table_S2_scenario_metrics.csv measures for the supplement
#   thresholds.csv                lowest density meeting the criteria per design x blocks (Figure 6)
#   approximation_check.csv       deviation of the Poisson and binomial CV approximations (nadir designs)
#
# Where dsims provides a measure it is taken from summary(sim): Truth (the study-area density from
# the population generator), mean.Estimate, sd.of.means, percent.bias, RMSE, mean.se,
# CI.coverage.prob (replicates with a defined CI) and mean.n. The rest is computed here over the
# same replicates dsims includes. Densities in individuals per 1,000 km^2 (dsims m^-2 x 1e9).
# Effort is reported as transect length (km); one block corresponds to one survey day.

library(here)
library(dsims)

cv_max <- 0.15              # criteria: empirical CV <= 0.15 and abs(mean bias) <= 10%
bias_max <- 10
min_detections <- 20        # gimbal scenarios below this mean are not reliably estimable
design_name <- c("65" = "65% overlap", "10" = "10% overlap", "0" = "0% overlap", "200" = "200 m", "gimbal" = "Gimbal")
dir.create(here("tables"), showWarnings = FALSE)


# ---- measures per scenario -----------------------------------------------------------------------
rows <- list()
for (f in list.files(here("output", "sims"), pattern = "\\.RData$", full.names = TRUE)) {
  load(f)                                                    # QC_sim, meta
  s <- summary(QC_sim, description.summary = FALSE)
  D <- s@individuals$D
  truth <- D$Truth * 1e9
  stopifnot(grepl("population.generator", s@truth.source[["D"]]), all.equal(truth, meta$abundance * 10))
  # replicates included by dsims: successful fits and zero-detection replicates
  used <- dsims:::get.rep.index(QC_sim@results$Detection, QC_sim@reps, 1, FALSE)
  est <- QC_sim@results$individuals$D[1, "Estimate", used] * 1e9
  cv <- QC_sim@results$individuals$D[1, "cv", used]
  all_reps <- seq_len(QC_sim@reps)
  n_in_blocks <- QC_sim@results$individuals$true.N[1, "N", all_reps]          # animals inside the blocks
  q <- quantile(est, c(0.5, 0.025, 0.975), names = FALSE)
  rows[[f]] <- data.frame(
    design = meta$design, blocks = meta$blocks, true_density = truth,
    n_reps = s@total.reps, n_success = s@total.reps - s@failures, n_failed = s@failures,
    n_zero_detections = s@zero.detection.reps,
    mean_estimate = D$mean.Estimate * 1e9, median_estimate = q[1], q025 = q[2], q975 = q[3],
    sd_estimate = D$sd.of.means * 1e9,
    emp_cv = D$sd.of.means / D$mean.Estimate,               # empirical CV
    mean_cv = mean(cv, na.rm = TRUE), n_cv = sum(!is.na(cv)),
    mean_se = D$mean.se * 1e9, ci_coverage = D$CI.coverage.prob,
    n_ci = sum(!is.na(QC_sim@results$individuals$D[1, "lcl", used]) & !is.na(QC_sim@results$individuals$D[1, "ucl", used])),
    bias_pct = D$percent.bias, median_bias_pct = 100 * (q[1] - truth) / truth,
    rmse = D$RMSE * 1e9, rrmse = D$RMSE * 1e9 / truth,      # relative RMSE
    mean_detections = s@individuals$summary$mean.n,
    mean_N_in_blocks = mean(n_in_blocks), sd_N_in_blocks = sd(n_in_blocks),
    mean_true_D_blocks = mean(n_in_blocks / QC_sim@results$individuals$summary[1, "Area", all_reps]) * 1e9,
    mean_covered_area_km2 = s@individuals$summary$mean.Cover.Area / 1e6)
}
m <- do.call(rbind, rows)
# Gimbal: a half-normal detection function cannot be fitted reliably from fewer than 20 detections
# (Buckland et al. 2015 recommend 60-80). Nadir designs use the uniform key and are never flagged.
m$reliably_estimable <- !(m$design == "gimbal" & m$mean_detections < min_detections)
m$rel_median <- m$median_estimate / m$true_density
m$rel_q025 <- m$q025 / m$true_density
m$rel_q975 <- m$q975 / m$true_density
m$failed_gt5pct <- m$n_failed / m$n_reps > 0.05
m$meets_criteria <- m$emp_cv <= cv_max & abs(m$bias_pct) <= bias_max & m$reliably_estimable
# Sampling-theory approximations of the CV for the nadir designs (count / covered area): with fixed N
# and every animal inside the strips with probability f, the count is binomial, CV = sqrt((1 - f) / n);
# for small f this is the Poisson 1 / sqrt(n). NA for the gimbal (fitted detection function).
m$covered_fraction <- m$mean_covered_area_km2 / 100                # study area 100 km^2
nadir <- m$design != "gimbal"
m$cv_poisson <- ifelse(nadir, 1 / sqrt(m$mean_detections), NA)
m$cv_binomial <- ifelse(nadir, sqrt((1 - m$covered_fraction) / m$mean_detections), NA)
m <- m[order(match(m$design, names(design_name)), m$blocks, m$true_density),
       c("design", "blocks", "true_density", "n_reps", "n_success", "n_failed", "n_zero_detections",
         "mean_estimate", "median_estimate", "q025", "q975", "sd_estimate", "emp_cv",
         "covered_fraction", "cv_poisson", "cv_binomial", "mean_cv", "n_cv",
         "mean_se", "ci_coverage", "n_ci", "bias_pct", "median_bias_pct", "rmse", "rrmse",
         "mean_detections", "reliably_estimable", "mean_N_in_blocks", "sd_N_in_blocks",
         "mean_true_D_blocks", "mean_covered_area_km2", "rel_median", "rel_q025", "rel_q975",
         "failed_gt5pct", "meets_criteria")]
rownames(m) <- NULL
write.csv(m, here("tables", "scenario_metrics.csv"), row.names = FALSE, na = "")


# ---- effort (Table S1) ---------------------------------------------------------------------------
e <- read.csv(here("output", "effort.csv"), colClasses = c(design = "character"))
eff <- data.frame(design = e$design, blocks = e$blocks, spacing_m = e$spacing_m,
                  lines_per_block = e$lines_per_block, n_transects = e$n_transects,
                  transect_length_km = e$line_length_km,
                  flight_path_km = e$within_block_trackline_km)   # lines plus the legs between them
# covered area: union of the strips (2w wide) in each 2 km block; overlapping strips count once
eff$covered_area_km2 <- (2 * e$truncation_m + (e$lines_per_block - 1) * pmin(e$spacing_m, 2 * e$truncation_m)) *
  2000 * e$blocks / 1e6
eff$survey_days <- eff$blocks                                    # one block per day
eff <- eff[order(match(eff$design, names(design_name)), eff$blocks), ]
# for the simulated designs this equals the dsims covered area
for (i in which(eff$design %in% m$design)) {
  sim_area <- m$mean_covered_area_km2[m$design == eff$design[i] & m$blocks == eff$blocks[i]]
  stopifnot(all.equal(sim_area, rep(eff$covered_area_km2[i], length(sim_area))))
}
s1 <- eff[, c("design", "blocks", "spacing_m", "lines_per_block", "n_transects", "transect_length_km",
              "covered_area_km2", "survey_days")]
s1$design <- design_name[s1$design]
write.csv(s1, here("tables", "Table_S1_effort.csv"), row.names = FALSE)


# ---- thresholds (Figure 6) -----------------------------------------------------------------------
# Lowest simulated density from which the criteria hold at that and every higher density.
# 65% and 10% overlap record the same animals as 0% overlap: they use the 0% measures with their own effort.
thr <- data.frame()
for (k in names(design_name)) for (b in c(4, 9, 16, 25)) {
  src <- if (k %in% c("65", "10")) "0" else k
  s <- m[m$design == src & m$blocks == b, ]
  met <- s$meets_criteria; dens <- s$true_density
  from_here <- which(sapply(seq_along(met), function(i) all(met[i:length(met)])))
  ef <- eff[eff$design == k & eff$blocks == b, ]
  thr <- rbind(thr, data.frame(design = k, blocks = b, metrics_from = src,
                               threshold_density = dens[from_here[1]], lowest_density_met = dens[met][1],
                               transect_length_km = ef$transect_length_km, flight_path_km = ef$flight_path_km))
}
write.csv(thr, here("tables", "thresholds.csv"), row.names = FALSE, na = "")


# ---- Table 2: smallest number of blocks meeting the criteria at both densities of a class ------------
classes <- list("50-100" = c(50, 100), "200-300" = c(200, 300), "400-500" = c(400, 500),
                "600-700" = c(600, 700), "800-900" = c(800, 900))
t2 <- data.frame(density_class = names(classes), check.names = FALSE)
for (k in names(design_name)) {
  src <- if (k %in% c("65", "10")) "0" else k
  t2[[design_name[[k]]]] <- sapply(classes, function(dens) {
    for (b in c(4, 9, 16, 25)) {
      s <- m[m$design == src & m$blocks == b & m$true_density %in% dens, ]
      km <- eff$transect_length_km[eff$design == k & eff$blocks == b]
      if (all(s$meets_criteria)) return(sprintf("%d blocks (%s km)", b, formatC(km, format = "d", big.mark = ",")))
    }
    "not met"
  })
}
write.csv(t2, here("tables", "Table_2_decision.csv"), row.names = FALSE)

near <- m[abs(m$emp_cv - cv_max) <= 0.01 + 1e-12, ]
writeLines(c(
  "Notes to Table 2:",
  "- Each cell gives the smallest number of 2 x 2 km blocks that meets empirical CV <= 0.15 and abs(mean bias) <= 10% at both densities of the class (individuals per 1,000 km^2), with the total transect length in brackets (one block = one survey day). 65% and 10% overlap use the 0% overlap estimates with their own effort.",
  "- Gimbal scenarios with fewer than 20 detections on average are not reliably estimable (a half-normal detection function cannot be fitted reliably from so few distances; Buckland et al. 2015 recommend 60-80) and never count as meeting the criteria. The nadir designs use a uniform key (count / covered area) and are not flagged.",
  "- Scenarios with an empirical CV within 0.01 of 0.15 (exact CV):",
  sprintf("  - %s, %d blocks, %.0f: CV = %.3f, mean bias = %+.1f%% (%s)", design_name[near$design], near$blocks,
          near$true_density, near$emp_cv, near$bias_pct,
          ifelse(near$meets_criteria, "meets the criteria", "does not meet the criteria"))),
  here("tables", "Table_2_notes.txt"))


# ---- Table S2 ------------------------------------------------------------------------------------
s2 <- m[, c("design", "blocks", "true_density", "n_success", "n_failed", "n_zero_detections",
            "mean_estimate", "median_estimate", "q025", "q975", "sd_estimate", "emp_cv",
            "covered_fraction", "cv_poisson", "cv_binomial", "mean_cv",
            "mean_se", "ci_coverage", "bias_pct", "median_bias_pct", "rmse", "rrmse", "mean_detections",
            "reliably_estimable")]
s2$design <- design_name[s2$design]
write.csv(s2, here("tables", "Table_S2_scenario_metrics.csv"), row.names = FALSE)


# ---- how well the approximations reproduce the empirical CV (nadir designs) -------------------------
nad <- m[m$design != "gimbal", ]
dev_binomial <- abs(nad$cv_binomial - nad$emp_cv) / nad$emp_cv
dev_poisson <- abs(nad$cv_poisson - nad$emp_cv) / nad$emp_cv
approx <- data.frame(approximation = c("cv_binomial", "cv_poisson"), n_scenarios = nrow(nad),
                     max_rel_deviation = c(max(dev_binomial), max(dev_poisson)),
                     median_rel_deviation = c(median(dev_binomial), median(dev_poisson)))
print(approx)
write.csv(approx, here("tables", "approximation_check.csv"), row.names = FALSE)
