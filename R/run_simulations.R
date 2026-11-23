# Quadcopter (QC) block simulations for "When is a quadcopter the right tool for wildlife surveys?
# A simulation of survey performance across designs and population densities"
#
# Simulates distance-sampling surveys of 4, 9, 16 or 25 blocks (2 x 2 km) in a 10 x 10 km study area
# with three designs, at 10 densities, 999 replicates each (120 scenarios), and saves one dsims
# Simulation object per scenario in output/sims/. Run from the repository root; about 3 h on 28 cores.
#
# Set-up (details in R/setup.R):
#  * Nadir camera: certain detection within w = 34.7 m, analysed with the uniform key. Simulated with
#    0% overlap (69.4 m spacing) and 200 m spacing. 10% and 65% overlap record the same animals as
#    0% overlap, so they are not simulated; only their effort is computed.
#  * Gimbal camera: half-normal detection (sigma 160 m), w = 289.7 m, spacing 2w = 579.4 m, analysed
#    with the half-normal key.
#  * In every replicate a new population of N animals is drawn over the whole study area and the
#    block grid gets a random origin; only the animals inside the blocks are surveyed. The true
#    value is the study-area density N / 100 km^2.
#  * Fixed transects with every strip inside its block (transects/).
#  * Replicates without detections are estimates of 0 (dsims zero.n = "zero").
#  * Seeds: the population is seeded by (N, replicate), so all designs and block counts see the same
#    animals in replicate i; the detection process is seeded per scenario. dsims draws a seed for
#    every replicate, so the results do not depend on the number of cores.

source(here::here("R", "setup.R"))

reps <- 999                                  # replicates per scenario (e.g. 10 for a quick test)
cores <- 28                                  # parallel workers; does not change the results
block_counts <- c(4, 9, 16, 25)
abundances <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90)   # animals in the study area (x 10 = density per 1,000 km^2)
simulated <- c("0", "200", "gimbal")

dir.create(here("output", "sims"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("transects"), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), here("output", "sessionInfo.txt"))


# ---- transects and effort ------------------------------------------------------------------------
effort <- data.frame()
for (k in designs$design) {
  spacing <- designs$spacing[designs$design == k]; w <- designs$w[designs$design == k]
  xs <- line_positions(spacing, w)
  for (b in block_counts) {
    if (k %in% simulated) {
      f <- here("transects", sprintf("QC_%s_%d.shp", k, b))
      st_write(make_transects(b, spacing, w), f, delete_dsn = file.exists(f), quiet = TRUE)
    }
    effort <- rbind(effort, data.frame(
      design = k, blocks = b, spacing_m = spacing, truncation_m = w,
      lines_per_block = length(xs), n_transects = length(xs) * b,
      line_length_km = length(xs) * block_side * b / 1000,
      # flight path within a block: the lines plus the legs between neighbouring lines (turns ignored)
      within_block_trackline_km = (length(xs) * block_side + (length(xs) - 1) * spacing) * b / 1000,
      edge_margin_m = xs[1] - w))
  }
}
write.csv(effort, here("output", "effort.csv"), row.names = FALSE)


# ---- simulations ---------------------------------------------------------------------------------
for (b in block_counts) {
  blocks <- make_blocks(b)
  # make.simulation() needs a population description on the blocks; the populations themselves
  # come from the generator, so its N and density are not used.
  pop_desc_blocks <- make.population.description(region = blocks, density = make_density(blocks),
                                                 N = rep(1, b), fixed.N = TRUE)
  for (A in abundances) {
    generator <- make_generator(A, b)
    for (k in simulated) {
      id <- sprintf("QC_%s_design_%d-A%d", k, b, A)
      message(format(Sys.time(), "%H:%M "), id)
      w <- designs$w[designs$design == k]
      design <- make.design(region = blocks, transect.type = "line", design = "systematic",
                            samplers = numeric(0), line.length = numeric(0),
                            spacing = designs$spacing[designs$design == k],
                            design.angle = 0, edge.protocol = "minus", truncation = w)
      analysis <- make.ds.analysis(dfmodel = ~1, key = if (k == "gimbal") "hn" else "unif",
                                   criteria = "AIC", truncation = w, zero.n = "zero",
                                   group.strata = data.frame(design.id = as.character(seq_len(b)), analysis.id = "A"))
      sim <- make.simulation(reps = reps, design = design, population.description = pop_desc_blocks,
                             detectability = if (k == "gimbal") detect_gimbal else detect_nadir,
                             ds.analysis = analysis, population.generator = generator)

      seed <- 30000000 + match(k, designs$design) * 100000 + b * 1000 + A   # detection seed of the scenario
      set.seed(seed)
      t0 <- Sys.time()
      sim <- run.simulation(simulation = sim, run.parallel = cores > 1, max.cores = cores, counter = FALSE,
                            transect.path = here("transects", sprintf("QC_%s_%d.shp", k, b)))
      sim@results$filename <- sub(paste0("^", here(), "/"), "", sim@results$filename)  # repo-relative path

      QC_sim <- sim
      meta <- data.frame(id = id, design = k, blocks = b, abundance = A, seed_simulation = seed,
                         runtime_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
      save(QC_sim, meta, file = here("output", "sims", paste0(id, ".RData")), compress = "xz")
    }
  }
}
message("Done")
