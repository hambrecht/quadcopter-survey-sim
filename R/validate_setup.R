# Checks of the simulation set-up in R/setup.R (optional; writes nothing). Run from the repository root.
#  1. Equal inclusion probability: the mean density inside the blocks equals the study-area density.
#  2. Common random numbers: nadir and gimbal designs get identical animals in the same replicate.
#  3. The generated populations pass the dsims checks.
#  4. Every strip (2w) lies inside its own block, for every design and block count.

source(here::here("R", "setup.R"))

A <- 90                    # animals in the study area
n_check <- 2000            # replicates for check 1

# a one-replicate simulation, only needed to call the generator
small_sim <- function(blocks, k) {
  w <- designs$w[designs$design == k]
  design <- make.design(region = blocks, transect.type = "line", design = "systematic", samplers = numeric(0),
                        line.length = numeric(0), spacing = designs$spacing[designs$design == k],
                        design.angle = 0, edge.protocol = "minus", truncation = w)
  make.simulation(reps = 1, design = design,
                  population.description = make.population.description(region = blocks, density = make_density(blocks),
                                                                        N = rep(1, length(blocks@strata.name)), fixed.N = TRUE),
                  detectability = if (k == "gimbal") detect_gimbal else detect_nadir,
                  ds.analysis = make.ds.analysis(dfmodel = ~1, key = if (k == "gimbal") "hn" else "unif", truncation = w))
}

results <- data.frame()
for (b in c(4, 9, 16, 25)) {
  blocks <- make_blocks(b)
  generator <- make_generator(A, b)
  sim_nadir <- small_sim(blocks, "0")
  sim_gimbal <- small_sim(blocks, "gimbal")

  # 1. animals inside the blocks over many replicates, against the expected number
  n_inside <- sapply(seq_len(n_check), function(r) nrow(generator(sim_nadir, r)$population@population))
  expected <- A * b * block_side^2 / study_side^2
  ratio <- mean(n_inside) / expected
  mc_se <- sd(n_inside) / sqrt(n_check) / expected

  # 2. and 3. the same animals for both designs, and valid populations
  same_animals <- TRUE
  for (r in 1:20) {
    p_nadir <- generator(sim_nadir, r)$population
    p_gimbal <- generator(sim_gimbal, r)$population
    dsims:::check.generated.population(p_nadir, sim_nadir, r)
    same_animals <- same_animals && identical(p_nadir@population[, c("x", "y", "Region.Label")],
                                              p_gimbal@population[, c("x", "y", "Region.Label")])
  }

  # 4. every strip inside its own block
  strips_inside <- TRUE
  for (k in designs$design) {
    w <- designs$w[designs$design == k]
    transects <- make_transects(b, designs$spacing[designs$design == k], w)
    strips <- st_buffer(transects, w * (1 - 1e-9), endCapStyle = "FLAT")
    for (i in seq_len(nrow(strips)))
      strips_inside <- strips_inside &&
        st_covered_by(strips[i, ], blocks@region[as.integer(transects$strata[i]), ], sparse = FALSE)[1, 1]
  }

  results <- rbind(results, data.frame(blocks = b, ratio_block_to_study_density = round(ratio, 4),
                                       mc_se = round(mc_se, 4), sd_animals_in_blocks = round(sd(n_inside), 2),
                                       same_animals = same_animals, strips_inside = strips_inside))
}
print(results)

stopifnot(abs(results$ratio_block_to_study_density - 1) < 4 * results$mc_se + 1e-9,
          results$same_animals, results$strips_inside,
          results$sd_animals_in_blocks[results$blocks == 25] == 0)   # 25 blocks cover the whole area
cat("All checks passed\n")
