# Definitions shared by R/run_simulations.R and R/validate_setup.R: camera geometry, designs, study
# area, blocks, transects, density surface, detection functions and the population generator.
# Run from the repository root. This script only defines objects; it writes nothing.

library(here)
library(dsims)
library(sf)

# dsims is the hambrecht/dsims fork pinned in renv.lock (install with renv::restore()). The fork adds
# the population.generator hook, zero.n = "zero" and the uniform key.
dsims_commit <- renv::lockfile_read(here("renv.lock"))$Packages$dsims$RemoteSha
if (!identical(packageDescription("dsims")$RemoteSha, dsims_commit))
  stop("Installed dsims is not the commit in renv.lock; run renv::restore()")


# ---- camera geometry and designs ------------------------------------------------------------
# Nadir camera: 640 x 512 thermal, HFOV 32.2 deg at 120 m: w = 120 * tan(16.1 deg) = 34.64 m,
# rounded to 34.7 m as in the manuscript.
w_nadir <- 34.7
# Gimbal camera: HFOV 25 deg at 55 deg lateral tilt: w = 120 * tan(55 + 12.5 deg) = 289.7 m.
w_gimbal <- round(120 * tan((55 + 25 / 2) * pi / 180), 1)
sigma_gimbal <- 160

# Line spacing: nadir designs with 65%, 10% and 0% overlap of the 69.4 m swath, 200 m spacing,
# and the gimbal with spacing 2w (no overlap). w is the truncation distance.
designs <- data.frame(
  design  = c("65", "10", "0", "200", "gimbal"),
  spacing = c(round(2 * w_nadir * (1 - 0.65), 1),   # 24.3 m
              round(2 * w_nadir * (1 - 0.10), 1),   # 62.5 m
              2 * w_nadir,                          # 69.4 m
              200,
              2 * w_gimbal),                        # 579.4 m
  w       = c(w_nadir, w_nadir, w_nadir, w_nadir, w_gimbal))


# ---- study area and blocks ----------------------------------------------------------------------
# 10 x 10 km study area. Blocks of 2 x 2 km at [j * pitch, j * pitch + 2000], pitch = 10000 / sqrt(n),
# so the gap between blocks equals the gap across the boundary. The random grid origin of each
# replicate is applied to the animals (population generator below), so the blocks never move.
# Blocks are numbered with x outer, y inner: block = ix * sqrt(n) + iy + 1.
study_side <- 10000
block_side <- 2000
study_square <- matrix(c(0, 0, study_side, 0, study_side, study_side, 0, study_side, 0, 0), ncol = 2, byrow = TRUE)
region <- make.region(region.name = "study area", shape = st_polygon(list(study_square)), units = "m")
st_crs(region@region) <- 32633

make_blocks <- function(n_blocks) {
  n_side <- sqrt(n_blocks); pitch <- study_side / n_side
  squares <- list()
  for (i in 0:(n_side - 1)) for (j in 0:(n_side - 1)) {
    x0 <- i * pitch; y0 <- j * pitch
    squares[[length(squares) + 1]] <- st_polygon(list(rbind(
      c(x0, y0), c(x0, y0 + block_side), c(x0 + block_side, y0 + block_side), c(x0 + block_side, y0), c(x0, y0))))
  }
  blocks <- st_sf(geometry = st_sfc(squares))
  st_crs(blocks) <- st_crs(region@region)
  make.region(region.name = "study area", strata.name = as.character(seq_len(n_blocks)), shape = blocks, units = "m")
}


# ---- fixed transects -----------------------------------------------------------------------------
# No random start: in every block the same north-south lines are flown at the design spacing, as
# many as fit with every strip (2w) inside the block, centred so the unsurveyed margin is split
# between the two edges. Nothing outside a block is viewed. The survey is still a random sample of
# the study area because the block grid origin is random in every replicate.
line_positions <- function(spacing, w) {          # x positions of the lines from the block's west edge
  room <- block_side - 2 * w
  k <- floor(room / spacing + 1e-9) + 1
  w + (room - (k - 1) * spacing) / 2 + (seq_len(k) - 1) * spacing
}

make_transects <- function(n_blocks, spacing, w) {
  n_side <- round(sqrt(n_blocks)); pitch <- study_side / n_side; xs <- line_positions(spacing, w)
  lines <- list(); block <- character(0)
  for (ix in 0:(n_side - 1)) for (iy in 0:(n_side - 1)) {
    x0 <- ix * pitch; y0 <- iy * pitch
    for (x in xs) lines[[length(lines) + 1]] <- st_linestring(rbind(c(x0 + x, y0), c(x0 + x, y0 + block_side)))
    block <- c(block, rep(as.character(ix * n_side + iy + 1), length(xs)))
  }
  transects <- st_sf(transect = seq_along(lines), strata = block, geometry = st_sfc(lines))
  st_crs(transects) <- st_crs(region@region)
  transects
}


# ---- density surface: ring of hotspots and one low-density area ---------------------------------
make_density <- function(reg) {
  centres <- list(
    c(5000, 3000), c(5789, 3108), c(6551, 3421), c(7236, 3906), c(7765, 4509),
    c(8039, 5165), c(7968, 5803), c(7506, 6354), c(6678, 6758), c(5591, 6972),
    c(4408, 6972), c(3321, 6758), c(2493, 6354), c(2031, 5803), c(1960, 5165),
    c(2234, 4509), c(2763, 3906), c(3448, 3421), c(4210, 3108), c(4999, 3000))
  d <- make.density(region = reg, x.space = 500, constant = 1)
  for (centre in centres) d <- add.hotspot(object = d, centre = centre, sigma = 1500, amplitude = 2)
  add.hotspot(object = d, centre = c(5000, 8000), sigma = 1500, amplitude = -4)
}
density_study <- make_density(region)


# ---- detection -------------------------------------------------------------------------------------
# Nadir: certain detection within the strip. Gimbal: half-normal.
detect_nadir  <- make.detectability(key.function = "uf", scale.param = 1, truncation = w_nadir)
detect_gimbal <- make.detectability(key.function = "hn", scale.param = sigma_gimbal, truncation = w_gimbal)


# ---- population generator ------------------------------------------------------------------------
# dsims calls generator(simulation, rep) in every replicate. It draws N = A animals over the whole
# study area (seeded by A and rep only, so every design and block count sees the same animals in
# replicate rep), shifts them by a random grid origin (wrapping at the boundary) and returns the
# animals inside the blocks. The truth is the study-area density A / 100 km^2.
# Everything the generator uses is copied into its own environment because dsims runs it in
# separate worker sessions.
make_generator <- function(A, n_blocks) {
  pop_desc_study <- make.population.description(region = region, density = density_study, N = A, fixed.N = TRUE)
  study_region <- region; side <- study_side; bside <- block_side
  study_area <- sum(region@area)
  n_side <- round(sqrt(n_blocks)); pitch <- side / n_side
  function(simulation, rep) {
    set.seed(50000000 + A * 10000 + rep)          # population seed; dsims restores its RNG state afterwards
    offset <- stats::runif(2, 0, side)            # random grid origin
    pop <- dsims::generate.population(object = pop_desc_study, detectability = simulation@detectability,
                                      region = study_region)
    df <- pop@population
    x <- (df$x - offset[1]) %% side; y <- (df$y - offset[2]) %% side
    ix <- floor(x / pitch); iy <- floor(y / pitch)
    inside <- (x - ix * pitch) < bside & (y - iy * pitch) < bside
    df <- df[inside, , drop = FALSE]
    df$x <- x[inside]; df$y <- y[inside]
    df$Region.Label <- as.character(ix[inside] * n_side + iy[inside] + 1)
    df$individual <- seq_len(nrow(df))
    rownames(df) <- NULL
    pop@population <- df
    pop@N <- as.numeric(table(factor(df$Region.Label, levels = simulation@design@region@strata.name)))
    pop@strata.names <- simulation@design@region@strata.name
    list(population = pop, truth = list(N = A, D = A / study_area))
  }
}
