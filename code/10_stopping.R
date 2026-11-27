# 10_stopping.R
# Stopping point: how many respondents each Phase 1 sample needed before its
# Penta3 estimate (and its OPV3 and MCV1 estimates) was both stable and precise
# Reproduces: no report item; writes data/clean/stopping.rds, the stopping
# points that 11_cost.R reports in Table 14, Figures 18 and 21 and Annex Table
# L2, and the RDS paths that 07_rds_vbs.R draws in Annex Figure E2
# Runtime: about 30 minutes on 15 cores (several hours on one)

# Setup ----

pacman::p_load(survey, RDS, parallel, tidyverse)

source("code/functions.R")

# Data ----

frames <- read_rds("data/clean/frames.rds")
estimates <- read_rds("data/clean/estimates.rds")

# the Phase 1 samples, the tracer antigen and the two other primary antigens
# followed along the same orders, and the outcomes traced along the realised
# RDS orders
samples <- c(
  "eth_phase1_hh", "eth_phase1_cbe", "eth_phase1_rds", "eth_phase1_vbs",
  "som_phase1_hh", "som_phase1_cbe", "som_phase1_rds", "som_phase1_vbs"
)
rds_samples <- c("eth_phase1_rds", "som_phase1_rds")
tracer <- "penta_3p_clean"
order_outcomes <- c(tracer, "opv_3p_clean", "measles_1p_clean")
realised_outcomes <- c(tracer, "opv_3p_clean", "measles_1p_clean", "zero_dose")

# stopping criterion: the estimate stays within 5 percentage points of its
# final value and the 95% CI is at most 10 percentage points wide, counted
# from 10% of the sample (at least 50 respondents)
max_drift <- 0.05
max_ci_width <- 0.10
floor_share <- 0.10
floor_min <- 50

# 200 random enrolment orders per sample; RDS, whose estimator is slow, uses
# 100 orders and refits every 10th respondent
n_orders <- c(hh = 200, cbe = 200, rds = 100, vbs = 200)
n_step <- c(hh = 1, cbe = 1, rds = 10, vbs = 1)

# a VBS interview site used on at most a quarter of the field days is an
# excursion, and its consecutive days move together
excursion_share <- 0.25

# a modelled stopping point is not estimated beyond 1.5 times the sample
extrapolation_cap <- 1.5

# fits run in parallel, one R session per core but one
cores <- max(1L, detectCores() - 1L)

# Analysis ----

## Enrolment order ----
# Directly reporting respondents in the order they were enrolled: by field
# day, then by the time consent was completed, then by record. CBE forms were
# filled on paper and keyed in later, so within a day their order is by record;
# CBE records with no recoverable field day come last.
# order the respondents of one sample
enrolment_order <- function(sample) {
  respondents <- analytic_respondents(frames[[sample]], sample)
  if (str_detect(sample, "cbe")) respondents$t_end_consent <- as.POSIXct(NA)
  arrange(respondents, is.na(date), date, is.na(t_end_consent), t_end_consent, key)
}

# order every sample
orders <- map(set_names(samples), enrolment_order)

# each sample's directly reported children, keyed by respondent
cohorts <- map(set_names(samples), \(sample) {
  cohort <- direct_cohort(frames[[sample]], sample)
  cohort$key <- if (str_detect(sample, "cbe")) cohort$record_id else cohort$record_ix
  cohort
})

## Prefix estimates ----
# One outcome among the children of a set of respondents, from each method's
# own estimator, returned as the estimate, its 95% CI and the CI width in
# percentage points. In RDS, a respondent whose recruiter is not in the set is
# treated as a seed, and zero-dose shares the Penta1 seed.
# estimate an outcome among the children of the respondents in `keys`
fit_prefix <- function(sample, keys, outcome = tracer) {
  cohort <- cohorts[[sample]][cohorts[[sample]]$key %in% keys, ]
  method <- str_extract(sample, "[a-z]+$")
  if (method %in% c("hh", "cbe")) {
    if (method == "cbe") cohort <- mutate(cohort, cluster_id = cluster, record_ix = record_id)
    design <- svydesign(id = ~ cluster_id + record_ix, weights = ~weight, nest = TRUE, data = mutate(cohort, weight = 1))
    result <- survey_proportion(design, outcome)
  } else if (method == "vbs") {
    cohort$weight <- vbs_weight(cohort, sample)
    design <- svydesign(id = ~record_ix, weights = ~weight, data = cohort)
    result <- survey_proportion(design, outcome)
  } else {
    orphan <- cohort$seed_participant == "No" & !cohort$id_referring %in% cohort$id_participant_calc
    cohort$seed_participant[orphan] <- "Yes"
    cohort$id_referring[orphan] <- NA
    nodes <- rds_network(cohort, outcome)
    design <- as.rds.data.frame(
      nodes,
      id = "id_participant_calc", recruiter.id = "id_referring",
      network.size = "degree", population.size = population_size[[sample]]
    )
    seed_field <- if (outcome == "zero_dose") "penta_1p_clean" else outcome
    seed <- derived_seed(sample_seed(sample), paste0("rds_outcome/", seed_field))
    result <- suppressWarnings(rds_proportion(design, cohort, outcome, population_size[[sample]], seed))
  }
  c(estimate = result$estimate, ci_lower = result$ci_lower, ci_upper = result$ci_upper, ci_width_pp = 100 * (result$ci_upper - result$ci_lower))
}

# sample sizes at which the estimate is refitted, from the floor to the full sample
prefix_grid <- function(n_final, step = 1) {
  start <- ceiling(max(floor_min, floor_share * n_final))
  unique(c(seq(start, n_final, by = step), n_final))
}

# refit along one enrolment order (a vector of respondent keys)
fit_path <- function(sample, keys, grid, outcome) {
  map(grid, \(n) fit_prefix(sample, keys[seq_len(n)], outcome)) |>
    reduce(rbind) |>
    as_tibble() |>
    mutate(n = grid, .before = 1)
}

# refit the first n respondents of the realised order, or along a random order
fit_realised <- function(n, sample, outcome) fit_prefix(sample, orders[[sample]]$key[seq_len(n)], outcome)
fit_order <- function(perm, sample, grid, outcome) fit_path(sample, orders[[sample]]$key[perm], grid, outcome)

## Stopping criteria ----
# The stability point is the first sample size after the last time the
# estimate was more than 5 percentage points from its final value; the
# precision point is the first after the last time the 95% CI was wider than
# 10 percentage points. The stopping point is the later of the two, and is
# missing when either criterion is still failing at the full sample.
first_stable <- function(n, estimate) {
  outside <- which(abs(estimate - estimate[length(estimate)]) > max_drift)
  if (length(outside) == 0) n[1] else if (max(outside) >= length(n)) NA_real_ else n[max(outside) + 1]
}

first_precise <- function(n, ci_width_pp) {
  above <- which(ci_width_pp / 100 > max_ci_width)
  if (length(above) == 0) n[1] else if (max(above) >= length(n)) NA_real_ else n[max(above) + 1]
}

stopping_points <- function(path) {
  stability_n <- first_stable(path$n, path$estimate)
  precision_n <- first_precise(path$n, path$ci_width_pp)
  tibble(stability_n, precision_n, stopping_n = if_else(is.na(stability_n) | is.na(precision_n), NA_real_, pmax(stability_n, precision_n)))
}

## Parallel sessions ----
# Each session loads the estimators and the cohorts once.
# start the sessions and give them what the fits need
sessions <- makeCluster(cores)
invisible(clusterEvalQ(sessions, {
  pacman::p_load(survey, RDS, tidyverse)
  source("code/functions.R")
}))
clusterExport(sessions, c("cohorts", "orders", "tracer", "fit_prefix", "fit_path", "fit_realised", "fit_order"))

## Realised enrolment order ----
# The estimate and its 95% CI refitted at every sample size along the order in
# which respondents were actually enrolled: Penta3 in every sample, and in RDS
# also OPV3, MCV1 and zero-dose.
# refit one sample and outcome along its realised order
realised_path <- function(sample, outcome) {
  grid <- prefix_grid(nrow(orders[[sample]]))
  fits <- parLapplyLB(sessions, grid, fit_realised, sample = sample, outcome = outcome)
  tibble(
    sample = sample, outcome = outcome, n = grid, estimate = map_dbl(fits, "estimate"),
    ci_lower = map_dbl(fits, "ci_lower"), ci_upper = map_dbl(fits, "ci_upper"), ci_width_pp = map_dbl(fits, "ci_width_pp")
  )
}

# refit Penta3 in every sample and the other outcomes in RDS
realised <-
  c(
    map(samples, \(sample) realised_path(sample, tracer)),
    map(rds_samples, \(sample) map(setdiff(realised_outcomes, tracer), \(outcome) realised_path(sample, outcome))) |> list_flatten()
  ) |>
  list_rbind()

## Random enrolment orders ----
# Which respondents were enrolled is held fixed and only their order is
# redrawn, in a way that keeps each method's field logistics:
# - HH: each cluster-day keeps its respondents, cluster-days are reassigned to
#   the realised daily slots (no cluster twice on one day), and respondents
#   are shuffled within each day;
# - CBE: respondents are shuffled within each field day;
# - RDS: seeds keep their place, and each later respondent is drawn at random
#   from those whose recruiter has already been enrolled;
# - VBS: field days are shuffled as whole blocks (the consecutive days of an
#   excursion site stay together, in order), respondents within each day.
# Each order has its own seed derived from the study seed.
# shuffle a vector of row positions
shuffle <- function(rows) rows[sample.int(length(rows))]

# field day of each respondent in realised order (undated CBE records last)
field_day <- function(order) {
  day <- as.integer(factor(order$date))
  day[is.na(day)] <- max(c(0L, day), na.rm = TRUE) + 1L
  day
}

# HH: cluster-days to daily slots, then respondents within days
permute_hh <- function(order) {
  day <- field_day(order)
  unit_key <- paste(day, order$cluster_id)
  unit <- match(unit_key, unique(unit_key))
  n_units <- max(unit)
  first_row <- match(seq_len(n_units), unit)
  unit_day <- day[first_row]
  unit_cluster <- order$cluster_id[first_row]
  unit_rows <- split(seq_along(unit), factor(unit, levels = seq_len(n_units)))
  slots <- sort(unit_day)
  for (attempt in 1:1000) {
    assigned <- integer(n_units)
    assigned[sample.int(n_units)] <- slots
    if (!anyDuplicated(paste(assigned, unit_cluster))) break
  }
  map(seq_len(max(day)), \(d) shuffle(unlist(unit_rows[which(assigned == d)], use.names = FALSE))) |> unlist()
}

# CBE: respondents within each field day
permute_cbe <- function(order) {
  day <- field_day(order)
  split(seq_along(day), factor(day, levels = seq_len(max(day)))) |>
    map(shuffle) |>
    unlist(use.names = FALSE)
}

# VBS: blocks of field days (a run of consecutive days at an excursion site
# is one block), then respondents within each day
venue_day_blocks <- function(day, site) {
  days <- sort(unique(day))
  block <- set_names(seq_along(days), days)
  for (s in unique(na.omit(site))) {
    site_days <- sort(unique(day[!is.na(site) & site == s]))
    if (length(site_days) / length(days) > excursion_share) next
    runs <- split(site_days, cumsum(c(1, diff(site_days) != 1)))
    for (run in runs[lengths(runs) >= 2]) {
      merged <- unique(block[as.character(run)])
      block[block %in% merged] <- min(merged)
    }
  }
  split(days, match(block, sort(unique(block))))
}

permute_vbs <- function(order) {
  day <- field_day(order)
  day_rows <- split(seq_along(day), factor(day, levels = seq_len(max(day))))
  blocks <- venue_day_blocks(day, order$int_location)
  perm <- integer(0)
  for (block in blocks[sample.int(length(blocks))]) {
    for (d in sort(block)) perm <- c(perm, shuffle(day_rows[[d]]))
  }
  perm
}

# RDS: seeds in place, recruits drawn from those whose recruiter is enrolled
permute_rds <- function(order) {
  n <- nrow(order)
  is_seed <- is.na(order$id_referring)
  recruiter <- match(order$id_referring, order$id_participant_calc)
  recruits <- split(which(!is_seed), factor(recruiter[!is_seed], levels = seq_len(n)))
  perm <- integer(n)
  available <- integer(0)
  for (slot in seq_len(n)) {
    if (is_seed[slot]) {
      pick <- slot
    } else {
      j <- sample.int(length(available), 1)
      pick <- available[j]
      available[j] <- available[length(available)]
      available <- available[-length(available)]
    }
    perm[slot] <- pick
    available <- c(available, recruits[[pick]])
  }
  perm
}

# draw the orders of one sample, each under its own seed
draw_orders <- function(sample) {
  method <- str_extract(sample, "[a-z]+$")
  permute <- switch(method, hh = permute_hh, cbe = permute_cbe, rds = permute_rds, vbs = permute_vbs)
  map(seq_len(n_orders[[method]]), \(b) {
    seed <- derived_seed(study_seed, sprintf("stopping_rule_resampling/%s/%04d", str_replace(sample, "phase", "phase_"), b))
    withr::with_seed(seed, permute(orders[[sample]]), .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  })
}

# draw every sample's orders
permutations <- map(set_names(samples), draw_orders)

# refit each sample along each of its random orders, for Penta3, OPV3 and MCV1
resampled <- imap(permutations, \(perms, sample) {
  method <- str_extract(sample, "[a-z]+$")
  grid <- prefix_grid(nrow(orders[[sample]]), n_step[[method]])
  map(order_outcomes, \(outcome) {
    paths <- parLapplyLB(sessions, perms, fit_order, sample = sample, grid = grid, outcome = outcome)
    list_rbind(set_names(paths, seq_along(paths)), names_to = "order") |> mutate(outcome = outcome, order = as.integer(order), .before = 1)
  }) |>
    list_rbind()
})

stopCluster(sessions)

## Stopping point over orders ----
# The stopping point of each random order, summarised by its median and its
# 10th and 90th percentiles. An order that never meets the criterion counts
# at the full sample, and a percentile that falls among those orders is
# missing (so the median is missing unless at least half of the orders meet
# the criterion).
# stopping point of every order and outcome
order_stopping <-
  resampled |>
  imap(\(paths, sample) {
    paths |>
      group_by(outcome, order) |>
      group_modify(\(path, key) stopping_points(path)) |>
      ungroup() |>
      mutate(sample = sample, .before = 1)
  }) |>
  list_rbind()

# median and range of the Penta3 stopping point over orders
observed <-
  order_stopping |>
  filter(outcome == tracer) |>
  left_join(map_dfr(orders, \(o) tibble(n_final = nrow(o)), .id = "sample"), by = "sample") |>
  summarise(
    stopping_n_q10 = order_percentile(stopping_n, 0.1, n_final[1]),
    stopping_n_q90 = order_percentile(stopping_n, 0.9, n_final[1]),
    stability_n_median = order_percentile(stability_n, 0.5, n_final[1]),
    stopping_n = order_percentile(stopping_n, 0.5, n_final[1]),
    .by = sample
  )

## Modelled stopping point ----
# Where the Penta3 interval was still wider than 10 percentage points at the
# end of fieldwork, the stopping point is projected: holding the estimate and
# the design effect fixed, the standard error falls with the square root of
# the number of sampling units (clusters for HH and CBE, venue-days for VBS,
# respondents for RDS), and the Wilson interval uses the design's own t
# critical value. The modelled stopping point is the smallest whole number of
# units at which the interval is at most 10 points wide (and at least the
# stability point), times the realised respondents per unit; beyond 1.5 times
# the sample it is not estimable.
# Wilson interval width at an effective sample size
wilson_width <- function(p, se, crit) {
  n_eff <- p * (1 - p) / se^2
  (crit / n_eff) * sqrt(4 * n_eff * p * (1 - p) + crit^2) / (1 + crit^2 / n_eff)
}

# sampling units: clusters, venue-days (any interview at a venue on a day,
# whether or not it found an eligible caregiver) or respondents
sampling_units <- function(sample) {
  method <- str_extract(sample, "[a-z]+$")
  switch(
    method,
    hh = n_distinct(orders[[sample]]$cluster_id),
    cbe = n_distinct(orders[[sample]]$cluster),
    rds = nrow(orders[[sample]]),
    vbs = frames[[sample]] |> distinct(record_ix, .keep_all = TRUE) |> filter(!is.na(date), !is.na(int_location)) |> distinct(int_location, date) |> nrow()
  )
}

# project the stopping point of one sample
model_stopping <- function(sample, stability_n) {
  method <- str_extract(sample, "[a-z]+$")
  primary <- filter(estimates, sample == .env$sample, antigen == "Penta3")
  units <- sampling_units(sample)
  per_unit <- nrow(orders[[sample]]) / units
  df <- function(u) switch(method, hh = , cbe = u - 1, vbs = u * per_unit - 1, rds = Inf)
  crit <- function(u) if (method == "rds") 1.96 else qt(0.975, pmax(df(u), 1))
  candidates <- (units + 1):floor(extrapolation_cap * units + 1e-9)
  width <- map_dbl(candidates, \(u) wilson_width(primary$estimate, primary$standard_error * sqrt(units / u), crit(u)))
  precise_units <- candidates[which(width <= max_ci_width + 1e-9)[1]]
  stable_units <- if (is.na(stability_n)) 0 else ceiling(stability_n / per_unit - 1e-9)
  model_units <- max(precise_units, stable_units)
  if (is.na(model_units) || model_units > extrapolation_cap * units + 1e-9) model_units <- NA_real_
  tibble(n_units_final = units, respondents_per_unit = per_unit, n_units_at_stopping = model_units, modelled_n = model_units * per_unit)
}

## Stopping point ----
# Observed where the median order met the criterion within the sample,
# modelled otherwise.
# combine the observed and modelled stopping points
stopping <-
  observed |>
  mutate(
    n_final = map_int(sample, \(s) nrow(orders[[s]])),
    map2_dfr(sample, stability_n_median, model_stopping),
    stopping_basis = case_when(!is.na(stopping_n) ~ "observed", !is.na(modelled_n) ~ "modelled", .default = "not_estimable"),
    n_units_at_stopping = if_else(stopping_basis == "observed", stopping_n / respondents_per_unit, n_units_at_stopping),
    stopping_n = if_else(stopping_basis == "modelled", modelled_n, stopping_n)
  ) |>
  left_join(
    estimates |>
      filter(antigen == "Penta3") |>
      transmute(sample, country, method = str_to_lower(method), final_estimate = estimate, final_ci_width_pp = 100 * (ci_upper - ci_lower)),
    by = "sample"
  ) |>
  select(
    sample, country, method, n_final, final_estimate, final_ci_width_pp, stopping_basis, stopping_n, stopping_n_q10,
    stopping_n_q90, stability_n_median, n_units_final, respondents_per_unit, n_units_at_stopping
  )

# Save ----

write_rds(
  list(stopping = stopping, realised = realised, order_stopping = order_stopping),
  "data/clean/stopping.rds"
)
