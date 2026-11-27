# 02_estimates.R
# Coverage estimates for every sample: Penta1, Penta3, OPV1, OPV3, MCV1, MCV2
# and zero-dose (no Penta1), from the children reported directly by respondents
# Reproduces: no report item; writes data/clean/estimates.rds, read by 04-11

# Setup ----

pacman::p_load(survey, RDS, tidyverse)

source("code/functions.R")

# Data ----

frames <- read_rds("data/clean/frames.rds")
tracking <- read_rds("data/clean/rds_tracking.rds")

# Analysis ----

## Analytic cohorts ----
# Enrolled respondents' children (every CBE record is a child), after the
# analytic exclusions, with zero-dose as no Penta1.
# keep each sample's directly reported children
cohorts <- imap(frames[study_samples], direct_cohort)

## Household survey and CBE ----
# Two-stage cluster design (cluster, then household or record) with equal
# weights: clusters were drawn with probability proportional to population and
# fully enumerated.
# estimate each antigen in each HH and CBE sample
cluster_estimates <- map(str_subset(names(cohorts), "_(hh|cbe)$"), \(sample) {
  cohort <- cohorts[[sample]]
  if (str_detect(sample, "cbe")) cohort <- mutate(cohort, cluster_id = cluster, record_ix = record_id)
  design <- svydesign(id = ~ cluster_id + record_ix, weights = ~weight, nest = TRUE, data = mutate(cohort, weight = 1))
  map(antigen_outcomes, \(field) {
    survey_proportion(design, field) |>
      mutate(n_children = nrow(cohort), n_observed = sum(!is.na(cohort[[field]])))
  }) |>
    list_rbind(names_to = "antigen") |>
    mutate(sample = sample, n_respondents = if (str_detect(.env$sample, "hh")) n_distinct(cohort$record_ix) else NA_integer_)
}) |>
  list_rbind()

## VBS ----
# Respondents weighted by the inverse of how often they travel to the venue
# town where they were met, each respondent a cluster of their children.
# weight respondents and estimate each antigen in each VBS sample
vbs_estimates <- map(str_subset(names(cohorts), "vbs"), \(sample) {
  cohort <- cohorts[[sample]]
  cohort$weight <- vbs_weight(cohort, sample)
  design <- svydesign(id = ~record_ix, weights = ~weight, data = cohort)
  map(antigen_outcomes, \(field) {
    survey_proportion(design, field) |>
      mutate(n_children = nrow(cohort), n_observed = sum(!is.na(cohort[[field]])))
  }) |>
    list_rbind(names_to = "antigen") |>
    mutate(sample = sample, n_respondents = n_distinct(cohort$record_ix))
}) |>
  list_rbind()

## RDS ----
# Successive-sampling estimator on the recruitment network, with the protocol
# population size and a seed for each sample and outcome. Somalia Phase 2
# completes the network from the supervisors' tracking sheet. Zero-dose shares
# the Penta1 seed, so the two estimates sum to 1.
# build each recruitment network and estimate each antigen
rds_estimates <- map(str_subset(names(cohorts), "rds"), \(sample) {
  cohort <- cohorts[[sample]]
  nodes <- rds_network(cohort, antigen_outcomes, tracking = if (sample == "som_phase2_rds") tracking)
  design <- as.rds.data.frame(
    nodes,
    id = "id_participant_calc", recruiter.id = "id_referring",
    network.size = "degree", population.size = population_size[[sample]]
  )
  map(antigen_outcomes, \(field) {
    seed_field <- if (field == "zero_dose") "penta_1p_clean" else field
    seed <- derived_seed(sample_seed(sample), paste0("rds_outcome/", seed_field))
    rds_proportion(design, cohort, field, population_size[[sample]], seed) |>
      mutate(n_children = nrow(cohort))
  }) |>
    list_rbind(names_to = "antigen") |>
    mutate(sample = sample, n_respondents = n_distinct(cohort$id_participant_calc), n_network_nodes = nrow(nodes))
}) |>
  list_rbind()

## Combine ----
# One row per sample and antigen, methods and antigens in report order.
# stack the three designs
estimates <-
  bind_rows(cluster_estimates, vbs_estimates, rds_estimates) |>
  separate_wider_regex(sample, c(country = "[a-z]+", "_phase", phase = "\\d", "_", method = "[a-z]+"), cols_remove = FALSE) |>
  mutate(
    phase = as.integer(phase),
    method = factor(method, levels = c("hh", "cbe", "rds", "vbs"), labels = c("HH", "CBE", "RDS", "VBS")),
    antigen = factor(antigen, levels = names(antigen_outcomes))
  ) |>
  arrange(phase, country, method, antigen) |>
  select(sample, country, phase, method, antigen, estimate, standard_error, ci_lower, ci_upper, n_children, n_observed, n_respondents, n_network_nodes)

# Save ----

write_rds(estimates, "data/clean/estimates.rds")
