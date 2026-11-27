# 08_sensitivity.R
# Sensitivity analyses: each analytic choice replaced by an alternative, and the
# Phase 1 and Phase 2 verdicts that change
# Reproduces: Tables 11 and 13, Annex Tables H1, H3-H6 and I1-I3
# Runtime: about 12 minutes on 15 cores, about an hour on one (the raking bootstraps)

# Setup ----

pacman::p_load(survey, RDS, mice, parallel, tidyverse)

source("code/functions.R")

# Data ----

frames <- read_rds("data/clean/frames.rds")
estimates <- read_rds("data/clean/estimates.rds")
tracking <- read_rds("data/clean/rds_tracking.rds")
card_photo_status <- read_rds("data/clean/card_photo_status.rds")

# the fourteen samples and the four outcomes compared
samples <- setdiff(names(frames), "som_phase2_rds_all_enumerators")
outcomes <- antigen_outcomes[primary_outcomes]

# Analysis ----

## Seeds ----
# Every sample draws from its own seed, derived from the study seed; each
# analysis derives its seeds from that, keyed on the sample name.
# name of a sample in the seed keys
sample_key <- function(sample) str_replace(sample, "_phase(\\d)_", "_phase_\\1_")

## Cohorts ----
# Each sample's analytic cohort, with the columns its design needs: the cluster
# and record of the HH and CBE designs (and the HH cluster label used as an
# imputation predictor), and the VBS mobility weights.
# weight VBS respondents by the inverse of their visits to the venue town;
# respondents who made no visit are scored as 1-4 visits (lowest_positive),
# scored 1 (floor) or left out (exclude); unweighted gives every respondent 1
weight_vbs <- function(cohort, sample, zero_trip = "lowest_positive", weighted = TRUE) {
  travel <- vbs_travel(cohort, sample)
  if (zero_trip == "exclude") {
    cohort <- cohort[!travel %in% "0", ]
    travel <- travel[!travel %in% "0"]
  }
  multiplier <- travel_multiplier(cohort, sample)
  zero_score <- if (zero_trip == "floor") 1 else 2.5
  score <- case_when(travel == "0" ~ zero_score, travel == "1-4" ~ 2.5, travel == "5-9" ~ 7, travel == "10+" ~ 10) * multiplier
  weight <- replace_na((1 / score) / mean(1 / score, na.rm = TRUE), 1)
  cohort$weight <- if (weighted) weight else rep(1, nrow(cohort))
  cohort
}

# build each sample's cohort
prepare_cohort <- function(sample) {
  cohort <- direct_cohort(frames[[sample]], sample)
  if (str_detect(sample, "hh")) cohort$cluster <- paste("Cluster", cohort$cluster_id)
  if (str_detect(sample, "cbe")) {
    cohort$cluster_id <- cohort$cluster
    cohort$record_ix <- cohort$record_id
  }
  if (str_detect(sample, "vbs")) cohort <- weight_vbs(cohort, sample)
  cohort
}
cohorts <- map(set_names(samples), prepare_cohort)

## Estimators ----
# The primary estimator of each method, with the options the analyses below
# vary. With imputed completions, each completion is estimated and the results
# pooled by Rubin's rules (Barnard-Rubin degrees of freedom).
# pool estimates over imputed completions, with a Wilson interval at the pooled degrees of freedom
pool_imputations <- function(estimate, standard_error, complete_df) {
  m <- length(estimate)
  within <- mean(standard_error^2)
  between <- var(estimate)
  total <- within + (1 + 1 / m) * between
  df <- complete_df
  if (between > 0) {
    lambda <- (1 + 1 / m) * between / total
    df_old <- (m - 1) / lambda^2
    df_observed <- ((complete_df + 1) / (complete_df + 3)) * complete_df * (1 - lambda)
    df <- if (is.infinite(complete_df)) df_old else df_old * df_observed / (df_old + df_observed)
  }
  interval <- wilson_interval(mean(estimate), sqrt(total), z = qt(0.975, df))
  tibble(estimate = mean(estimate), standard_error = sqrt(total), ci_lower = interval[[1]], ci_upper = interval[[2]])
}

# HH and CBE: clusters then households or records, equal weights; VBS: respondents as clusters, mobility weights
cluster_design <- function(cohort) svydesign(id = ~ cluster_id + record_ix, weights = ~w, nest = TRUE, data = mutate(cohort, w = 1))
vbs_design <- function(cohort) svydesign(id = ~record_ix, weights = ~weight, data = cohort)

# design-based estimate of each outcome, pooled over completions when given
survey_estimates <- function(cohort, build_design, fields = outcomes, completions = NULL) {
  if (is.null(completions)) {
    design <- build_design(cohort)
    return(map(fields, \(field) survey_proportion(design, field)) |> list_rbind(names_to = "antigen"))
  }
  designs <- map(completions, build_design)
  complete_df <- degf(designs[[1]])
  map(fields, \(field) {
    fits <- map(designs, \(design) survey_proportion(design, field)) |> list_rbind()
    pool_imputations(fits$estimate, fits$standard_error, complete_df)
  }) |>
    list_rbind(names_to = "antigen")
}

# recruitment network; Somalia Phase 2 completed from the tracking sheet
rds_nodes <- function(cohort, sample, fields = outcomes) {
  rds_network(cohort, fields, tracking = if (sample == "som_phase2_rds") tracking) |>
    as.data.frame()
}

rds_design <- function(nodes, population) {
  as.rds.data.frame(nodes, id = "id_participant_calc", recruiter.id = "id_referring", network.size = "degree", population.size = population)
}

# successive-sampling (or Volz-Heckathorn) estimate of each outcome; burn_in
# drops respondents recruited before that wave, their recruits becoming the
# start of a chain
rds_estimates <- function(cohort, sample, estimator = "ss", population = population_size[[sample]], burn_in = 0, completions = NULL) {
  nodes <- rds_nodes(cohort, sample)
  design <- rds_design(nodes, population)
  if (burn_in > 0) {
    waves <- get.wave(design)
    nodes <- nodes[waves >= burn_in, ]
    nodes$id_referring[!nodes$id_referring %in% c(nodes$id_participant_calc, "0")] <- "0"
    cohort <- cohort[as.character(cohort$id_participant_calc) %in% nodes$id_participant_calc, ]
    design <- rds_design(nodes, population)
  }
  node_ids <- as.character(design$id_participant_calc)
  seed <- sample_seed(sample)
  map(outcomes, \(field) {
    seed_field <- if (field == "zero_dose") "penta_1p_clean" else field
    fit <- withr::with_seed(derived_seed(seed, paste0("rds_outcome/", seed_field)), {
      if (estimator == "vh") RDS.II.estimates(design, outcome.variable = field, N = population) else RDS.SS.estimates(design, outcome.variable = field, N = population)
    })
    weights <- as.numeric(fit$weights)
    if (is.null(completions)) {
      result <- rds_child_mean(weights * node_children(cohort, node_ids, field), as.numeric(design[[field]]), population)
      interval <- wilson_interval(result$estimate, result$standard_error)
      return(tibble(estimate = result$estimate, standard_error = result$standard_error, ci_lower = interval[[1]], ci_upper = interval[[2]]))
    }
    fits <- map(completions, \(completion) {
      rds_child_mean(weights * node_children(completion, node_ids, field), node_mean(completion, node_ids, field), population)
    })
    pool_imputations(map_dbl(fits, "estimate"), map_dbl(fits, "standard_error"), Inf)
  }) |>
    list_rbind(names_to = "antigen")
}

# the primary estimator of a sample, applied to a cohort
estimate_sample <- function(sample, cohort, completions = NULL, ...) {
  if (str_detect(sample, "hh|cbe")) return(survey_estimates(cohort, cluster_design, completions = completions))
  if (str_detect(sample, "vbs")) return(survey_estimates(cohort, vbs_design, completions = completions))
  rds_estimates(cohort, sample, completions = completions, ...)
}

## Multiple imputation ----
# Missing antigen values imputed five times by predictive mean matching,
# jointly over the six antigen fields, from the respondent and child
# characteristics available in each sample (predictors missing for more than
# half the children, or constant, are dropped). An imputed first dose of 0
# sets the later dose to 0.
imputation_predictors <- c(
  "cluster_id", "cluster", "residence", "has_card", "respondent_age", "respondent_sex", "education", "religion",
  "relation_child", "hh_size", "nsum_child_total", "child_sex", "int_location"
)

# impute the antigen fields of one cohort
impute_antigens <- function(cohort, seed, m = 5) {
  if (!anyNA(cohort[antigen_fields])) return(rep(list(cohort), m))
  predictors <- keep(intersect(imputation_predictors, names(cohort)), \(p) {
    values <- cohort[[p]]
    mean(is.na(values)) <= 0.5 && length(unique(values[!is.na(values)])) > 1
  })
  data <- as.data.frame(cohort)[c(antigen_fields, predictors)]
  data <- mutate(data, across(where(is.character), factor))
  method <- set_names(rep("", ncol(data)), names(data))
  method[antigen_fields] <- "pmm"
  predictor_matrix <- matrix(0L, ncol(data), ncol(data), dimnames = list(names(data), names(data)))
  predictor_matrix[antigen_fields, ] <- 1L
  diag(predictor_matrix) <- 0L
  fit <- withr::with_seed(seed, mice(data, m = m, maxit = 10, method = method, predictorMatrix = predictor_matrix, printFlag = FALSE))
  map(seq_len(m), \(i) {
    completed <- mice::complete(fit, i)
    completion <- cohort
    for (field in antigen_fields) completion[[field]] <- completed[[field]]
    completion$opv_3p_clean[completion$opv_1p_clean == 0] <- 0
    completion$penta_3p_clean[completion$penta_1p_clean == 0] <- 0
    completion$measles_2p_clean[completion$measles_1p_clean == 0] <- 0
    completion$zero_dose <- 1 - completion$penta_1p_clean
    completion
  })
}

# impute and estimate each sample
imputation <- map(samples, \(sample) {
  seed <- derived_seed(sample_seed(sample), paste0(sample_key(sample), "/imputation"))
  estimate_sample(sample, cohorts[[sample]], completions = impute_antigens(cohorts[[sample]], seed)) |>
    mutate(analysis = "imputation", sample = sample, variant = "multiple")
}) |>
  list_rbind()

## RDS estimator and population size ----
# The Volz-Heckathorn estimator in place of successive sampling, and the
# population size halved and doubled.
# estimate each RDS sample under the alternatives
rds_samples <- str_subset(samples, "rds")
rds_alternatives <- map(rds_samples, \(sample) {
  population <- population_size[[sample]]
  bind_rows(
    rds_estimates(cohorts[[sample]], sample, estimator = "vh") |> mutate(analysis = "rds_estimator", variant = "vh"),
    rds_estimates(cohorts[[sample]], sample, population = floor(population * 0.5 + 0.5)) |> mutate(analysis = "rds_population_size", variant = "half"),
    rds_estimates(cohorts[[sample]], sample, population = floor(population * 2 + 0.5)) |> mutate(analysis = "rds_population_size", variant = "double")
  ) |>
    mutate(sample = sample)
}) |>
  list_rbind()

## RDS wave burn-in ----
# Respondents recruited before wave k left out, for k from 1 to 10.
# estimate each RDS sample at each burn-in
burn_in <- map(rds_samples, \(sample) {
  map(1:10, \(k) rds_estimates(cohorts[[sample]], sample, burn_in = k) |> mutate(variant = paste0("wave_", k))) |>
    list_rbind() |>
    mutate(analysis = "rds_burn_in", sample = sample)
}) |>
  list_rbind()

## Interview period ----
# Somalia Phase 2 RDS without the interviews from 30 May, when the team was
# told to probe dose counts further.
# drop the late interviews and re-estimate
late_cohort <- filter(cohorts$som_phase2_rds, as.Date(date) < as.Date("2026-05-30"))
late_period <-
  rds_estimates(late_cohort, "som_phase2_rds") |>
  mutate(analysis = "late_period_cut", sample = "som_phase2_rds", variant = "exclude_from_30_may")

## Card seen or never vaccinated ----
# The household survey restricted to children whose card was seen, plus those
# reported as never vaccinated (strict), or also those vaccinated without a
# card or with an electronic card only; and to card holders alone.
ascertainment_reasons <- list(
  verified = c(
    "Child has never been vaccinated", "Child was vaccinated but was never given a card",
    "An RFID card was presented but no paper card was available."
  ),
  verified_strict = "Child has never been vaccinated",
  card_only = character(0)
)

# restrict each household sample and re-estimate
ascertainment <- map(c("eth_phase1_hh", "som_phase1_hh"), \(sample) {
  cohort <- cohorts[[sample]]
  imap(ascertainment_reasons, \(reasons, variant) {
    keep <- as.character(cohort$has_card) %in% "Yes" | as.character(cohort$reason_no_card) %in% reasons
    estimate_sample(sample, cohort[keep, ]) |> mutate(variant = variant)
  }) |>
    list_rbind() |>
    mutate(analysis = "ascertainment", sample = sample)
}) |>
  list_rbind()

## Card-documented coverage ----
# Bounds on each estimate from card documentation alone: the lower bound counts
# only doses copied from a card, the upper bound counts every child except those
# whose card shows the dose was not given. Estimated over the children with an
# observed value, within the full design.
card_samples <- c("eth_phase1_hh", "eth_phase1_cbe", "som_phase1_hh", "eth_phase2_cbe", "som_phase2_cbe")
card_spec <- tibble(
  antigen = c("Penta3", "OPV3", "MCV1", "Zero-dose"),
  field = c("penta_3p_clean", "opv_3p_clean", "measles_1p_clean", "penta_1p_clean"),
  family = c("penta", "opv", "measles", "penta"),
  reverse = c(FALSE, FALSE, FALSE, TRUE)
)

# derive the bound columns and estimate them
card_documented <- map(card_samples, \(sample) {
  cohort <- cohorts[[sample]]
  map(c("documented_lower", "documented_upper"), \(variant) {
    pmap(card_spec, \(antigen, field, family, reverse) {
      y <- cohort[[field]]
      from_card <- as.character(cohort[[paste0(family, "_source")]]) %in% "Vaccination card"
      card_yes <- from_card & y %in% 1
      card_no <- from_card & y %in% 0
      value <- if (variant == "documented_lower") as.numeric(if (reverse) card_no else card_yes) else as.numeric(!(if (reverse) card_yes else card_no))
      cohort$bound <- if_else(is.na(y), NA_real_, value)
      result <- survey_estimates(cohort, cluster_design, fields = c(bound = "bound"))
      result$antigen <- antigen
      result
    }) |>
      list_rbind() |>
      mutate(variant = variant)
  }) |>
    list_rbind() |>
    mutate(analysis = "card_documented", sample = sample)
}) |>
  list_rbind()

# third-dose completion among children with a first pentavalent dose, all of
# them and those whose pentavalent record came from a card, and the share of
# children with a recorded Penta3 status whose record came from a card
card_completion <- map(card_samples, \(sample) {
  cohort <- cohorts[[sample]]
  from_card <- as.character(cohort$penta_source) %in% "Vaccination card"
  starters <- cohort$penta_1p_clean %in% 1 & !is.na(cohort$penta_3p_clean)
  observed <- !is.na(cohort$penta_3p_clean)
  tibble(
    sample = sample,
    card_completion = mean(cohort$penta_3p_clean[starters & from_card]), n_card_starters = sum(starters & from_card),
    all_completion = mean(cohort$penta_3p_clean[starters]), n_starters = sum(starters),
    card_share = mean(from_card[observed])
  )
}) |>
  list_rbind()

## Photo-corrected household survey ----
# The recorded status of a card-sourced child replaced by what the card
# photographs show, where the photographs give a definite answer for Penta1,
# Penta3, OPV3 or MCV1; Penta3 = 1 implies Penta1 = 1.
photo_columns <- c(penta1 = "penta_1p_clean", penta3 = "penta_3p_clean", opv3 = "opv_3p_clean", mcv1 = "measles_1p_clean")
photo_sources <- c(penta1 = "penta_source", penta3 = "penta_source", opv3 = "opv_source", mcv1 = "measles_source")

# apply the photo status to one household cohort
photo_correct <- function(cohort, country) {
  status <- card_photo_status[card_photo_status$country == country, ]
  at <- match(cohort$record_ix_child, status$record_ix_child)
  changed <- list()
  for (antigen in names(photo_columns)) {
    recorded <- cohort[[photo_columns[[antigen]]]]
    photo <- status[[antigen]][at]
    from_card <- as.character(cohort[[photo_sources[[antigen]]]]) %in% "Vaccination card"
    change <- from_card & !is.na(recorded) & photo %in% c("yes", "no") & (photo == "yes") != (recorded == 1)
    changed[[antigen]] <- change %in% TRUE
    cohort[[photo_columns[[antigen]]]][changed[[antigen]]] <- as.numeric(photo[changed[[antigen]]] == "yes")
  }
  implied_first <- changed$penta3 & cohort$penta_3p_clean %in% 1 & cohort$penta_1p_clean %in% 0
  implied_third <- changed$penta1 & !implied_first & cohort$penta_1p_clean %in% 0 & cohort$penta_3p_clean %in% 1
  cohort$penta_1p_clean[implied_first] <- 1
  cohort$penta_3p_clean[implied_third] <- 0
  cohort$zero_dose <- 1 - cohort$penta_1p_clean
  cohort
}

# correct and re-estimate each household sample, with the corrected counts
photo_cohorts <- list(eth_phase1_hh = photo_correct(cohorts$eth_phase1_hh, "eth"), som_phase1_hh = photo_correct(cohorts$som_phase1_hh, "som"))
photo_corrected <- imap(photo_cohorts, \(cohort, sample) {
  estimate_sample(sample, cohort) |> mutate(analysis = "photo_corrected", sample = sample, variant = "photo_corrected")
}) |>
  list_rbind()
photo_counts <- imap(photo_cohorts, \(cohort, sample) {
  map(primary_outcomes, \(antigen) {
    values <- cohort[[antigen_outcomes[[antigen]]]]
    tibble(sample = sample, antigen = antigen, n_domain = sum(!is.na(values)), n_events = sum(values == 1, na.rm = TRUE))
  }) |>
    list_rbind()
}) |>
  list_rbind()

# bounds on Penta3 and OPV3: every card-sourced child with a card photograph
# whose photographs could not show the third dose counted as unvaccinated
# (lower) and as vaccinated (upper)
photo_bound <- function(cohort, country, bound) {
  status <- card_photo_status[card_photo_status$country == country, ]
  at <- match(cohort$record_ix_child, status$record_ix_child)
  has_card_photo <- (status$n_card_photos[at] >= 1) %in% TRUE
  for (antigen in c("penta3", "opv3")) {
    field <- photo_columns[[antigen]]
    from_card <- as.character(cohort[[photo_sources[[antigen]]]]) %in% "Vaccination card"
    value <- cohort[[field]]
    unverifiable <- from_card & has_card_photo & !is.na(value) & status[[antigen]][at] %in% "unknown"
    cohort[[field]][unverifiable] <- if (bound == "upper") 1 else 0
  }
  cohort
}

# estimate both bounds in each household sample
photo_bounds <- imap(photo_cohorts, \(cohort, sample) {
  country <- str_sub(sample, 1, 3)
  map(c("lower", "upper"), \(bound) estimate_sample(sample, photo_bound(cohort, country, bound)) |> mutate(bound = bound)) |>
    list_rbind() |>
    filter(antigen %in% c("Penta3", "OPV3")) |>
    mutate(sample = sample)
}) |>
  list_rbind()

## Unknown dose counts ----
# Children reported vaccinated with an unknown number of doses, or of unknown
# vaccination status, counted as below the target dose (0) and then as at it
# (1), for Penta3 and OPV3.
dose_targets <- list(
  penta_3p_clean = list(ever = "penta_ever", n = "penta_n", known = c("1", "2", "3")),
  opv_3p_clean = list(ever = "opv_ever", n = "opv_n", known = c("1", "2", "3", "4+"))
)

# children of unknown vaccination status or dose count for one target
unknown_children <- function(cohort, field) {
  spec <- dose_targets[[field]]
  ever <- as.character(cohort[[spec$ever]])
  doses <- as.character(cohort[[spec$n]])
  missing <- is.na(cohort[[field]])
  unknown_status <- missing & !is.na(ever) & !ever %in% c("Yes", "No")
  unknown_doses <- missing & ever %in% "Yes" & !doses %in% spec$known
  unknown_status | unknown_doses
}

# recode the unknown children of one cohort
recode_unknown <- function(cohort, value) {
  for (field in names(dose_targets)) {
    cohort[[field]][unknown_children(cohort, field)] <- value
  }
  cohort
}

# children unknown for Penta3 or OPV3 in each sample
unknown_counts <- map(samples, \(sample) {
  cohort <- cohorts[[sample]]
  tibble(sample = sample, n_unknown = sum(unknown_children(cohort, "penta_3p_clean") | unknown_children(cohort, "opv_3p_clean")))
}) |>
  list_rbind()

# estimate each sample under both bounds
dose_unknown <- map(samples, \(sample) {
  bind_rows(
    estimate_sample(sample, recode_unknown(cohorts[[sample]], 0)) |> mutate(variant = "dose_unknown_below_target"),
    estimate_sample(sample, recode_unknown(cohorts[[sample]], 1)) |> mutate(variant = "dose_unknown_at_target")
  ) |>
    mutate(analysis = "dose_unknown_bounds", sample = sample)
}) |>
  list_rbind()

## VBS weights and venues ----
# Zero-trip respondents scored at the largest weight (floor) or left out, the
# weights dropped altogether, and each Phase 1 venue town left out in turn.
vbs_samples <- str_subset(samples, "vbs")
analytic_vbs <- map(set_names(vbs_samples), \(sample) direct_cohort(frames[[sample]], sample))

# re-weight each VBS sample
vbs_weighting <- map(vbs_samples, \(sample) {
  cohort <- analytic_vbs[[sample]]
  bind_rows(
    estimate_sample(sample, weight_vbs(cohort, sample, zero_trip = "floor")) |> mutate(analysis = "never_travels", variant = "floor"),
    estimate_sample(sample, weight_vbs(cohort, sample, zero_trip = "exclude")) |> mutate(analysis = "never_travels", variant = "exclude"),
    estimate_sample(sample, weight_vbs(cohort, sample, weighted = FALSE)) |> mutate(analysis = "mobility_weighting", variant = "unweighted")
  ) |>
    mutate(sample = sample)
}) |>
  list_rbind()

# children in each re-weighted analysis
vbs_children <- map(vbs_samples, \(sample) {
  cohort <- analytic_vbs[[sample]]
  n_excluded <- nrow(weight_vbs(cohort, sample, zero_trip = "exclude"))
  tibble(sample = sample, variant = c("weighted", "unweighted", "exclude", "floor"), n_children = c(rep(nrow(cohort), 2), n_excluded, nrow(cohort)))
}) |>
  list_rbind()

# leave out one venue town at a time
site_exclusions <- tribble(
  ~sample, ~town, ~analysis, ~variant,
  "som_phase1_vbs", "Adado", "site_inclusion", "exclude_adado",
  "eth_phase1_vbs", "Worabe", "site_inclusion_eth", "exclude_worabe",
  "eth_phase1_vbs", "Inseno", "site_inclusion_eth", "exclude_inseno"
)
site_inclusion <- pmap(site_exclusions, \(sample, town, analysis, variant) {
  cohort <- filter(analytic_vbs[[sample]], !str_detect(int_location, town))
  estimate_sample(sample, weight_vbs(cohort, sample)) |> mutate(analysis = analysis, sample = sample, variant = variant, n_children = nrow(cohort))
}) |>
  list_rbind()

## Recruitment-tree standard errors ----
# Phase 1 RDS standard errors inflated by the ratio of a recruitment-tree
# bootstrap standard error (2,000 resamples) to the empirical-likelihood
# standard error, both for the caregiver-level indicator that all the
# caregiver's children have the outcome; ratios below 1 are not applied.
tree_keys <- c(Penta3 = "penta_3", OPV3 = "opv_3", MCV1 = "measles_1", "Zero-dose" = "zero_dose")

# bootstrap each Phase 1 RDS sample and outcome
tree_ratios <- map(c("eth_phase1_rds", "som_phase1_rds"), \(sample) {
  nodes <- rds_nodes(cohorts[[sample]], sample)
  population <- population_size[[sample]]
  map(primary_outcomes, \(antigen) {
    field <- outcomes[[antigen]]
    seed_field <- if (field == "zero_dose") "penta_1p_clean" else field
    tree_nodes <- nodes
    tree_nodes[[field]] <- if_else(is.na(nodes[[field]]), NA_real_, as.numeric(nodes[[field]] >= 1 - 1e-12))
    design <- rds_design(tree_nodes, population)
    el <- withr::with_seed(derived_seed(sample_seed(sample), paste0("rds_outcome/", seed_field)), RDS.SS.estimates(design, outcome.variable = field, N = population))
    boot_seed <- derived_seed(study_seed, paste0("rds_se_bootstrap/", sample_key(sample), "/", tree_keys[[antigen]], "/recruitment_tree"))
    boot <- withr::with_seed(boot_seed, {
      capture.output(result <- suppressMessages(suppressWarnings(RDS.bootstrap.intervals(
        design, outcome.variable = field, weight.type = "Gile's SS", uncertainty = "Gile", N = population, number.of.bootstrap = 2000
      ))))
      result
    })
    tibble(
      sample = sample, antigen = antigen,
      se_ratio = as.data.frame(boot$interval)[["s.e."]][[1]] / as.data.frame(el$interval)[["s.e."]][[1]]
    )
  }) |>
    list_rbind()
}) |>
  list_rbind()

## Verdicts ----
# Each variant replaces its sample's primary estimates, and the verdicts are
# recomputed: the 24 Phase 1 calibration cells (Newcombe against the
# household survey), the order of the three methods on RMSE (per country and
# pooled), and the 8 Phase 2 agreement cells (inter-method range above 10 pp).
# Unweighted counts are the primary ones, except for the photo-corrected
# household survey, which uses its own.
# count children with an observed outcome and with the outcome
counts <- event_counts(frames)

# order of the three methods on RMSE, per country and pooled
rmse_order <- function(summary) {
  summary |>
    filter(scope == "primary_antigens") |>
    arrange(country, rmse_pp, method) |>
    group_by(country) |>
    summarise(order = paste(method, collapse = " "), .groups = "drop") |>
    deframe()
}

# the three verdict layers for one set of estimates
verdict_layers <- function(estimates, counts) {
  performance <- phase1_performance(estimates, counts)
  ranges <- intermethod_range(estimates)
  list(
    calibration = set_names(performance$not_different, paste(performance$country, performance$method, performance$antigen, sep = "/")),
    rank = rmse_order(phase1_summary(performance)),
    agreement = set_names(ranges$range_pp > 10, paste(ranges$country, ranges$antigen, sep = "/"))
  )
}
baseline <- verdict_layers(estimates, counts)

# swap a variant's estimates into the primary set and list what changes
compare_variant <- function(rows, counts) {
  swapped <-
    estimates |>
    rows_update(
      rows |> mutate(antigen = factor(antigen, levels = levels(estimates$antigen))) |> select(sample, antigen, estimate, standard_error, ci_lower, ci_upper),
      by = c("sample", "antigen")
    )
  layers <- verdict_layers(swapped, counts)
  list(
    calibration = names(baseline$calibration)[layers$calibration != baseline$calibration],
    agreement = names(baseline$agreement)[layers$agreement != baseline$agreement],
    rank = names(baseline$rank)[layers$rank[names(baseline$rank)] != baseline$rank]
  )
}

# every variant, with its change from the primary estimate
variants <-
  bind_rows(
    imputation, rds_alternatives, burn_in, late_period, ascertainment, card_documented, photo_corrected,
    dose_unknown, vbs_weighting, site_inclusion
  ) |>
  left_join(
    estimates |> transmute(sample, antigen = as.character(antigen), primary_estimate = estimate),
    by = c("sample", "antigen")
  ) |>
  mutate(delta_pp = 100 * (estimate - primary_estimate))

# score each variant
photo_variant_counts <- rows_update(counts, photo_counts, by = c("sample", "antigen"))
variant_verdicts <-
  variants |>
  group_by(analysis, sample, variant) |>
  group_modify(\(rows, key) {
    changed <- compare_variant(mutate(rows, sample = key$sample), if (key$analysis == "photo_corrected") photo_variant_counts else counts)
    top <- which.max(abs(rows$delta_pp))
    tibble(
      max_abs_delta_pp = abs(rows$delta_pp[[top]]), max_delta_antigen = rows$antigen[[top]],
      calibration = list(changed$calibration), agreement = list(changed$agreement), rank = list(changed$rank)
    )
  }) |>
  ungroup()

## Recruitment-tree verdicts ----
# The Phase 1 RDS calibration verdicts recomputed with the inflated standard
# error (Wilson interval at 1.96; count-based where the arm is a boundary arm).
# recompute the RDS calibration cells with the inflated standard error
rds_cells <-
  phase1_performance(estimates, counts) |>
  filter(method == "RDS") |>
  mutate(antigen = as.character(antigen), sample = paste0(country, "_phase1_rds")) |>
  left_join(tree_ratios, by = c("sample", "antigen"))
inflated_se <- pmax(1, rds_cells$se_ratio) * rds_cells$standard_error
inflated <- map2(rds_cells$estimate, inflated_se, wilson_interval)
arm <- boundary_arm(rds_cells$estimate, inflated_se, map_dbl(inflated, 1), map_dbl(inflated, 2), rds_cells$n_domain, rds_cells$n_events)
ref <- boundary_arm(rds_cells$ref_estimate, rds_cells$ref_se, rds_cells$ref_lower, rds_cells$ref_upper, rds_cells$ref_n_domain, rds_cells$ref_n_events)
tree_interval <- newcombe_interval(rds_cells$estimate, arm$ci_lower, arm$ci_upper, rds_cells$ref_estimate, ref$ci_lower, ref$ci_upper)
tree_changed <- sum((tree_interval$lower <= 0 & tree_interval$upper >= 0) != rds_cells$not_different)

## Tipping multiplier ----
# The multiplier of the RDS standard error at which a distinguishable Phase 1
# RDS cell would first become not distinguishable (none for a cell already not
# distinguishable or with a count-based interval).
# calibration interval of one cell at multiplier m
calibration_at <- function(cell, m) {
  se <- m * cell$standard_error
  interval <- wilson_interval(cell$estimate, se)
  arm <- boundary_arm(cell$estimate, se, interval[[1]], interval[[2]], cell$n_domain, cell$n_events)
  ref <- boundary_arm(cell$ref_estimate, cell$ref_se, cell$ref_lower, cell$ref_upper, cell$ref_n_domain, cell$ref_n_events)
  c(newcombe_interval(cell$estimate, arm$ci_lower, arm$ci_upper, cell$ref_estimate, ref$ci_lower, ref$ci_upper), boundary = arm$flag)
}

# solve for the multiplier at which the interval first reaches zero
tipping_multiplier <- function(cell) {
  at_one <- calibration_at(cell, 1)
  if (at_one$boundary || (at_one$lower <= 0 && at_one$upper >= 0)) return(NA_real_)
  bound_at <- \(m) {
    interval <- calibration_at(cell, m)
    if (cell$estimate >= cell$ref_estimate) interval$lower else -interval$upper
  }
  high <- 4
  while (bound_at(high) > 0) high <- high * 2
  uniroot(bound_at, interval = c(1, high), tol = 1e-10)$root
}
rds_cells$tipping <- map_dbl(seq_len(nrow(rds_cells)), \(i) tipping_multiplier(rds_cells[i, ]))
rds_cells$tree_not_different <- tree_interval$lower <= 0 & tree_interval$upper >= 0

## Adjustment to reference margins ----
# The Phase 1 CBE, RDS and VBS samples adjusted to the household survey of
# their country, and the Phase 2 RDS and VBS samples to the CBE sample of their
# country, on the caregiver's residence, age band, education, asset ownership
# and (Ethiopia) religion: by raking, with the raking factors bounded to
# [0.2, 5], and by regression standardisation, a quasi-binomial model of the
# outcome predicted onto the reference children.
# the samples and their references, the bounds on the raking factors, the
# order in which characteristics are dropped when the bounds cannot be met,
# and the characteristics of each country
adjust_samples <- tribble(
  ~sample, ~reference,
  "eth_phase1_cbe", "eth_phase1_hh", "eth_phase1_rds", "eth_phase1_hh", "eth_phase1_vbs", "eth_phase1_hh",
  "som_phase1_cbe", "som_phase1_hh", "som_phase1_rds", "som_phase1_hh", "som_phase1_vbs", "som_phase1_hh",
  "eth_phase2_rds", "eth_phase2_cbe", "eth_phase2_vbs", "eth_phase2_cbe",
  "som_phase2_rds", "som_phase2_cbe", "som_phase2_vbs", "som_phase2_cbe"
)
adjust_bounds <- c(0.2, 5)
drop_order <- c("age_band", "religion_grp", "wealth_grp", "education_stage", "residence")
adjust_covariates <- function(country) c("residence", "age_band", "education_stage", "wealth_grp", if (country == "eth") "religion_grp")
bootstrap_reps <- 500

# child cohort of a sample with its base weight (per outcome for RDS: the
# successive-sampling weight of the child's respondent)
adjust_cohort <- function(sample) {
  cohort <- direct_cohort(frames[[sample]], sample)
  if (str_detect(sample, "cbe")) {
    cohort$cluster_id <- cohort$cluster
    cohort$record_ix <- cohort$record_id
  }
  cohort$w_base <- if (str_detect(sample, "vbs")) vbs_weight(cohort, sample) else rep(1, nrow(cohort))
  weights <- NULL
  if (str_detect(sample, "rds")) {
    nodes <- rds_nodes(cohort, sample)
    population <- population_size[[sample]]
    design <- rds_design(nodes, population)
    weights <- map(outcomes, \(field) {
      seed_field <- if (field == "zero_dose") "penta_1p_clean" else field
      fit <- withr::with_seed(derived_seed(sample_seed(sample), paste0("rds_outcome/", seed_field)), RDS.SS.estimates(design, outcome.variable = field, N = population))
      as.numeric(fit$weights)[match(as.character(cohort$id_participant_calc), as.character(design$id_participant_calc))]
    })
  }
  list(cohort = cohort, weights = weights)
}
adjust_cohorts <- map(set_names(unique(c(adjust_samples$sample, adjust_samples$reference))), adjust_cohort)

# caregiver characteristics, with asset ownership scored on the assets asked
# about in every sample of the country and phase
religion_groups <- c("Muslim" = "Muslim", "Orthodox" = "Christian", "Protestant" = "Christian", "Catholic" = "Christian")
caregiver_field <- \(sample) if (str_detect(sample, "rds")) "id_participant_calc" else "record_ix"

# derive the characteristics of each caregiver and attach them to the children
for (group in unique(str_sub(names(adjust_cohorts), 1, 10))) {
  members <- str_subset(names(adjust_cohorts), paste0("^", group))
  respondents <- map(set_names(members), \(sample) {
    cohort <- adjust_cohorts[[sample]]$cohort
    cohort[!duplicated(cohort[[caregiver_field(sample)]]), ]
  })
  assets <- intersect(asset_candidates, reduce(map(respondents, names), intersect))
  for (sample in members) {
    d <- respondents[[sample]]
    age <- d$respondent_age
    d$age_band <- case_when(age <= 24 ~ "under 25", age <= 34 ~ "25-34", age >= 35 ~ "35+")
    d$education_stage <- if ("education" %in% names(d)) unname(education_stages[as.character(d$education)]) else NA_character_
    d$religion_grp <- if ("religion" %in% names(d)) unname(religion_groups[as.character(d$religion)]) else NA_character_
    d$wealth_grp <- wealth_group(d, assets)
    at <- match(adjust_cohorts[[sample]]$cohort[[caregiver_field(sample)]], d[[caregiver_field(sample)]])
    for (field in c("age_band", "education_stage", "religion_grp", "wealth_grp")) adjust_cohorts[[sample]]$cohort[[field]] <- d[[field]][at]
    adjust_cohorts[[sample]]$cohort$residence <- as.character(adjust_cohorts[[sample]]$cohort$residence)
  }
}

# weighted mean over the children with an observed outcome
weighted_share <- function(y, w) {
  keep <- !is.na(y) & !is.na(w)
  if (!any(keep)) return(NA_real_)
  sum(y[keep] * w[keep]) / sum(w[keep])
}

# design of a sample for the regression model
adjust_design <- function(sample, data) {
  if (str_detect(sample, "cbe")) return(svydesign(id = ~ cluster_id + record_ix, weights = ~w_base, nest = TRUE, data = data))
  if (str_detect(sample, "vbs")) return(svydesign(id = ~record_ix, weights = ~w_base, data = data))
  svydesign(ids = ~1, weights = ~w_base, data = data)
}

## Raking ----
# Raking by iterative proportional fitting, accepted when it meets the margins
# with every factor inside the bounds; otherwise the bounded raking problem
# solved through its dual. A characteristic with one level in the sample and
# reference combined is left out.
# bounded raking by its dual: the raking distance continued linearly outside the bounds
bounded_raking <- function(x, w, totals, bounds) {
  low <- log(bounds[1])
  high <- log(bounds[2])
  phi <- \(u) if_else(u < low, exp(low) * (u - low) + exp(low), if_else(u > high, exp(high) * (u - high) + exp(high), exp(u)))
  factor_at <- \(u) pmin(pmax(exp(u), exp(low)), exp(high))
  objective <- \(lambda) sum(w * phi(as.numeric(x %*% lambda))) - sum(lambda * totals)
  gradient <- \(lambda) as.numeric(t(x) %*% (w * factor_at(as.numeric(x %*% lambda)))) - totals
  best <- NULL
  for (i in 1:5) {
    start <- if (is.null(best)) rep(0, ncol(x)) else best$par
    best <- optim(start, objective, gradient, method = "BFGS", control = list(maxit = 5000, reltol = 1e-16))
  }
  factor_at(as.numeric(x %*% best$par))
}

# calibrate a sample to the reference margins on one set of characteristics;
# children with a missing characteristic keep a factor of 1
calibrate_set <- function(sample_df, reference_df, covariates, bounds = adjust_bounds) {
  n <- nrow(sample_df)
  failed <- \(complete, dropped = character(0)) list(ok = FALSE, g = rep(1, n), complete = complete, dropped = dropped)
  w_base <- sample_df$w_base / mean(sample_df$w_base)
  if (length(covariates) == 0) return(failed(rep(FALSE, n)))
  complete_sample <- complete.cases(sample_df[covariates])
  complete_reference <- complete.cases(reference_df[covariates])
  if (sum(complete_sample) == 0 || sum(complete_reference) == 0) return(failed(complete_sample))
  sample_c <- sample_df[complete_sample, , drop = FALSE]
  reference_c <- reference_df[complete_reference, , drop = FALSE]
  levels_of <- map(set_names(covariates), \(cv) sort(unique(c(as.character(sample_c[[cv]]), as.character(reference_c[[cv]])))))
  degenerate <- covariates[lengths(levels_of) < 2]
  covariates <- setdiff(covariates, degenerate)
  if (length(covariates) == 0) return(failed(complete_sample, degenerate))
  for (cv in covariates) {
    sample_c[[cv]] <- factor(as.character(sample_c[[cv]]), levels = levels_of[[cv]])
    reference_c[[cv]] <- factor(as.character(reference_c[[cv]]), levels = levels_of[[cv]])
  }
  sample_c$w_base <- w_base[complete_sample]
  total_mass <- sum(sample_c$w_base)
  reference_c$.ref_w <- reference_c$w_base_ref / sum(reference_c$w_base_ref) * total_mass
  formula <- reformulate(covariates)
  x_reference <- model.matrix(formula, data = reference_c)
  totals <- as.numeric(t(x_reference) %*% reference_c$.ref_w)
  x_sample <- model.matrix(formula, data = sample_c)
  margin_error <- \(w) max(abs(as.numeric(t(x_sample) %*% w) - totals)) / total_mass

  calibrated <- NULL
  raked <- tryCatch(
    withCallingHandlers(
      {
        capture.output(result <- rake(
          svydesign(ids = ~1, weights = ~w_base, data = sample_c),
          sample.margins = map(covariates, \(cv) reformulate(cv)),
          population.margins = map(covariates, \(cv) xtabs(as.formula(paste0(".ref_w ~ ", cv)), data = reference_c)),
          control = list(maxit = 500, epsilon = 1e-9, verbose = FALSE)
        ))
        result
      },
      warning = \(w) invokeRestart("muffleWarning")
    ),
    error = \(e) NULL
  )
  if (!is.null(raked)) {
    w_raked <- weights(raked)
    g_raked <- w_raked / sample_c$w_base
    if (is.finite(margin_error(w_raked)) && margin_error(w_raked) <= 1e-6 && all(g_raked >= bounds[1] - 1e-9 & g_raked <= bounds[2] + 1e-9)) calibrated <- w_raked
  }
  if (is.null(calibrated)) {
    g_dual <- tryCatch(bounded_raking(x_sample, sample_c$w_base, totals, bounds), error = \(e) NULL)
    if (!is.null(g_dual)) {
      w_dual <- sample_c$w_base * g_dual
      if (is.finite(margin_error(w_dual)) && margin_error(w_dual) <= 1e-6) calibrated <- w_dual
    }
  }
  if (is.null(calibrated)) return(failed(complete_sample, degenerate))
  g <- rep(1, n)
  g[complete_sample] <- calibrated / sample_c$w_base
  list(ok = TRUE, g = g, complete = complete_sample, dropped = degenerate)
}

# calibrate on as many characteristics as the bounds allow: drop them in a
# fixed order until a calibration succeeds, then add back any that fit
calibrate <- function(sample_df, reference_df, covariates) {
  remaining <- covariates
  dropped <- character(0)
  result <- calibrate_set(sample_df, reference_df, remaining)
  while (!result$ok && length(remaining) > 0) {
    next_drop <- drop_order[drop_order %in% remaining][1]
    dropped <- c(dropped, next_drop)
    remaining <- setdiff(remaining, next_drop)
    result <- calibrate_set(sample_df, reference_df, remaining)
  }
  changed <- TRUE
  while (changed && length(dropped) > 0) {
    changed <- FALSE
    for (cv in rev(dropped)) {
      trial <- calibrate_set(sample_df, reference_df, union(remaining, cv))
      if (trial$ok) {
        remaining <- union(remaining, cv)
        dropped <- setdiff(dropped, cv)
        result <- trial
        changed <- TRUE
      }
    }
  }
  g <- result$g[result$complete]
  list(
    g = result$g,
    covariates_used = if (result$ok) setdiff(remaining, result$dropped) else character(0),
    at_lower = mean(abs(g - adjust_bounds[1]) < 1e-9), at_upper = mean(abs(g - adjust_bounds[2]) < 1e-9)
  )
}

## Regression standardisation ----
# A quasi-binomial model of the outcome on the characteristics, its predictions
# averaged over the reference children; not converged where a reference
# category is absent from the sample.
regression_standardise <- function(design, reference_df, covariates, field) {
  formula <- reformulate(covariates, response = field)
  fit <- tryCatch(
    svyglm(formula, design = design, family = quasibinomial()),
    error = \(e) NULL, warning = \(w) suppressWarnings(svyglm(formula, design = design, family = quasibinomial()))
  )
  if (is.null(fit)) return(list(estimate = NA_real_, converged = FALSE))
  reference_c <- reference_df[complete.cases(reference_df[covariates]), , drop = FALSE]
  for (cv in covariates) reference_c[[cv]] <- factor(as.character(reference_c[[cv]]), levels = fit$xlevels[[cv]])
  unmatched <- reduce(map(covariates, \(cv) is.na(reference_c[[cv]])), `|`, .init = rep(FALSE, nrow(reference_c)))
  reference_c <- reference_c[!unmatched, , drop = FALSE]
  prediction <- as.numeric(suppressWarnings(predict(fit, newdata = reference_c, type = "response")))
  list(estimate = sum(prediction * reference_c$w_base_ref) / sum(reference_c$w_base_ref), converged = isTRUE(fit$converged) && !any(unmatched))
}

## Bootstrap ----
# For CBE and VBS, 500 resamples of the sample's clusters (VBS: respondents)
# and of the reference clusters, each adjusted on the characteristics of the
# point estimate; an interval is given only where at least 90% of the
# resamples produce an estimate. Verdicts on the difference from the reference
# use the percentile interval of the paired differences.
# one resample: adjusted, unadjusted and reference estimates
bootstrap_replicate <- function(sample, reference, field, weight, covariates, analysis) {
  cohort <- adjust_cohorts[[sample]]$cohort
  reference_cohort <- adjust_cohorts[[reference]]$cohort
  unit_field <- if ("cluster_id" %in% names(cohort)) "cluster_id" else "record_ix"
  units <- unique(cohort[[unit_field]])
  rows <- unlist(map(sample(units, length(units), replace = TRUE), \(u) which(cohort[[unit_field]] == u)))
  reference_units <- unique(reference_cohort$cluster_id)
  reference_rows <- unlist(map(sample(reference_units, length(reference_units), replace = TRUE), \(u) which(reference_cohort$cluster_id == u)))
  boot <- cohort[rows, , drop = FALSE]
  boot_reference <- reference_cohort[reference_rows, , drop = FALSE]
  boot_reference$w_base_ref <- boot_reference$w_base

  sample_df <- boot[covariates]
  sample_df[[field]] <- boot[[field]]
  sample_df$w_base <- weight[rows]
  if ("cluster_id" %in% names(boot)) sample_df$cluster_id <- boot$cluster_id
  sample_df$record_ix <- seq_len(nrow(boot))
  reference_df <- boot_reference[c(covariates, "w_base_ref", field)]

  adjusted <- tryCatch({
    if (analysis == "raking") {
      fit <- calibrate_set(sample_df, reference_df, covariates)
      if (fit$ok) weighted_share(sample_df[[field]], sample_df$w_base / mean(sample_df$w_base) * fit$g) else NA_real_
    } else {
      complete <- sample_df[complete.cases(sample_df[c(covariates, field)]), , drop = FALSE]
      if (nrow(complete) == 0) NA_real_ else {
        fit <- regression_standardise(adjust_design(sample, complete), reference_df, covariates, field)
        if (fit$converged) fit$estimate else NA_real_
      }
    }
  }, error = \(e) NA_real_)
  c(adjusted = adjusted, unadjusted = weighted_share(sample_df[[field]], sample_df$w_base), reference = weighted_share(boot_reference[[field]], boot_reference$w_base_ref))
}

# summarise the resamples of one estimate
bootstrap_summary <- function(boot, point) {
  ok <- complete.cases(boot)
  if (sum(ok) < 0.9 * nrow(boot)) return(list(se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, diff_lower = NA_real_, diff_upper = NA_real_))
  se <- sd(boot[ok, "adjusted"])
  interval <- if (se > 0 && point > 0 && point < 1) wilson_interval(point, se, qnorm(0.975)) else c(max(0, point - qnorm(0.975) * se), min(1, point + qnorm(0.975) * se))
  difference <- quantile(boot[ok, "adjusted"] - boot[ok, "reference"], c(0.025, 0.975), names = FALSE)
  list(se = se, ci_lower = interval[[1]], ci_upper = interval[[2]], diff_lower = difference[[1]], diff_upper = difference[[2]])
}

## Adjusted estimates ----
# Each sample and outcome under both adjustments. RDS raking takes its interval
# from the standard error of the RDS estimator at the raked weights (Wilson
# interval); RDS regression standardisation has no interval. The CBE and VBS
# bootstraps run in parallel, one R session per core but one, each resampled
# estimate from its own seed.
# point estimates of one sample on every outcome
adjusted_points <- pmap(adjust_samples, \(sample, reference) {
  country <- str_sub(sample, 1, 3)
  covariates <- adjust_covariates(country)
  entry <- adjust_cohorts[[sample]]
  reference_df <- adjust_cohorts[[reference]]$cohort |> mutate(w_base_ref = w_base)
  reference_df <- as.data.frame(reference_df)[c(covariates, "w_base_ref")]
  map(names(outcomes), \(antigen) {
    field <- outcomes[[antigen]]
    sample_df <- as.data.frame(entry$cohort)[unique(c(intersect(c("cluster_id", "record_ix", "id_participant_calc"), names(entry$cohort)), covariates, field))]
    sample_df$w_base <- if (is.null(entry$weights)) entry$cohort$w_base else entry$weights[[antigen]]

    # raking
    cal <- calibrate(sample_df, reference_df, covariates)
    w_raked <- sample_df$w_base / mean(sample_df$w_base) * cal$g
    raked <- weighted_share(sample_df[[field]], w_raked)

    # regression standardisation on the characteristics used by the raking
    used <- if (length(cal$covariates_used) > 0) cal$covariates_used else covariates
    complete <- sample_df[complete.cases(sample_df[c(used, field)]), , drop = FALSE]
    regression <- regression_standardise(adjust_design(sample, complete), reference_df, used, field)$estimate

    # RDS standard error at the raked weights
    rds_se <- NA_real_
    if (str_detect(sample, "rds")) {
      keep <- !is.na(sample_df[[field]]) & !is.na(w_raked)
      rds_se <- rds_child_mean(w_raked[keep], sample_df[[field]][keep], population_size[[sample]])$standard_error
    }
    tibble(
      sample = sample, reference = reference, antigen = antigen, analysis = c("raking", "regression"), estimate = c(raked, regression),
      rds_se = c(rds_se, NA), used = list(used), covariates_used = paste(cal$covariates_used, collapse = "+"), at_lower = cal$at_lower, at_upper = cal$at_upper
    )
  }) |>
    list_rbind()
}) |>
  list_rbind()

# resample one CBE or VBS estimate
bootstrap_cell <- function(sample, reference, antigen, analysis, used) {
  entry <- adjust_cohorts[[sample]]
  seed <- derived_seed(study_seed, paste0("poststrat/", sample_key(sample), "/poststrat_", if (analysis == "raking") "raking" else "regression_std"))
  withr::with_seed(seed, {
    draws <- matrix(NA_real_, bootstrap_reps, 3, dimnames = list(NULL, c("adjusted", "unadjusted", "reference")))
    for (b in seq_len(bootstrap_reps)) {
      draws[b, ] <- tryCatch(bootstrap_replicate(sample, reference, outcomes[[antigen]], entry$cohort$w_base, used, analysis), error = \(e) rep(NA_real_, 3))
    }
    draws
  })
}

# run the bootstraps in parallel sessions
bootstrap_tasks <- adjusted_points |> filter(!str_detect(sample, "rds"))
sessions <- makeCluster(max(1L, detectCores() - 1L))
invisible(clusterEvalQ(sessions, {
  pacman::p_load(survey, tidyverse)
  source("code/functions.R")
}))
clusterExport(sessions, c(
  "adjust_cohorts", "outcomes", "bootstrap_reps", "adjust_bounds", "sample_key", "weighted_share", "adjust_design",
  "calibrate_set", "bounded_raking", "regression_standardise", "bootstrap_replicate", "bootstrap_cell"
))
bootstrap_draws <- clusterMap(
  sessions, bootstrap_cell, bootstrap_tasks$sample, bootstrap_tasks$reference, bootstrap_tasks$antigen,
  bootstrap_tasks$analysis, bootstrap_tasks$used, .scheduling = "dynamic", SIMPLIFY = FALSE
)
stopCluster(sessions)

# intervals of every adjusted estimate
adjusted <-
  bind_rows(
    bootstrap_tasks |> mutate(summary = map2(bootstrap_draws, estimate, bootstrap_summary)) |> unnest_wider(summary),
    adjusted_points |>
      filter(str_detect(sample, "rds")) |>
      mutate(
        se = rds_se,
        ci_lower = map2_dbl(estimate, rds_se, \(p, se) if (is.na(se)) NA_real_ else wilson_interval(p, se)[[1]]),
        ci_upper = map2_dbl(estimate, rds_se, \(p, se) if (is.na(se)) NA_real_ else wilson_interval(p, se)[[2]])
      )
  ) |>
  select(-rds_se, -used) |>
  left_join(sample_parts(adjust_samples$sample), by = "sample")

## Accuracy and agreement after adjustment ----
# Phase 1: each adjusted estimate against the household survey's primary
# estimate. Phase 2: the inter-method range with RDS and VBS adjusted, and the
# RDS minus CBE difference with Newcombe's interval where RDS has one.
# verdict of each Phase 1 cell
adjusted_phase1 <-
  adjusted |>
  filter(phase == 1) |>
  left_join(
    estimates |> filter(method == "HH") |> transmute(country, antigen = as.character(antigen), ref_estimate = estimate, ref_se = standard_error, ref_lower = ci_lower, ref_upper = ci_upper, ref_sample = sample),
    by = c("country", "antigen")
  ) |>
  left_join(counts, by = c("sample", "antigen")) |>
  left_join(counts |> rename(ref_sample = sample, ref_n_domain = n_domain, ref_n_events = n_events), by = c("ref_sample", "antigen")) |>
  mutate(difference_pp = 100 * (estimate - ref_estimate), has_interval = !is.na(ci_lower))
adjusted_phase1$not_different <- map_lgl(seq_len(nrow(adjusted_phase1)), \(i) {
  cell <- adjusted_phase1[i, ]
  if (!cell$has_interval) return(NA)
  if (cell$method != "rds") return(cell$diff_lower <= 0 && cell$diff_upper >= 0)
  arm <- boundary_arm(cell$estimate, cell$se, cell$ci_lower, cell$ci_upper, cell$n_domain, cell$n_events)
  ref <- boundary_arm(cell$ref_estimate, cell$ref_se, cell$ref_lower, cell$ref_upper, cell$ref_n_domain, cell$ref_n_events)
  interval <- newcombe_interval(cell$estimate, arm$ci_lower, arm$ci_upper, cell$ref_estimate, ref$ci_lower, ref$ci_upper)
  interval$lower <= 0 && interval$upper >= 0
})

# RMSE and verdict counts over the three primary antigens
adjusted_summary <-
  adjusted_phase1 |>
  filter(antigen != "Zero-dose") |>
  summarise(rmse_pp = sqrt(mean(difference_pp^2)), n_intervals = sum(has_interval), n_not_different = sum(not_different, na.rm = TRUE), .by = c(country, method, analysis))
unadjusted_summary <-
  phase1_summary(phase1_performance(estimates, counts)) |>
  filter(scope == "primary_antigens", country != "pooled") |>
  transmute(country = as.character(country), method = str_to_lower(method), rmse_pp, n_not_different, n_cells)

# Phase 2 ranges and RDS minus CBE, unadjusted and adjusted
cbe_phase2 <- estimates |> filter(phase == 2, method == "CBE") |> transmute(country, antigen = as.character(antigen), cbe = estimate, cbe_lower = ci_lower, cbe_upper = ci_upper)
phase2_versions <-
  bind_rows(
    estimates |> filter(phase == 2, method %in% c("RDS", "VBS"), antigen %in% primary_outcomes) |>
      transmute(country, method = str_to_lower(method), antigen = as.character(antigen), analysis = "unadjusted", estimate, ci_lower, ci_upper),
    adjusted |> filter(phase == 2) |> select(country, method, antigen, analysis, estimate, ci_lower, ci_upper)
  ) |>
  left_join(cbe_phase2, by = c("country", "antigen"))
adjusted_phase2 <-
  phase2_versions |>
  summarise(range_pp = 100 * (max(c(estimate, first(cbe))) - min(c(estimate, first(cbe)))), .by = c(country, antigen, analysis)) |>
  left_join(
    phase2_versions |>
      filter(method == "rds") |>
      mutate(
        difference_pp = 100 * (estimate - cbe),
        distinguishable = !((estimate - cbe) - sqrt((estimate - ci_lower)^2 + (cbe_upper - cbe)^2) <= 0 & (estimate - cbe) + sqrt((ci_upper - estimate)^2 + (cbe - cbe_lower)^2) >= 0)
      ) |>
      select(country, antigen, analysis, difference_pp, distinguishable),
    by = c("country", "antigen", "analysis")
  )

# Tables ----

## Table 11 ----
# Somalia Phase 2 RDS: Penta3 and OPV3 with the children of unknown dose count
# excluded (the primary estimate) and multiply imputed.
# complete-case (primary) and multiply imputed Somalia Phase 2 RDS estimates
table_11 <-
  bind_rows(
    estimates |> filter(sample == "som_phase2_rds") |> mutate(antigen = as.character(antigen), treatment = "Excluded (complete case)"),
    imputation |> filter(sample == "som_phase2_rds") |> mutate(treatment = "Multiple imputation")
  ) |>
  filter(antigen %in% c("Penta3", "OPV3")) |>
  mutate(value = format_ci(estimate, ci_lower, ci_upper), antigen = paste(antigen, "(%)")) |>
  select(`Treatment of the unknown children` = treatment, antigen, value) |>
  pivot_wider(names_from = antigen, values_from = value) |>
  select(1, `Penta3 (%)`, `OPV3 (%)`)

save_table(table_11, "output/tables/table_11.csv")

## Table 13 ----
# One row per sensitivity analysis: the samples it was run on, the largest
# change in a primary estimate and where it occurred, and the verdicts that
# change.
# the analyses in the order of Methods Table 6, with the analyses each row gathers
table_13_rows <- tribble(
  ~label, ~analyses,
  "Multiple imputation of missing antigen values", "imputation",
  "Estimator (RDS)", "rds_estimator",
  "Card seen or never vaccinated (household survey, Phase 1)", "ascertainment",
  "Card-documented coverage (household survey, Phase 1; CBE)", "card_documented",
  "Photo-corrected household survey (household survey, Phase 1)", "photo_corrected",
  "RDS population-size input", "rds_population_size",
  "RDS wave burn-in", "rds_burn_in",
  "RDS standard error, recruitment-tree inflation (Phase 1 calibration)", "tree",
  "Unknown dose counts, bounds", "dose_unknown_bounds",
  "Interview period (RDS, Somalia Phase 2)", "late_period_cut",
  "Mobility weights (VBS)", "never_travels;mobility_weighting",
  "Venue site inclusion (VBS, Phase 1)", "site_inclusion;site_inclusion_eth"
)

# wording of the samples, variants and countries in the Where and RMSE columns
variant_words <- c(
  multiple = "multiple imputation", vh = "Volz-Heckathorn estimator", half = "population size halved",
  double = "population size doubled", exclude_from_30_may = "interviews from 30 May excluded",
  verified = "card seen or never vaccinated", verified_strict = "card seen or never vaccinated, strict", card_only = "card only",
  documented_lower = "documentary lower bound", documented_upper = "documentary upper bound", photo_corrected = "photo-corrected",
  dose_unknown_below_target = "unknown children counted as below the target dose",
  dose_unknown_at_target = "unknown children counted as at the target dose",
  floor = "zero-trip respondents at the largest weight", exclude = "zero-trip respondents excluded", unweighted = "unweighted",
  exclude_adado = "Adado excluded", exclude_worabe = "Worabe excluded", exclude_inseno = "Inseno excluded"
)
number_words <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")
describe_variant <- function(variant) {
  if_else(
    str_detect(variant, "^wave_"),
    paste("first", number_words[as.integer(str_extract(variant, "(?<=^wave_)\\d+"))], if_else(variant == "wave_1", "wave", "waves"), "excluded"),
    unname(variant_words[variant])
  )
}
describe_sample <- function(sample) {
  sprintf("%s Phase %s %s", country_labels[str_sub(sample, 1, 3)], str_extract(sample, "(?<=phase)\\d"), str_to_upper(str_extract(sample, "[a-z]+$")))
}
rank_words <- c(eth = "Ethiopia", som = "Somalia", pooled = "pooled")

# summarise one row of the table
summarise_analysis <- function(label, analyses) {
  rows <- filter(variant_verdicts, analysis %in% str_split_1(analyses, ";"))
  top <- rows[which.max(rows$max_abs_delta_pp), ]
  rank <- intersect(names(rank_words), unique(unlist(rows$rank)))
  tibble(
    `Analysis (Methods Table 6)` = label,
    Branches = as.character(n_distinct(rows$sample)),
    `Largest change in a primary estimate (pp)` = format_number(top$max_abs_delta_pp),
    Where = paste(describe_sample(top$sample), describe_variant(top$variant), if_else(top$max_delta_antigen == "Zero-dose", "zero-dose prevalence", top$max_delta_antigen), sep = ", "),
    `Phase 1 calibration verdicts changed` = as.character(n_distinct(unlist(rows$calibration))),
    `Phase 2 agreement verdicts changed` = if (label == table_13_rows$label[[5]]) "–" else as.character(n_distinct(unlist(rows$agreement))),
    `Order of point RMSEs changed` = if (length(rank) > 0) sprintf("Yes (%s)", paste(rank_words[rank], collapse = ", ")) else "No"
  )
}

# the recruitment-tree row: a variance-only change, Phase 1 only
tree_row <- tibble(
  `Analysis (Methods Table 6)` = table_13_rows$label[[8]], Branches = "2",
  `Largest change in a primary estimate (pp)` = "0.0", Where = "–",
  `Phase 1 calibration verdicts changed` = as.character(tree_changed),
  `Phase 2 agreement verdicts changed` = "–", `Order of point RMSEs changed` = "No"
)

# one row per analysis
table_13 <-
  pmap(table_13_rows, \(label, analyses) if (analyses == "tree") tree_row else summarise_analysis(label, analyses)) |>
  list_rbind()

save_table(table_13, "output/tables/table_13.csv")

## Annex Table H1 ----
# Third-dose completion among children with a first pentavalent dose, and the
# documentary brackets for Penta3 and zero-dose prevalence, in the samples that
# recorded where each answer came from.
# documentary bracket of each sample and outcome
brackets <-
  card_documented |>
  select(sample, antigen, variant, estimate) |>
  pivot_wider(names_from = variant, values_from = estimate) |>
  mutate(bracket = paste(format_number(100 * documented_lower), "to", format_number(100 * documented_upper))) |>
  select(sample, antigen, bracket) |>
  pivot_wider(names_from = antigen, values_from = bracket)

# one row per sample, household surveys first
table_H01 <-
  card_completion |>
  left_join(sample_parts(card_completion$sample), by = "sample") |>
  arrange(phase, method != "hh", country) |>
  left_join(brackets, by = "sample") |>
  left_join(estimates |> filter(antigen == "Penta3") |> select(sample, primary = estimate), by = "sample") |>
  transmute(
    Sample = sprintf("%s %s, Phase %s", str_to_upper(str_extract(sample, "[a-z]+$")), country_labels[str_sub(sample, 1, 3)], str_extract(sample, "(?<=phase)\\d")),
    `Third dose, card-documented starters, % (n)` = sprintf("%s (%s)", format_number(100 * card_completion), format_count(n_card_starters)),
    `Third dose, all starters, % (n)` = sprintf("%s (%s)", format_number(100 * all_completion), format_count(n_starters)),
    `Primary Penta3, %` = format_number(100 * primary),
    `Penta3 documentary bracket, %` = Penta3,
    `Share based on a card, %` = format_number(100 * card_share),
    `Zero-dose documentary bracket, %` = `Zero-dose`
  )

save_table(table_H01, "output/tables/table_H01.csv")

## Annex Table H3 ----
# Household-survey coverage as recorded and after the confirmed photograph
# corrections, with the bounds for Penta3 and OPV3.
# bounds as a range of point estimates
bound_ranges <-
  photo_bounds |>
  select(sample, antigen, bound, estimate) |>
  pivot_wider(names_from = bound, values_from = estimate) |>
  mutate(bounds = paste(format_number(100 * lower), "to", format_number(100 * upper)))

# one row per country and outcome
table_H03 <-
  estimates |>
  filter(sample %in% names(photo_cohorts), antigen %in% primary_outcomes) |>
  mutate(antigen = as.character(antigen)) |>
  select(sample, country, antigen, estimate, ci_lower, ci_upper) |>
  left_join(photo_corrected |> select(sample, antigen, corrected = estimate, corrected_lower = ci_lower, corrected_upper = ci_upper), by = c("sample", "antigen")) |>
  left_join(bound_ranges |> select(sample, antigen, bounds), by = c("sample", "antigen")) |>
  arrange(country, factor(antigen, levels = primary_outcomes)) |>
  transmute(
    `Country and antigen` = paste0(country_labels[country], ", ", if_else(antigen == "Zero-dose", "zero-dose", antigen)),
    Recorded = format_ci(estimate, ci_lower, ci_upper),
    `Photo-corrected` = format_ci(corrected, corrected_lower, corrected_upper),
    Bounds = replace_na(bounds, "–")
  )

save_table(table_H03, "output/tables/table_H03.csv")

## Annex Table H4 ----
# Penta3 and OPV3 with the unknown children excluded (the primary estimate),
# counted as below the third dose, and counted as at it, by sample; Somalia
# Phase 2 RDS with its 95% CI.
# the three treatments of each sample and antigen
unknown_values <-
  bind_rows(
    estimates |> filter(antigen %in% c("Penta3", "OPV3")) |> mutate(antigen = as.character(antigen), variant = "complete"),
    dose_unknown |> filter(antigen %in% c("Penta3", "OPV3"))
  ) |>
  filter(sample %in% samples) |>
  mutate(
    value = if_else(sample == "som_phase2_rds", format_ci(estimate, ci_lower, ci_upper), format_number(100 * estimate)),
    column = paste0(antigen, ", ", c(complete = "complete case", dose_unknown_below_target = "below", dose_unknown_at_target = "at")[variant])
  ) |>
  select(sample, column, value) |>
  pivot_wider(names_from = column, values_from = value)

# one row per sample in report order
table_H04 <-
  sample_parts(intersect(study_samples, unknown_values$sample)) |>
  left_join(unknown_counts, by = "sample") |>
  left_join(unknown_values, by = "sample") |>
  transmute(
    `Phase, country` = paste0(phase, ", ", country_labels[country]),
    Method = str_to_upper(method),
    `Unknown children` = format_count(n_unknown),
    `Penta3, complete case`, `Penta3, below`, `Penta3, at`,
    `OPV3, complete case`, `OPV3, below`, `OPV3, at`
  )

save_table(table_H04, "output/tables/table_H04.csv")

## Annex Table H5 ----
# Phase 1 RDS calibration verdicts on the primary antigens at the primary
# standard error and at the recruitment-tree ratio, with the multiplier at
# which each verdict would change.
# one row per country and antigen
table_H05 <-
  rds_cells |>
  filter(antigen != "Zero-dose") |>
  transmute(
    Country = country_labels[country],
    Antigen = antigen,
    `Difference from HH (pp)` = format_signed(difference_pp),
    `Recruitment-tree ratio` = if_else(se_ratio < 1, paste(format_number(se_ratio, 2), "(1 applied)"), format_number(se_ratio, 2)),
    `Verdict, primary standard error` = verdict_label(!not_different),
    `Verdict at the recruitment-tree ratio` = verdict_label(!tree_not_different),
    `Multiplier at which the verdict would change` = format_number(tipping, 2)
  )

save_table(table_H05, "output/tables/table_H05.csv")

## Annex Table H6 ----
# Phase 1 VBS coverage under the alternative treatments of the mobility
# weights and of the venue towns, with the household survey as reference.
# estimates and children of each analysis
vbs_rows <-
  bind_rows(
    estimates |> filter(sample %in% vbs_samples) |> mutate(antigen = as.character(antigen), variant = "weighted"),
    vbs_weighting,
    site_inclusion
  ) |>
  select(sample, antigen, variant, estimate) |>
  left_join(bind_rows(vbs_children, site_inclusion |> distinct(sample, variant, n_children)), by = c("sample", "variant"))
reference_rows <-
  estimates |>
  filter(sample %in% c("eth_phase1_hh", "som_phase1_hh")) |>
  transmute(sample, antigen = as.character(antigen), estimate, variant = "reference", n_children)

# one row per country and analysis
analysis_labels <- c(
  reference = "Household survey (reference)", weighted = "VBS, primary (mobility-weighted)", unweighted = "VBS, unweighted",
  exclude = "VBS, zero-trip caregivers excluded", floor = "VBS, zero-trip caregivers given the largest weight",
  exclude_worabe = "VBS, Worabe excluded", exclude_inseno = "VBS, Inseno excluded", exclude_adado = "VBS, Adado excluded"
)
table_H06 <-
  bind_rows(reference_rows, vbs_rows) |>
  filter(str_detect(sample, "phase1"), antigen %in% primary_outcomes) |>
  mutate(country = str_sub(sample, 1, 3), value = format_number(100 * estimate)) |>
  select(country, variant, n_children, antigen, value) |>
  pivot_wider(names_from = antigen, values_from = value) |>
  arrange(country, factor(variant, levels = names(analysis_labels))) |>
  transmute(
    Country = country_labels[country], Analysis = analysis_labels[variant], Children = format_count(n_children),
    Penta3, OPV3, MCV1, `Zero-dose prevalence` = `Zero-dose`
  )

save_table(table_H06, "output/tables/table_H06.csv")

## Annex Table I1 ----
# Phase 1 RMSE over the primary antigens against the household survey,
# unadjusted and adjusted (with the change from unadjusted), and the primary
# antigens not distinguishable from it; a double dagger where no adjusted
# primary antigen has an interval.
# a signed change with a plus sign when positive
format_change <- function(x) if_else(round_half_up(x) > 0, paste0("+", format_number(x)), format_signed(x))

# adjusted cells across columns
adjusted_cells <-
  adjusted_summary |>
  mutate(
    marker = if_else(n_intervals == 0, "‡", ""),
    rmse = sprintf("%s%s (%s)", format_number(rmse_pp), marker, format_change(rmse_pp - unadjusted_summary$rmse_pp[match(paste(country, method), paste(unadjusted_summary$country, unadjusted_summary$method))])),
    verdict = if_else(n_intervals == 0, "‡", paste(n_not_different, "of", n_intervals))
  ) |>
  select(country, method, analysis, rmse, verdict) |>
  pivot_wider(names_from = analysis, values_from = c(rmse, verdict))

# one row per country and method
table_I01 <-
  unadjusted_summary |>
  filter(method != "hh") |>
  left_join(adjusted_cells, by = c("country", "method")) |>
  arrange(country, method) |>
  transmute(
    Country = country_labels[country], Method = str_to_upper(method),
    `RMSE, unadjusted` = format_number(rmse_pp), `RMSE, raking` = rmse_raking, `RMSE, regression standardisation` = rmse_regression,
    `Not distinguishable, unadjusted` = paste(n_not_different, "of", n_cells),
    `Not distinguishable, raking` = verdict_raking, `Not distinguishable, regression standardisation` = verdict_regression
  )

save_table(table_I01, "output/tables/table_I01.csv")

## Annex Table I2 ----
# The characteristics each calibration used, the share of raking factors at
# each bound, and the outcomes with an adjusted interval, by sample.
# describe a set of characteristics or outcomes
margin_labels <- c(residence = "residence", age_band = "caregiver age", education_stage = "caregiver education", wealth_grp = "asset ownership", religion_grp = "religion")
describe_margins <- function(used, country) {
  all <- adjust_covariates(country)
  used <- str_split_1(used, "\\+")
  if (setequal(used, all)) return(paste("All", if (length(all) == 5) "five" else "four"))
  str_to_sentence(paste(margin_labels[intersect(all, used)], collapse = ", "))
}
describe_intervals <- function(antigen, has_interval) {
  if (all(has_interval)) return("All")
  if (!any(has_interval)) return("None")
  if (sum(!has_interval) == 1) return(paste("All but", if_else(antigen[!has_interval] == "Zero-dose", "zero-dose", antigen[!has_interval])))
  paste(outcome_label(antigen[has_interval]), collapse = ", ")
}

# one row per sample
table_I02 <-
  adjusted |>
  summarise(
    used = first(covariates_used), at_lower = first(at_lower), at_upper = first(at_upper),
    raking = describe_intervals(antigen[analysis == "raking"], !is.na(ci_lower[analysis == "raking"])),
    regression = describe_intervals(antigen[analysis == "regression"], !is.na(ci_lower[analysis == "regression"])),
    .by = c(sample, country, phase, method)
  ) |>
  arrange(phase, country, method) |>
  mutate(first = !duplicated(paste(phase, country))) |>
  transmute(
    `Country, phase` = if_else(first, paste0(country_labels[country], ", ", phase), ""), Method = str_to_upper(method),
    `Margins used` = map2_chr(used, country, describe_margins),
    `Raking factors at the lower bound (%)` = format_number(100 * at_lower), `Raking factors at the upper bound (%)` = format_number(100 * at_upper),
    `Intervals, raking` = raking, `Intervals, regression standardisation` = regression
  )

save_table(table_I02, "output/tables/table_I02.csv")

## Annex Table I3 ----
# Phase 2 inter-method range and RDS minus CBE, unadjusted and with RDS and
# VBS adjusted to the CBE margins; an asterisk marks RDS distinguishable from
# CBE (no interval, so no verdict, under regression standardisation).
# one row per country and outcome
table_I03 <-
  adjusted_phase2 |>
  mutate(
    range = format_number(range_pp),
    difference = paste0(format_change(difference_pp), if_else(distinguishable %in% TRUE, "*", ""))
  ) |>
  select(country, antigen, analysis, range, difference) |>
  pivot_wider(names_from = analysis, values_from = c(range, difference)) |>
  arrange(country, factor(antigen, levels = primary_outcomes)) |>
  mutate(first = !duplicated(country)) |>
  transmute(
    Country = if_else(first, country_labels[country], ""), Outcome = outcome_label(antigen),
    `Range, unadjusted` = range_unadjusted, `Range, raking` = range_raking, `Range, regression standardisation` = range_regression,
    `RDS minus CBE, unadjusted` = difference_unadjusted, `RDS minus CBE, raking` = difference_raking,
    `RDS minus CBE, regression standardisation` = difference_regression
  )

save_table(table_I03, "output/tables/table_I03.csv")
