# 09_network.R
# The network module: what it collected, and network scale-up estimators
# against the household survey (Phase 1) and the directly reported estimates
# Reproduces: Figure 16, Annex Tables K1-K8

# Setup ----

pacman::p_load(survey, RDS, patchwork, tidyverse)

source("code/functions.R")

# bootstrap replicates
replicates <- 2000

# Data ----

frames <- read_rds("data/clean/frames.rds")
network <- read_rds("data/clean/network.rds")
tracking <- read_rds("data/clean/rds_tracking.rds")
estimates <- read_rds("data/clean/estimates.rds")

# the Phase 1 samples with a network module, and the three primary antigens
samples <- c("eth_phase1_hh", "som_phase1_hh", "eth_phase1_rds", "som_phase1_rds", "eth_phase1_vbs", "som_phase1_vbs")
outcomes <- antigen_outcomes[c("Penta3", "OPV3", "MCV1")]

# Analysis ----

## Direct arm ----
# Each respondent's weighted sum of children with the outcome (num) and of
# children with an observed outcome (den), so that sum(num) / sum(den) is the
# sample's directly reported estimate: equal weights for HH, mobility weights
# for VBS, and for RDS the successive-sampling weight of the respondent times
# their number of children.
# respondent sums of the direct arm, with each RDS respondent's recruitment
# chain (the seed it descends from)
direct_arm <- function(sample, cohort, outcomes) {
  if (str_detect(sample, "hh|vbs")) {
    weight <- if (str_detect(sample, "vbs")) vbs_weight(cohort, sample) else rep(1, nrow(cohort))
    ids <- as.character(cohort$record_ix)
    num <- map(outcomes, \(field) if_else(is.na(cohort[[field]]), 0, cohort[[field]] * weight)) |> bind_cols() |> as.matrix()
    den <- map(outcomes, \(field) if_else(is.na(cohort[[field]]), 0, weight)) |> bind_cols() |> as.matrix()
    num <- rowsum(num, ids, reorder = FALSE)
    den <- rowsum(den, ids, reorder = FALSE)
    cluster <- if (str_detect(sample, "hh")) as.character(cohort$cluster_id[match(rownames(num), ids)]) else NULL
    return(list(ids = rownames(num), num = num, den = den, cluster = cluster, chain = NULL))
  }
  nodes <- rds_network(cohort, antigen_outcomes, tracking = if (sample == "som_phase2_rds") tracking) |> as.data.frame()
  population <- population_size[[sample]]
  design <- as.rds.data.frame(nodes, id = "id_participant_calc", recruiter.id = "id_referring", network.size = "degree", population.size = population)
  ids <- as.character(nodes$id_participant_calc)
  seed <- sample_seed(sample)
  num <- matrix(0, length(ids), length(outcomes), dimnames = list(ids, names(outcomes)))
  den <- num
  for (antigen in names(outcomes)) {
    field <- outcomes[[antigen]]
    seed_field <- if (field == "zero_dose") "penta_1p_clean" else field
    fit <- withr::with_seed(derived_seed(seed, paste0("rds_outcome/", seed_field)), RDS.SS.estimates(design, outcome.variable = field, N = population))
    child_weight <- as.numeric(fit$weights) * node_children(cohort, ids, field)
    use <- !is.na(child_weight) & !is.na(nodes[[field]])
    num[use, antigen] <- nodes[[field]][use] * child_weight[use]
    den[use, antigen] <- child_weight[use]
  }
  list(ids = ids, num = num, den = den, cluster = NULL, chain = chain_root(nodes))
}

## Network arm ----
# For each respondent: the network children reported on (r), the stated number
# of network children (dg), and per antigen the children reported vaccinated
# (y) and those whose status was known (k). For the visibility adjustment, the
# respondent's own children with the outcome (kcd1) and those weighted by the
# share of the respondent's answered network-child rows where the child's
# caregiver was said to know the respondent's own child's status (kc1). Every
# enrolled respondent is marked as a network respondent (is_network).
# assemble one respondent matrix per sample
respondent_matrix <- function(sample, outcomes) {
  id_field <- if (str_detect(sample, "rds")) "id_participant_calc" else "record_ix"
  cohort <- direct_cohort(frames[[sample]], sample)
  enrolled <- network[[sample]] |> filter(consent %in% "Enrolled")
  respondents <- enrolled |> distinct(record_ix, .keep_all = TRUE)
  respondent_ids <- as.character(respondents[[id_field]])
  reports <- enrolled |> filter(!is.na(record_ix_nchild)) |> mutate(zero_dose = 1 - penta_1p_clean)
  report_ids <- as.character(reports[[id_field]])

  direct <- direct_arm(sample, cohort, outcomes)
  pool <- if (str_detect(sample, "vbs")) union(respondent_ids, direct$ids) else direct$ids

  columns <- c("r", "dg", "is_network", as.vector(outer(names(outcomes), c("d_num", "d_den", "y", "k", "kc1", "kcd1"), paste, sep = ".")))
  x <- matrix(0, length(pool), length(columns), dimnames = list(pool, columns))
  at <- match(direct$ids, pool)
  for (antigen in names(outcomes)) {
    x[at, paste0(antigen, ".d_num")] <- direct$num[, antigen]
    x[at, paste0(antigen, ".d_den")] <- direct$den[, antigen]
  }

  # network children reported and the stated degree of respondents with at least one
  rows_by <- rowsum(rep(1, nrow(reports)), report_ids, reorder = FALSE)
  x[match(rownames(rows_by), pool), "r"] <- rows_by[, 1]
  degree <- replace_na(as.numeric(respondents$nsum_child_total), 0)
  degree_at <- match(respondent_ids, pool)
  x[degree_at, "dg"] <- if_else(x[degree_at, "r"] > 0, degree, 0)
  x[degree_at, "is_network"] <- 1

  # vaccinated and known-status network children
  for (antigen in names(outcomes)) {
    values <- reports[[outcomes[[antigen]]]]
    y_by <- rowsum(as.numeric(values %in% 1), report_ids, reorder = FALSE)
    k_by <- rowsum(as.numeric(!is.na(values)), report_ids, reorder = FALSE)
    x[match(rownames(y_by), pool), paste0(antigen, ".y")] <- y_by[, 1]
    x[match(rownames(k_by), pool), paste0(antigen, ".k")] <- k_by[, 1]
  }

  # share of answered network-child rows where the caregiver knows, applied to own vaccinated children
  knows <- as.character(reports$nsum_caregiver_knows)
  yes_by <- rowsum(as.numeric(knows %in% "Yes"), report_ids, reorder = FALSE)
  answered_by <- rowsum(as.numeric(!is.na(knows)), report_ids, reorder = FALSE)
  share <- set_names(if_else(answered_by[, 1] > 0, yes_by[, 1] / answered_by[, 1], NA_real_), rownames(answered_by))
  own_ids <- as.character(cohort[[id_field]])
  for (antigen in names(outcomes)) {
    vaccinated_by <- rowsum(as.numeric(cohort[[outcomes[[antigen]]]] %in% 1), own_ids, reorder = FALSE)
    ids <- rownames(vaccinated_by)
    respondent_share <- unname(share[ids])
    x[match(ids, pool), paste0(antigen, ".kc1")] <- if_else(is.na(respondent_share), 0, vaccinated_by[, 1] * respondent_share)
    x[match(ids, pool), paste0(antigen, ".kcd1")] <- if_else(is.na(respondent_share), 0, vaccinated_by[, 1])
  }

  unit <- if (str_detect(sample, "hh")) set_names(direct$cluster, direct$ids)[pool] else pool
  chain <- if (str_detect(sample, "rds")) unname(set_names(direct$chain, direct$ids)[pool])
  list(x = x, unit = unname(unit), chain = chain)
}

## Estimators ----
# Known-status prevalence (y / k), the NSUM estimator (y / dg) and the
# visibility-adjusted NSUM estimator (y / dg divided by kc1 / kcd1, undefined
# above 1), each against the direct estimate sum(num) / sum(den).
# the three estimators and the direct estimate from column sums
estimators <- function(sums, antigen) {
  col <- \(name) sums[, paste0(antigen, ".", name)]
  ratio <- \(num, den) if_else(den > 0, num / den, NA_real_)
  nsum <- ratio(col("y"), sums[, "dg"])
  visibility <- ratio(col("kc1"), col("kcd1"))
  adjusted <- nsum / visibility
  adjusted[is.na(visibility) | visibility <= 0 | is.na(adjusted) | adjusted < 0 | adjusted > 1] <- NA_real_
  list(
    known_status = ratio(col("y"), col("k")), nsum = nsum, visibility_adjusted = adjusted,
    direct = ratio(col("d_num"), col("d_den"))
  )
}

## Bootstrap ----
# Units (HH clusters, VBS and RDS respondents) resampled with replacement,
# 2,000 times, re-estimating every quantity in each replicate; the interval for
# the difference from the direct estimate is the 2.5th to 97.5th percentile of
# the paired replicate differences (missing when fewer than 95% are finite).
# percentile interval
percentile <- function(x) {
  finite <- x[is.finite(x)]
  if (length(x) == 0 || length(finite) < 0.95 * length(x)) return(c(NA_real_, NA_real_))
  unname(quantile(finite, c(0.025, 0.975), type = 7))
}

# estimate each sample, antigen and estimator with its bootstrap interval
network_estimates <- map(samples, \(sample) {
  built <- respondent_matrix(sample, outcomes)
  unit_sums <- rowsum(built$x, built$unit, reorder = FALSE)
  n_units <- nrow(unit_sums)
  scheme <- if (str_detect(sample, "hh")) "cluster" else "respondent"
  draws <- withr::with_seed(
    derived_seed(study_seed, paste0("proxy_network_estimators/", str_replace(sample, "_phase(\\d)_", "_phase_\\1_"), "/", scheme)),
    rmultinom(replicates, n_units, rep(1 / n_units, n_units))
  )
  replicate_sums <- crossprod(draws, unit_sums)
  full_sums <- matrix(colSums(built$x), nrow = 1, dimnames = list(NULL, colnames(built$x)))
  map(names(outcomes), \(antigen) {
    point <- estimators(full_sums, antigen)
    replicate <- estimators(replicate_sums, antigen)
    vaccinated_own <- full_sums[1, paste0(antigen, ".kcd1")]
    map(c("known_status", "nsum", "visibility_adjusted"), \(estimator) {
      interval <- percentile(replicate[[estimator]] - replicate$direct)
      tibble(
        estimator = estimator, difference_pp = 100 * (point[[estimator]] - point$direct),
        lower_pp = 100 * interval[[1]], upper_pp = 100 * interval[[2]],
        small_arm = estimator == "visibility_adjusted" & vaccinated_own < 30
      )
    }) |>
      list_rbind() |>
      mutate(antigen = antigen)
  }) |>
    list_rbind() |>
    mutate(sample = sample)
}) |>
  list_rbind()

## Network-module descriptions ----
# What the network module collected in each sample: the respondents (every
# enrolled respondent) and their network-child reports, the share of reports
# with a known status by antigen, the certainty ratings, the women reported by
# relationship category, network size, and reports that match another report
# on the child's name, place and sex.
network_samples <- intersect(study_samples, names(network))
phase1_samples <- intersect(study_samples, samples)

# enrolled respondents and their network-child reports in each sample
network_respondents <- map(set_names(network_samples), \(sample) network[[sample]] |> filter(consent %in% "Enrolled") |> distinct(record_ix, .keep_all = TRUE))
network_reports <- map(set_names(network_samples), \(sample) network[[sample]] |> filter(consent %in% "Enrolled", !is.na(record_ix_nchild)))

# share of reports with a known status, by antigen
known_status <- map(network_samples, \(sample) {
  reports <- network_reports[[sample]]
  tibble(
    sample = sample, n_respondents = nrow(network_respondents[[sample]]), n_reports = nrow(reports),
    antigen = c("Penta3", "OPV3", "MCV1", "Penta1"),
    known = map_dbl(antigen_outcomes[antigen], \(field) mean(!is.na(reports[[field]]))),
    no_penta_never = !any(as.character(reports$nsum_penta_ever) %in% "No")
  )
}) |>
  list_rbind()

# certainty of the reported status of each vaccine; not asked where the
# respondent did not know whether the child had received it
certainty_levels <- c("Very certain", "Somewhat certain", "Not certain")
certainty <- map(network_samples, \(sample) {
  reports <- network_reports[[sample]]
  map(c(Pentavalent = "penta", OPV = "opv", Measles = "measles"), \(family) {
    rating <- as.character(reports[[paste0("nsum_", family, "_certainty")]])
    ever <- as.character(reports[[paste0("nsum_", family, "_ever")]])
    shares <- map_dbl(certainty_levels, \(level) mean(rating %in% level))
    tibble(level = c(certainty_levels, "Not asked"), share = c(shares, mean(is.na(rating) & (is.na(ever) | ever == "Do not know"))))
  }) |>
    list_rbind(names_to = "vaccine") |>
    mutate(sample = sample)
}) |>
  list_rbind()

# women reported by relationship category, Phase 1
relationship_labels <- c(
  closefam = "Close family", extfam = "Extended family", inlaw = "In-laws", friends = "Friends", cowork = "Co-workers",
  relig = "Religious community", neighbors = "Neighbours", services = "Service providers", other = "Other"
)
relationships <- map(phase1_samples, \(sample) {
  respondents <- network_respondents[[sample]]
  women <- map_dbl(names(relationship_labels), \(category) sum(as.numeric(respondents[[paste0("nsum_", category)]]), na.rm = TRUE))
  tibble(sample = sample, category = names(relationship_labels), share = women / sum(as.numeric(respondents$nsum_total), na.rm = TRUE), total = sum(as.numeric(respondents$nsum_total), na.rm = TRUE))
}) |>
  list_rbind()

# network size and matching reports, Phase 1
network_size <- map(phase1_samples, \(sample) {
  size <- as.numeric(network_respondents[[sample]]$nsum_total)
  key <- network_reports[[sample]]$nsum_child_key
  tibble(
    sample = sample, n_respondents = length(size), no_member = sum(size == 0, na.rm = TRUE) / length(size),
    median = median(size, na.rm = TRUE), q1 = quantile(size, 0.25, na.rm = TRUE, names = FALSE), q3 = quantile(size, 0.75, na.rm = TRUE, names = FALSE),
    maximum = max(size, na.rm = TRUE), matching = mean(duplicated(key) | duplicated(key, fromLast = TRUE))
  )
}) |>
  list_rbind()

# known status by whether the network child's caregiver knew the status of
# the respondent's own child, Phase 1
caregiver_knows <- map(phase1_samples, \(sample) {
  reports <- network_reports[[sample]] |> mutate(answer = replace_na(as.character(nsum_caregiver_knows), "Missing"))
  reports |>
    summarise(
      n_reports = n(),
      Penta3 = mean(!is.na(penta_3p_clean)), OPV3 = mean(!is.na(opv_3p_clean)), MCV1 = mean(!is.na(measles_1p_clean)),
      `Zero-dose` = mean(!is.na(penta_1p_clean)),
      .by = answer
    ) |>
    mutate(sample = sample)
}) |>
  list_rbind() |>
  mutate(answer = factor(answer, levels = c("Yes", "No", "Do not know", "Missing"))) |>
  arrange(factor(sample, levels = phase1_samples), answer)

## Network estimators for the annex ----
# Each estimator for every network sample of Phase 1 and the Phase 2 RDS
# samples, on the four outcomes (zero-dose from the reports of no pentavalent
# dose), with its own interval: a linearised ratio standard error with a Wilson
# interval for known-status prevalence and the NSUM estimator (t on the design
# degrees of freedom for HH and VBS, 1.96 for RDS), and the bootstrap
# percentile interval for the visibility-adjusted estimator. RDS estimates also
# carry intervals from resampling whole recruitment chains.
annex_outcomes <- antigen_outcomes[primary_outcomes]
annex_samples <- c(samples, "eth_phase2_rds", "som_phase2_rds")
estimator_names <- c("known_status", "nsum", "visibility_adjusted")

# linearised standard error of a ratio of totals over sampling units
ratio_se <- function(num, den) {
  m <- length(num)
  if (m < 2 || sum(den) <= 0) return(NA_real_)
  ratio <- sum(num) / sum(den)
  sqrt(m / (m - 1) * sum((num - ratio * den)^2)) / sum(den)
}

# Wilson interval at a critical value, count-based at an estimate of 0 or 1
# with a positive standard error
network_wilson <- function(estimate, se, crit, n_events, n_domain) {
  if (is.na(estimate) || is.na(se) || estimate < 0 || estimate > 1) return(c(NA_real_, NA_real_))
  if (se > 0 && estimate %in% c(0, 1)) return(count_wilson(round(n_events), round(n_domain), crit))
  wilson_interval(estimate, se, crit)
}

# replicate column sums from resampling the rows of a unit matrix
resample_sums <- function(unit_sums, seed) {
  n_units <- nrow(unit_sums)
  draws <- withr::with_seed(seed, rmultinom(replicates, n_units, rep(1 / n_units, n_units)))
  crossprod(draws, unit_sums)
}

# estimate every outcome and estimator of one sample
network_cells <- map(annex_samples, \(sample) {
  built <- respondent_matrix(sample, annex_outcomes)
  x <- built$x
  key <- paste0("proxy_network_estimators/", str_replace(sample, "_phase(\\d)_", "_phase_\\1_"), "/")
  scheme <- if (str_detect(sample, "hh")) "cluster" else "respondent"
  replicate_sums <- resample_sums(rowsum(x, built$unit, reorder = FALSE), derived_seed(study_seed, paste0(key, scheme)))
  chain_sums <- if (!is.null(built$chain)) resample_sums(rowsum(x, built$chain, reorder = FALSE), derived_seed(study_seed, paste0(key, "seed_chain")))
  full_sums <- matrix(colSums(x), nrow = 1, dimnames = list(NULL, colnames(x)))
  network_rows <- x[, "is_network"] > 0
  network_units <- rowsum(x[network_rows, , drop = FALSE], built$unit[network_rows], reorder = FALSE)
  crit <- if (str_detect(sample, "rds")) 1.96 else qt(0.975, length(unique(built$unit[network_rows])) - 1)
  total <- \(name) full_sums[1, name]

  map(names(annex_outcomes), \(antigen) {
    point <- estimators(full_sums, antigen)
    replicate <- estimators(replicate_sums, antigen)
    chain <- if (!is.null(chain_sums)) estimators(chain_sums, antigen)
    y <- total(paste0(antigen, ".y"))
    k <- total(paste0(antigen, ".k"))
    r <- total("r")
    unbounded_adjusted <- point$nsum / (total(paste0(antigen, ".kc1")) / total(paste0(antigen, ".kcd1")))
    map(estimator_names, \(estimator) {
      estimate <- point[[estimator]]
      n_domain <- if (estimator == "known_status") k else total("dg")
      if (estimator == "visibility_adjusted") {
        se <- sd(replicate[[estimator]][is.finite(replicate[[estimator]])])
        interval <- percentile(replicate[[estimator]])
      } else {
        se <- ratio_se(network_units[, paste0(antigen, ".y")], network_units[, if (estimator == "known_status") paste0(antigen, ".k") else "dg"])
        interval <- network_wilson(estimate, se, crit, y, n_domain)
      }
      chain_interval <- if (is.null(chain)) c(NA_real_, NA_real_) else percentile(chain[[estimator]])
      paired <- percentile(replicate[[estimator]] - replicate$direct)
      chain_paired <- if (is.null(chain)) c(NA_real_, NA_real_) else percentile(chain[[estimator]] - chain$direct)
      tibble(
        sample = sample, antigen = antigen, estimator = estimator,
        estimate = estimate, se = se, lower = interval[[1]], upper = interval[[2]], n_domain = n_domain, n_events = y,
        chain_se = if (is.null(chain)) NA_real_ else sd(chain[[estimator]][is.finite(chain[[estimator]])]),
        chain_lower = chain_interval[[1]], chain_upper = chain_interval[[2]],
        direct = point$direct,
        paired_lower = paired[[1]], paired_upper = paired[[2]], chain_paired_lower = chain_paired[[1]], chain_paired_upper = chain_paired[[2]],
        bound_lower = y / r, bound_upper = (r - (k - y)) / r,
        above_one = estimator == "visibility_adjusted" & is.finite(unbounded_adjusted) & unbounded_adjusted > 1,
        small_arm = estimator == "visibility_adjusted" & total(paste0(antigen, ".kcd1")) < 30
      )
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  left_join(sample_parts(annex_samples), by = "sample") |>
  mutate(interpretable = !is.na(estimate) & !small_arm & !(country == "eth" & antigen == "Zero-dose"))

## Comparison with the household survey ----
# Phase 1: within the household survey, the network estimate against the
# directly reported estimate of the same sample, by the paired bootstrap; for
# RDS and VBS, against the household survey by Newcombe's interval (count-based
# for a boundary arm), RDS called distinguishable only when the interval from
# resampling recruitment chains also excludes zero.
# household-survey reference with its counts
counts <- event_counts(frames)
reference <-
  estimates |>
  filter(method == "HH") |>
  transmute(country, antigen = as.character(antigen), ref_estimate = estimate, ref_se = standard_error, ref_lower = ci_lower, ref_upper = ci_upper, sample) |>
  left_join(counts |> rename(ref_n_domain = n_domain, ref_n_events = n_events), by = c("sample", "antigen")) |>
  select(-sample)

# verdict of each Phase 1 cell
contains_zero <- \(lower, upper) lower <= 0 & upper >= 0
phase1_cells <-
  network_cells |>
  filter(phase == 1) |>
  left_join(reference, by = c("country", "antigen")) |>
  mutate(error_pp = 100 * (estimate - ref_estimate))
phase1_cells$respondent_nd <- NA
phase1_cells$chain_nd <- NA
for (i in seq_len(nrow(phase1_cells))) {
  cell <- phase1_cells[i, ]
  if (is.na(cell$estimate)) next
  if (cell$method == "hh") {
    phase1_cells$respondent_nd[i] <- contains_zero(cell$paired_lower, cell$paired_upper)
    next
  }
  if (is.na(cell$lower)) next
  ref <- boundary_arm(cell$ref_estimate, cell$ref_se, cell$ref_lower, cell$ref_upper, cell$ref_n_domain, cell$ref_n_events)
  arm <- boundary_arm(cell$estimate, cell$se, cell$lower, cell$upper, cell$n_domain, cell$n_events)
  interval <- newcombe_interval(cell$estimate, arm$ci_lower, arm$ci_upper, cell$ref_estimate, ref$ci_lower, ref$ci_upper)
  phase1_cells$respondent_nd[i] <- contains_zero(interval$lower, interval$upper)
  if (cell$method == "rds" && !is.na(cell$chain_lower)) {
    chain_arm <- boundary_arm(cell$estimate, cell$chain_se, cell$chain_lower, cell$chain_upper, cell$n_domain, cell$n_events)
    chain_interval <- newcombe_interval(cell$estimate, chain_arm$ci_lower, chain_arm$ci_upper, cell$ref_estimate, ref$ci_lower, ref$ci_upper)
    phase1_cells$chain_nd[i] <- contains_zero(chain_interval$lower, chain_interval$upper)
  }
}
phase1_cells$not_different <- if_else(phase1_cells$method == "rds", phase1_cells$respondent_nd | phase1_cells$chain_nd %in% TRUE, phase1_cells$respondent_nd)

# the directly reported estimate of each method against the household survey
direct_performance <- phase1_performance(estimates, counts)

## Agreement with direct report, Phase 2 RDS ----
# Each network estimate against the sample's own directly reported estimate:
# agreement when the paired respondent-resampled or the chain-resampled
# interval of the difference contains zero.
# score each Phase 2 cell
phase2_cells <-
  network_cells |>
  filter(phase == 2, interpretable) |>
  mutate(
    respondent_agrees = contains_zero(paired_lower, paired_upper),
    agrees = respondent_agrees | contains_zero(chain_paired_lower, chain_paired_upper) %in% TRUE,
    bounds_contain = direct >= bound_lower & direct <= bound_upper,
    difference_pp = 100 * (estimate - direct)
  )

# Tables ----

## Annex Table K1 ----
# Respondents, network-child reports and the share of reports with a known
# status, by antigen; a double bar marks a sample with no report that a child
# never received a pentavalent dose.
# label of a sample: country, phase
sample_label <- \(sample) paste0(country_labels[str_sub(sample, 1, 3)], ", ", str_extract(sample, "(?<=phase)\\d"))
method_label <- \(sample) str_to_upper(str_extract(sample, "[a-z]+$"))

# one row per sample, antigens across columns
table_K01 <-
  known_status |>
  mutate(value = paste0(format_number(100 * known), if_else(antigen == "Penta1" & no_penta_never, "‖", ""))) |>
  select(sample, n_respondents, n_reports, antigen, value) |>
  pivot_wider(names_from = antigen, values_from = value) |>
  transmute(
    `Country, phase` = sample_label(sample), Method = method_label(sample),
    Respondents = format_count(n_respondents), `Network-child reports` = format_count(n_reports), Penta3, OPV3, MCV1, Penta1
  )

save_table(table_K01, "output/tables/table_K01.csv")

## Annex Table K2 ----
# Certainty rating of the reported status of network children, by vaccine (%
# of network-child reports).
# one row per sample and vaccine
table_K02 <-
  certainty |>
  mutate(value = format_number(100 * share)) |>
  select(sample, vaccine, level, value) |>
  pivot_wider(names_from = level, values_from = value) |>
  mutate(first = vaccine == "Pentavalent") |>
  transmute(
    `Country, phase` = if_else(first, sample_label(sample), ""), Method = if_else(first, method_label(sample), ""), Vaccine = vaccine,
    `Very certain`, `Somewhat certain`, `Not certain`, `Not asked`
  )

save_table(table_K02, "output/tables/table_K02.csv")

## Annex Table K3 ----
# Women reported in the network module by relationship category (% of all
# women reported), with the number reported, Phase 1.
# one row per category, samples across columns
column_label <- \(sample) paste(country_labels[str_sub(sample, 1, 3)], method_label(sample))
table_K03 <-
  bind_rows(
    relationships |> transmute(sample, row = relationship_labels[category], value = format_number(100 * share)),
    relationships |> distinct(sample, total) |> transmute(sample, row = "All categories (members reported)", value = format_count(total))
  ) |>
  mutate(sample = column_label(sample)) |>
  pivot_wider(names_from = sample, values_from = value) |>
  rename(`Relationship category` = row)

save_table(table_K03, "output/tables/table_K03.csv")

## Annex Table K4 ----
# Network size, respondents reporting no network member, and network-child
# reports matching another report, Phase 1.
# one row per sample
table_K04 <-
  network_size |>
  transmute(
    Country = country_labels[str_sub(sample, 1, 3)], Method = method_label(sample), Respondents = format_count(n_respondents),
    `Reporting no network member (%)` = format_number(100 * no_member),
    `Network size, median (IQR)` = sprintf("%s (%s–%s)", median, q1, q3),
    `Network size, maximum` = as.character(maximum),
    `Network-child reports matching another report (%)` = format_number(100 * matching)
  )

save_table(table_K04, "output/tables/table_K04.csv")

## Annex Table K5 ----
# Network-child reports with a known status by the answer on whether the
# child's caregiver knew the status of the respondent's own child, Phase 1 (%).
# one row per sample and answer given
table_K05 <-
  caregiver_knows |>
  mutate(first = !duplicated(sample)) |>
  transmute(
    Country = if_else(first, country_labels[str_sub(sample, 1, 3)], ""), Method = if_else(first, method_label(sample), ""),
    `Caregiver knew the status of the respondent's child` = as.character(answer),
    `Network-child reports` = format_count(n_reports),
    Penta3 = format_number(100 * Penta3), OPV3 = format_number(100 * OPV3), MCV1 = format_number(100 * MCV1),
    `Zero-dose status` = format_number(100 * `Zero-dose`)
  )

save_table(table_K05, "output/tables/table_K05.csv")

## Annex Table K6 ----
# RMSE over the three primary antigens against the household-survey estimate
# (pp), Phase 1, by estimator, method and country; within the household survey
# the reference is the directly reported estimate of the same sample. (a) marks
# a cell left out because the visibility-adjusted estimate exceeded 100%, (b) a
# cell entered as a point without an interval.
# row labels of the estimators
estimator_labels_annex <- c(
  direct = "Directly reported (method's own)", known_status = "Known-status prevalence", nsum = "NSUM estimator",
  visibility_adjusted = "Visibility-adjusted NSUM estimator"
)
primary_cells <- phase1_cells |> filter(antigen %in% names(outcomes))

# RMSE of each estimator over its interpretable cells
network_rmse <-
  primary_cells |>
  summarise(
    rmse = sqrt(mean(error_pp[interpretable]^2)),
    marker = case_when(any(above_one) ~ " (a)", any(interpretable & is.na(lower)) ~ " (b)", .default = ""),
    .by = c(sample, country, method, estimator)
  ) |>
  mutate(value = paste0(format_number(rmse, 2), marker))
direct_rmse <-
  phase1_summary(direct_performance) |>
  filter(scope == "primary_antigens", country != "pooled", method != "CBE") |>
  transmute(sample = paste0(country, "_phase1_", str_to_lower(method)), estimator = "direct", value = format_number(rmse_pp, 2))

# one row per estimator, samples across columns
column_order <- paste(rep(c("HH", "RDS", "VBS"), each = 2), rep(c("Ethiopia", "Somalia"), times = 3))
sample_column <- \(sample) paste(method_label(sample), country_labels[str_sub(sample, 1, 3)])
estimator_table <- function(rows) {
  rows |>
    mutate(column = factor(sample_column(sample), levels = column_order), Estimator = factor(estimator_labels_annex[estimator], levels = estimator_labels_annex)) |>
    select(Estimator, column, value) |>
    arrange(column) |>
    pivot_wider(names_from = column, values_from = value, values_fill = "reference") |>
    arrange(Estimator) |>
    mutate(Estimator = as.character(Estimator))
}
table_K06 <- estimator_table(bind_rows(direct_rmse, network_rmse))

save_table(table_K06, "output/tables/table_K06.csv")

## Annex Table K7 ----
# Phase 1 primary-antigen comparisons not distinguishable from the household
# survey, by estimator, method and country, and those whose bounds contain the
# household-survey estimate; for RDS, in brackets, those not distinguishable on
# the respondent-resampled interval alone.
# count the comparisons with a verdict
not_distinguishable <-
  primary_cells |>
  filter(interpretable, !is.na(not_different)) |>
  summarise(
    value = paste0(
      sum(not_different), " of ", n(),
      if_else(first(method) == "rds", paste0(" (", sum(respondent_nd), ")"), "")
    ),
    .by = c(sample, estimator)
  )
direct_verdicts <-
  direct_performance |>
  filter(antigen != "Zero-dose", method != "CBE") |>
  summarise(value = paste(sum(not_different), "of", n()), .by = c(country, method)) |>
  transmute(sample = paste0(country, "_phase1_", str_to_lower(method)), estimator = "direct", value)
bounds_contain <-
  primary_cells |>
  filter(estimator == "known_status") |>
  summarise(value = paste(sum(ref_estimate >= bound_lower & ref_estimate <= bound_upper), "of", n()), .by = sample) |>
  mutate(estimator = "bounds")

# one row per estimator and the bounds, samples across columns
estimator_labels_annex <- c(estimator_labels_annex, bounds = "Bounds contain the household-survey estimate")
table_K07 <- estimator_table(bind_rows(direct_verdicts, not_distinguishable, bounds_contain))

save_table(table_K07, "output/tables/table_K07.csv")

## Annex Table K8 ----
# Phase 2 RDS comparisons in which the network estimate agreed with the
# directly reported estimate (Penta3, OPV3 and MCV1 in each country and
# zero-dose prevalence in Somalia), with the mean absolute difference over the
# primary-antigen comparisons.
# one row per estimator
table_K08 <-
  phase2_cells |>
  summarise(
    `Agreed (either interval)` = paste(sum(agrees), "of", n()),
    `Agreed (respondent-resampled interval alone)` = paste(sum(respondent_agrees), "of", n()),
    `Bounds contained the direct estimate` = paste(sum(bounds_contain), "of", n()),
    `Mean absolute difference, primary antigens (pp)` = format_number(mean(abs(difference_pp[antigen != "Zero-dose"])), 2),
    .by = estimator
  ) |>
  transmute(Estimator = estimator_labels_annex[estimator], `Agreed (either interval)`, `Agreed (respondent-resampled interval alone)`,
            `Bounds contained the direct estimate`, `Mean absolute difference, primary antigens (pp)`)

save_table(table_K08, "output/tables/table_K08.csv")

# Figures ----

## Figure 16 ----
# Difference of each network estimator from the directly reported estimate,
# with its bootstrap interval, by country, method and antigen. Estimates shown
# are those that are defined and, for the visibility adjustment, rest on at
# least 30 of the respondents' own children.
# keep the interpretable estimates
estimator_labels <- c(known_status = "Known-status prevalence", nsum = "NSUM estimator", visibility_adjusted = "Visibility-adjusted NSUM estimator")
methods <- c("hh", "rds", "vbs")
plotted <-
  network_estimates |>
  filter(!is.na(difference_pp), !small_arm) |>
  mutate(
    country = factor(country_labels[str_sub(sample, 1, 3)], levels = c("Ethiopia", "Somalia")),
    method = factor(str_extract(sample, "[a-z]+$"), levels = methods),
    estimator = factor(estimator_labels[estimator], levels = estimator_labels),
    antigen = factor(antigen, levels = names(outcomes))
  )

# each method a block of three antigen rows under a header, the estimators dodged within a row
stride <- 4.6
y_limits <- c(-(stride * 2 + 3.6), 0.6)
dodge <- 0.03 * diff(y_limits)
plotted <- mutate(
  plotted,
  y = -(stride * (as.integer(method) - 1) + as.integer(antigen)) - (as.integer(estimator) - 2) * dodge
)
headers <- expand_grid(country = factor(c("Ethiopia", "Somalia"), levels = c("Ethiopia", "Somalia")), method = factor(methods, levels = methods)) |>
  mutate(y = -(stride * (as.integer(method) - 1)))

# one horizontal scale for both countries, rounded out to 10 pp
x_limits <- c(
  floor(min(plotted$lower_pp, plotted$difference_pp, na.rm = TRUE) / 10) * 10,
  ceiling(max(plotted$upper_pp, plotted$difference_pp, na.rm = TRUE) / 10) * 10
)

# draw the differences
figure_16 <-
  ggplot(plotted, aes(y = y)) +
  geom_vline(xintercept = 0, colour = ink[["secondary"]], linewidth = 0.4) +
  geom_point(
    data = headers, aes(x = x_limits[[1]] + 0.04 * diff(x_limits), y = y, fill = method), inherit.aes = FALSE,
    shape = 21, colour = "white", size = 2.6, stroke = 0.6, show.legend = FALSE
  ) +
  geom_text(
    data = headers, aes(x = x_limits[[1]] + 0.075 * diff(x_limits), y = y, label = method_labels[as.character(method)]), inherit.aes = FALSE,
    hjust = 0, vjust = 0.5, fontface = "bold", size = 2.9, family = report_font, colour = ink[["primary"]]
  ) +
  geom_linerange(aes(xmin = lower_pp, xmax = upper_pp, colour = method), linewidth = 0.45, na.rm = TRUE) +
  geom_point(aes(x = difference_pp, shape = estimator, colour = method), size = 2) +
  facet_grid(country ~ ., switch = "y") +
  scale_colour_manual(values = method_colours, guide = "none") +
  scale_fill_manual(values = method_colours, guide = "none") +
  scale_shape_manual(
    values = set_names(c(17, 15, 18), estimator_labels), name = "Estimator", drop = FALSE,
    guide = guide_legend(nrow = 1, override.aes = list(colour = ink[["secondary"]]))
  ) +
  scale_y_continuous(breaks = -(rep(stride * 0:2, each = 3) + 1:3), labels = rep(names(outcomes), times = 3), limits = y_limits, expand = c(0, 0)) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = "Difference from direct (pp)", y = NULL) +
  theme_report() +
  theme(
    panel.grid.major.y = element_blank(), panel.spacing.y = unit(2, "lines"),
    panel.widths = unit(10, "cm"), axis.text.y = element_text(size = 7.5)
  )

# collect the legend below the panel and draw the country divider
figure_16 <- wrap_plots(figure_16, guides = "collect") & theme(legend.position = "bottom")
save_figure(country_divider(figure_16, height_cm = 13), "output/figures/figure_16.png", height_cm = 13)
