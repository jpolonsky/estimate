# 04_coverage.R
# Phase 1 coverage and zero-dose prevalence by country and method
# Reproduces: Table 8, Figure 11, Annex Tables G1-G10, Annex Figures G1-G5

# Setup ----

pacman::p_load(survey, RDS, patchwork, tidyverse)

source("code/functions.R")

# Data ----

estimates <- read_rds("data/clean/estimates.rds")
frames <- read_rds("data/clean/frames.rds")

# the fourteen samples' directly reported children
cohorts <- imap(frames[study_samples], direct_cohort)

# the six samples drawn in clusters, Ethiopia first
cluster_samples <- c("eth_phase1_hh", "eth_phase1_cbe", "eth_phase2_cbe", "som_phase1_hh", "som_phase1_cbe", "som_phase2_cbe")

# Analysis ----

## Phase 1 estimates ----
# The three primary antigens and zero-dose prevalence for the four methods,
# in report order.
# keep the Phase 1 primary outcomes
phase1 <-
  estimates |>
  filter(phase == 1, antigen %in% primary_outcomes) |>
  mutate(antigen = fct_drop(fct_relevel(antigen, primary_outcomes)))

## Antigen missingness ----
# Children with a recorded status for each antigen, and the two reasons a
# status is missing: the caregiver did not know whether the child had received
# the vaccine, or did not know the number of doses.
# count recorded values and the reasons for missing ones
missingness <- imap(cohorts, \(cohort, sample) {
  recorded <- map_int(antigen_outcomes[c("Penta1", "Penta3", "OPV1", "OPV3", "MCV1", "MCV2")], \(field) sum(!is.na(cohort[[field]])))
  not_known <- map_int(c(penta = "penta_n", opv = "opv_n", measles = "measles_n"), \(field) sum(cohort[[field]] %in% "Vacc Unk"))
  doses_unknown <- map_int(c(penta = "penta_n", opv = "opv_n", measles = "measles_n"), \(field) sum(cohort[[field]] %in% "Doses Unk"))
  tibble(sample, children = nrow(cohort), antigen = names(recorded), recorded, not_known = list(not_known), doses_unknown = list(doses_unknown))
}) |>
  list_rbind()

## Timeliness ----
# Age at each dose against the schedule (6, 10 and 14 weeks for Penta and OPV,
# 9 and 15 months for MCV): on time within 2 months either side, late beyond.
# The k-th dose of a vaccine is the child's k-th dated dose (for OPV, the
# birth dose included).
# The forms give the age in days at each dose, to the day in CBE and Phase 2
# and to the month in the Phase 1 household survey, RDS and VBS (birth and dose
# months read at their middle in Ethiopia and at the 15th in Somalia).
# the doses, their schedule age (months) and the receipt indicator
doses <- tibble(
  dose = c("penta_1", "penta_2", "penta_3", "opv_1", "opv_2", "opv_3", "mcv_1", "mcv_2"),
  family = c("penta", "penta", "penta", "opv", "opv", "opv", "measles", "measles"),
  slot = c(1, 2, 3, 1, 2, 3, 1, 2),
  target = c(6, 10, 14, 6, 10, 14, NA, NA) * 7 / 30.4375,
  received = c("penta_1p_clean", "penta_2p", "penta_3p_clean", "opv_1p_clean", "opv_2p", "opv_3p_clean", "measles_1p_clean", "measles_2p_clean")
) |>
  mutate(target = if_else(family == "measles", c(9, 15)[slot], target))

# age at each dose (days) of each child who received it, with the date's source:
# the k-th dated dose, its age missing when the date of birth is not known
dose_ages <- imap(cohorts, \(cohort, sample) {
  parts <- sample_parts(sample)
  by_day <- parts$method == "cbe" | parts$phase == 2
  map(seq_len(nrow(doses)), \(k) {
    family <- doses$family[[k]]
    slots <- str_subset(names(cohort), str_c("^", family, "_d\\d_age_days$"))
    ages <- matrix(sapply(slots, \(s) cohort[[s]]), nrow = nrow(cohort))
    n_dated <- rowSums(matrix(sapply(str_replace(slots, "_age_days$", "_dated"), \(s) cohort[[s]]), nrow = nrow(cohort)))
    sorted <- t(apply(ages, 1, sort, na.last = TRUE))
    tibble(
      sample, dose = doses$dose[[k]], target = doses$target[[k]], by_day = by_day,
      received = cohort[[doses$received[[k]]]] %in% 1,
      age_days = sorted[, doses$slot[[k]]],
      dated = n_dated >= doses$slot[[k]],
      source = as.character(cohort[[str_c(family, "_source")]] %||% NA)
    )
  }) |>
    list_rbind()
}) |>
  list_rbind()

# share of dated doses on time (early folded in) and late, and of received doses undated
timeliness <-
  dose_ages |>
  filter(received) |>
  mutate(
    age_months = age_days / 30.4375,
    implausible = dated & (is.na(age_days) | if_else(by_day, age_days < 0, age_months < -1)),
    age_months = if_else(!by_day & age_months < 0, 0, age_months),
    classified = dated & !implausible,
    late = classified & age_months - target > 2
  ) |>
  summarise(
    n_received = n(), n_dated = sum(dated), n_classified = sum(classified), n_late = sum(late),
    .by = c(sample, dose)
  ) |>
  left_join(sample_parts(study_samples), by = "sample") |>
  mutate(
    undated = (n_received - n_dated) / n_received,
    late = if_else(n_classified > 0, n_late / n_classified, 0) * (1 - undated),
    timely = if_else(n_classified > 0, 1 - n_late / n_classified, 0) * (1 - undated)
  )

# age at each dose by source of the date, Ethiopia CBE, in 15-day bins
source_offsets <-
  dose_ages |>
  filter(sample %in% c("eth_phase1_cbe", "eth_phase2_cbe"), received, dated, !is.na(age_days), source %in% c("Vaccination card", "Respondent recall")) |>
  mutate(arm = if_else(source == "Vaccination card", "card", "recall"), bin = floor(age_days / 15) * 15) |>
  count(sample, dose, arm, bin) |>
  mutate(pct = n / sum(n), .by = c(sample, dose, arm))

## Per-cluster coverage ----
# Coverage in each visited cluster of the household survey and CBE, the
# estimate with each cluster left out, and whether the cluster's count lies
# outside the 95% and 99.8% binomial limits around the sample estimate.
# HH and CBE use equal weights, so each is a share of children.
cluster_coverage <- map(cluster_samples, \(sample) {
  cohort <- cohorts[[sample]]
  cohort$cluster <- if (str_detect(sample, "hh")) cohort$cluster_id else cohort$cluster
  map(primary_outcomes, \(antigen) {
    reference <- filter(estimates, .data$sample == .env$sample, .data$antigen == .env$antigen)
    y <- cohort[[antigen_outcomes[[antigen]]]]
    cohort |>
      mutate(y = .env$y) |>
      summarise(
        n_children = n(), n_observed = sum(!is.na(y)), n_events = sum(y == 1, na.rm = TRUE),
        .by = cluster
      ) |>
      arrange(cluster) |>
      mutate(
        sample = sample, antigen = antigen,
        estimate = if_else(n_observed > 0, n_events / n_observed, NA_real_),
        small_n = n_observed < 10,
        branch_estimate = reference$estimate, branch_ci_lower = reference$ci_lower, branch_ci_upper = reference$ci_upper,
        loo_delta_pp = 100 * ((sum(n_events) - n_events) / (sum(n_observed) - n_observed) - branch_estimate),
        outside_95 = n_events < qbinom(0.025, n_observed, branch_estimate) | n_events > qbinom(0.975, n_observed, branch_estimate),
        outside_998 = n_events < qbinom(0.001, n_observed, branch_estimate) | n_events > qbinom(0.999, n_observed, branch_estimate)
      )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  left_join(sample_parts(cluster_samples), by = "sample")

## Dropout between doses ----
# Retention on each step is the share of children with the first dose who
# also had the later one, estimated by each sample's own estimator among the
# children with the first dose; dropout is 100% minus retention. Penta1 to
# MCV1 is not nested, so the share of children with MCV1 but no Penta1 is
# estimated beside it. A child whose dose count is unknown on an antigen of a
# step (Penta3, OPV3 or MCV2) is left out of that step.
# the three steps and the doses they join
dropout_steps <- tibble(
  step = c("penta1_to_penta3", "opv1_to_opv3", "penta1_to_mcv1"),
  first = c("penta_1p_clean", "opv_1p_clean", "penta_1p_clean"),
  later = c("penta_3p_clean", "opv_3p_clean", "measles_1p_clean")
)
step_labels <- c(penta1_to_penta3 = "Penta1 to Penta3", opv1_to_opv3 = "OPV1 to OPV3", penta1_to_mcv1 = "Penta1 to MCV1", mcv1_without_penta1 = "MCV1 without Penta1")

# children reported as vaccinated whose dose count is not known
count_unknown <- \(cohort, family, field) as.character(cohort[[paste0(family, "_ever")]]) %in% "Yes" & is.na(cohort[[field]])

# add the step columns to a cohort
add_dropout_columns <- function(cohort) {
  unknown <- list(
    penta1_to_penta3 = count_unknown(cohort, "penta", "penta_3p_clean"),
    opv1_to_opv3 = count_unknown(cohort, "opv", "opv_3p_clean"),
    penta1_to_mcv1 = count_unknown(cohort, "penta", "penta_3p_clean") | count_unknown(cohort, "measles", "measles_2p_clean")
  )
  for (i in seq_len(nrow(dropout_steps))) {
    first <- cohort[[dropout_steps$first[[i]]]]
    column <- if_else(first %in% 1, cohort[[dropout_steps$later[[i]]]], NA_real_)
    column[unknown[[dropout_steps$step[[i]]]]] <- NA_real_
    cohort[[paste0("dropoff_ret_", dropout_steps$step[[i]])]] <- column
  }
  both <- !is.na(cohort$penta_1p_clean) & !is.na(cohort$measles_1p_clean) & !unknown$penta1_to_mcv1
  cohort$dropoff_mcv1_without_penta1 <- if_else(both, as.numeric(cohort$measles_1p_clean == 1 & cohort$penta_1p_clean == 0), NA_real_)
  cohort
}
dropout_fields <- c(paste0("dropoff_ret_", dropout_steps$step), "dropoff_mcv1_without_penta1")

# estimate each step in each sample, with the counts behind it
tracking <- read_rds("data/clean/rds_tracking.rds")
dropout <- imap(cohorts, \(cohort, sample) {
  cohort <- add_dropout_columns(cohort)
  estimate_fields(cohort, sample, dropout_fields, tracking = if (sample == "som_phase2_rds") tracking) |>
    mutate(
      sample = sample,
      n_domain = map_int(field, \(f) sum(!is.na(cohort[[f]]))),
      n_events = map_int(field, \(f) sum(cohort[[f]] == 1, na.rm = TRUE))
    )
}) |>
  list_rbind() |>
  left_join(sample_parts(study_samples), by = "sample") |>
  mutate(
    step = if_else(field == "dropoff_mcv1_without_penta1", "mcv1_without_penta1", str_remove(field, "dropoff_ret_")),
    metric = if_else(step == "mcv1_without_penta1", "mcv1_without_penta1", "retention")
  )

# take a count-based interval for a boundary arm
arm <- boundary_arm(dropout$estimate, dropout$standard_error, dropout$ci_lower, dropout$ci_upper, dropout$n_domain, dropout$n_events)
dropout <- mutate(dropout, boundary = arm$flag, ci_lower_used = arm$ci_lower, ci_upper_used = arm$ci_upper)

## Dropout compared between methods ----
# Phase 1 compares each method with the household survey and Phase 2 each pair
# of methods, with Newcombe's interval for the difference. The dropout
# difference is minus the retention difference.
# the comparisons, in report order
dropout_pairs <- tibble(
  country = rep(c("eth", "som", "eth", "som"), each = 3),
  phase = rep(c(1L, 1L, 2L, 2L), each = 3),
  method_a = c("cbe", "rds", "vbs", "cbe", "rds", "vbs", "cbe", "cbe", "rds", "cbe", "cbe", "rds"),
  method_b = c("hh", "hh", "hh", "hh", "hh", "hh", "rds", "vbs", "vbs", "rds", "vbs", "vbs")
)

# compare the two methods on each step
dropout_comparisons <-
  dropout_pairs |>
  inner_join(dropout, by = c("country", "phase", method_a = "method"), relationship = "many-to-many") |>
  inner_join(dropout, by = c("country", "phase", method_b = "method", "step", "metric"), suffix = c("_a", "_b"))
interval <- with(dropout_comparisons, newcombe_interval(estimate_a, ci_lower_used_a, ci_upper_used_a, estimate_b, ci_lower_used_b, ci_upper_used_b))
retention <- dropout_comparisons$metric == "retention"
dropout_comparisons <-
  dropout_comparisons |>
  mutate(
    diff_pp = if_else(retention, -1, 1) * 100 * (estimate_a - estimate_b),
    diff_lower_pp = if_else(retention, -100 * interval$upper, 100 * interval$lower),
    diff_upper_pp = if_else(retention, -100 * interval$lower, 100 * interval$upper),
    not_different = diff_lower_pp <= 0 & diff_upper_pp >= 0,
    boundary = boundary_a | boundary_b
  )

# mean absolute difference over the three steps, beside that in end-dose
# coverage of Penta3, OPV3 and MCV1
end_dose <-
  estimates |>
  filter(antigen %in% c("Penta3", "OPV3", "MCV1")) |>
  transmute(country, phase, method = str_to_lower(method), antigen, coverage = estimate)
dropout_summary <-
  dropout_comparisons |>
  filter(metric == "retention") |>
  summarise(mean_dropout = mean(abs(diff_pp)), not_different = sum(not_different), boundary = sum(boundary), .by = c(country, phase, method_a, method_b)) |>
  left_join(
    dropout_pairs |>
      inner_join(end_dose, by = c("country", "phase", method_a = "method"), relationship = "many-to-many") |>
      inner_join(end_dose, by = c("country", "phase", method_b = "method", "antigen"), suffix = c("_a", "_b")) |>
      summarise(mean_coverage = mean(abs(100 * (coverage_a - coverage_b))), .by = c(country, phase, method_a, method_b)),
    by = c("country", "phase", "method_a", "method_b")
  )

## Subgroups and places ----
# Coverage within groups of children (residence, child sex, caregiver
# education and, for HH and CBE, the place of the cluster; the settlement
# groups of Somalia Phase 2 CBE), each estimated by the sample's own estimator
# with the other children left out inside the design. A missing or declined
# answer is a group of its own.
cluster_places <- read_rds("data/clean/cluster_places.rds")

# the groups each variable can take, in display order
group_levels <- list(
  residence = c("Urban", "Rural", "Do not know", "Missing"),
  child_sex = c("Female", "Male", "Missing"),
  education = c("No formal education", "Primary", "Secondary", "Higher", "Missing"),
  settlement = c("Dhobley", "Afmadow town", "Other villages"),
  settlement_excl_cluster5 = c("Dhobley", "Afmadow town", "Other villages")
)
settlement_groups <- c(Dhobley = "Dhobley", `Afmadow town` = "Afmadow town", Dif = "Other villages", `Bilis Qooqaani` = "Other villages", Tabta = "Other villages")

# the group of each child for one variable
group_labels <- function(cohort, sample, variable) {
  parts <- sample_parts(sample)
  if (variable %in% c("place", "settlement", "settlement_excl_cluster5")) {
    cluster <- if (parts$method == "hh") cohort$cluster_id else cohort$cluster
    places <- filter(cluster_places, country == parts$country, phase == parts$phase, method == parts$method)
    place <- places$place_name[match(cluster, places$cluster_id)]
    stopifnot(!anyNA(place))
    if (variable == "place") return(place)
    group <- unname(settlement_groups[place])
    if (variable == "settlement_excl_cluster5") group[cluster == 5] <- NA
    return(group)
  }
  raw <- as.character(cohort[[variable]])
  if (variable == "education") raw <- unname(education_stages[raw])
  replace_na(raw, "Missing")
}

# the variables of each sample
sample_variables <- function(sample) {
  c(
    "residence", "child_sex", "education",
    if (str_detect(sample, "hh|cbe")) "place",
    if (sample == "som_phase2_cbe") c("settlement", "settlement_excl_cluster5")
  )
}

# estimate each outcome in each group of one variable; RDS seeds follow the
# column names, one column per group and outcome
group_outcomes <- c(opv_3 = "OPV3", penta_3 = "Penta3", measles_1 = "MCV1", zero_dose = "Zero-dose")
group_estimates <- function(cohort, sample, variable) {
  labels <- group_labels(cohort, sample, variable)
  levels <- if (variable == "place") sort(unique(labels)) else intersect(group_levels[[variable]], labels)
  columns <- expand_grid(j = seq_along(levels), outcome = names(group_outcomes)) |>
    mutate(level = levels[j], field = str_c("sg_", outcome, "_", j), antigen = group_outcomes[outcome])
  for (k in seq_len(nrow(columns))) {
    values <- cohort[[antigen_outcomes[[columns$antigen[[k]]]]]]
    cohort[[columns$field[[k]]]] <- if_else(labels %in% columns$level[[k]], values, NA_real_)
  }
  columns <- filter(columns, map_lgl(field, \(f) any(!is.na(cohort[[f]]))))
  cluster <- if (str_detect(sample, "hh")) cohort$cluster_id else if (str_detect(sample, "cbe")) cohort$cluster
  estimate_fields(cohort, sample, columns$field, tracking = if (sample == "som_phase2_rds") tracking) |>
    left_join(columns, by = "field") |>
    mutate(
      sample = sample, variable = variable,
      n_children = map_int(level, \(l) sum(labels %in% l)),
      n_events = map_int(field, \(f) sum(cohort[[f]] == 1, na.rm = TRUE)),
      n_clusters = if (is.null(cluster)) NA_integer_ else map_int(field, \(f) n_distinct(cluster[!is.na(cohort[[f]])]))
    )
}

# estimate every group of every sample
subgroups <-
  imap(cohorts, \(cohort, sample) map(sample_variables(sample), \(variable) group_estimates(cohort, sample, variable)) |> list_rbind()) |>
  list_rbind() |>
  left_join(sample_parts(study_samples), by = "sample")

## Suppression and count-based intervals ----
# An estimate is suppressed below 25 children with the outcome recorded or, for
# HH and CBE, fewer than 2 clusters, and cautioned at 25 to 49 children or 2 to
# 4 clusters. An estimate of 0% or 100% takes the Wilson interval on the
# unweighted counts. With fewer than 5 children with or without the outcome,
# or 2 to 4 clusters, the interval is the wider of the design-based and the
# count-based one; a weighted estimate outside the count-based interval keeps
# its own.
# flag each estimate and choose its interval
subgroups <-
  subgroups |>
  mutate(
    suppressed = n_observed < 25 | (!is.na(n_clusters) & n_clusters < 2),
    caution = !suppressed & (n_observed < 50 | (!is.na(n_clusters) & n_clusters < 5)),
    count_lower = map2_dbl(n_events, n_observed, \(x, n) count_wilson(x, n)[[1]]),
    count_upper = map2_dbl(n_events, n_observed, \(x, n) count_wilson(x, n)[[2]]),
    boundary = !suppressed & (estimate < 1e-9 | estimate > 1 - 1e-9),
    few = n_events < 5 | n_observed - n_events < 5 | (!is.na(n_clusters) & n_clusters >= 2 & n_clusters < 5),
    widened = !suppressed & !boundary & few & count_upper - count_lower > ci_upper - ci_lower &
      estimate >= count_lower & estimate <= count_upper,
    estimate = if_else(boundary, pmin(1, pmax(0, estimate)), estimate),
    ci_lower = if_else(boundary | widened, count_lower, ci_lower),
    ci_upper = if_else(boundary | widened, count_upper, ci_upper),
    across(c(estimate, ci_lower, ci_upper), \(x) if_else(suppressed, NA_real_, x))
  )

## Differences between groups ----
# Urban minus rural, female minus male, each pair of education levels and each
# pair of settlement groups, within a sample and outcome, with Newcombe's
# interval on the intervals shown; both groups must have an estimate.
# the pairs compared
group_pairs <- bind_rows(
  tibble(variable = "residence", level_a = "Urban", level_b = "Rural"),
  tibble(variable = "child_sex", level_a = "Female", level_b = "Male"),
  as_tibble(t(combn(group_levels$education[1:4], 2)), .name_repair = \(x) c("level_a", "level_b")) |> mutate(variable = "education"),
  as_tibble(t(combn(group_levels$settlement, 2)), .name_repair = \(x) c("level_a", "level_b")) |> mutate(variable = "settlement")
)

# compare each pair
comparisons <-
  group_pairs |>
  inner_join(filter(subgroups, !suppressed), by = c("variable", level_a = "level"), relationship = "many-to-many") |>
  inner_join(filter(subgroups, !suppressed), by = c("sample", "country", "phase", "method", "variable", "outcome", "antigen", level_b = "level"), suffix = c("_a", "_b"))
interval <- with(comparisons, newcombe_interval(estimate_a, ci_lower_a, ci_upper_a, estimate_b, ci_lower_b, ci_upper_b))
comparisons <-
  comparisons |>
  mutate(
    difference_pp = 100 * (estimate_a - estimate_b),
    distinguishable = !(interval$lower <= 0 & interval$upper >= 0),
    # higher coverage, or lower zero-dose prevalence, in the urban, female or
    # more educated group (level_b is the more educated in an education pair)
    better = if_else(variable == "education", -1, 1) * if_else(antigen == "Zero-dose", -1, 1) * difference_pp > 0
  )

# Tables ----

## Table 8 ----
# Percentage (95% CI) by country and outcome, one column per method.
# format each estimate and spread the methods across columns
table_08 <-
  phase1 |>
  mutate(
    Country = country_labels[country],
    Outcome = fct_recode(antigen, `Zero-dose prevalence` = "Zero-dose"),
    value = format_ci(estimate, ci_lower, ci_upper)
  ) |>
  arrange(country, Outcome, method) |>
  select(Country, Outcome, method, value) |>
  pivot_wider(names_from = method, values_from = value) |>
  mutate(Country = if_else(Outcome == "Penta3", Country, ""))

save_table(table_08, "output/tables/table_08.csv")

## Annex Table G1 ----
# Children, children with a recorded value n (%) per antigen, and the counts
# of the two reasons for a missing value, one row per sample.
# label each sample as country, phase and method
sample_label <- \(sample) {
  parts <- sample_parts(sample)
  str_c(country_labels[parts$country], ", ", parts$phase, ", ", method_labels[parts$method])
}

# format the counts and percentages
table_G01 <-
  missingness |>
  mutate(value = str_c(format_count(recorded), " (", format_number(100 * recorded / children), ")")) |>
  select(sample, children, antigen, value, not_known, doses_unknown) |>
  pivot_wider(names_from = antigen, values_from = value) |>
  transmute(
    `Country, phase, method` = sample_label(sample),
    Children = format_count(children),
    Penta1, Penta3, OPV1, OPV3, MCV1, MCV2,
    `Not known whether received (Penta / OPV / MCV)` = map_chr(not_known, \(x) str_c(x, collapse = " / ")),
    `Dose count not known (Penta3 / OPV3 / MCV2)` = map_chr(doses_unknown, \(x) str_c(x, collapse = " / "))
  )

save_table(table_G01, "output/tables/table_G01.csv")

## Annex Table G2 ----
# Visited clusters, children, the median and range of children per cluster,
# and the clusters with fewer than 10 children.
# summarise the children per cluster
table_G02 <-
  cluster_coverage |>
  filter(antigen == "Penta3") |>
  summarise(
    clusters = n(), children = sum(n_children),
    median = median(n_children), low = min(n_children), high = max(n_children), small = sum(n_children < 10),
    .by = sample
  ) |>
  transmute(
    `Country, phase, method` = sample_label(sample),
    Clusters = clusters,
    Children = children,
    `Children per cluster, median (range)` = str_c(median, " (", low, "-", high, ")"),
    `Clusters with fewer than 10 children` = small
  )

save_table(table_G02, "output/tables/table_G02.csv")

## Annex Table G3 ----
# Per sample and outcome: the sample estimate, the median and range of
# cluster coverage, the clusters outside the binomial limits, and the cluster
# whose omission changes the estimate most.
# signed change in pp to two decimals, with the cluster
loo_label <- \(delta, cluster) sprintf("%s%.2f pp (cluster %s)", if_else(delta < 0, "-", "+"), abs(delta), cluster)

# summarise each sample and outcome
table_G03 <-
  cluster_coverage |>
  filter(n_observed > 0) |>
  summarise(
    branch_estimate = first(branch_estimate),
    median = median(estimate), low = min(estimate), high = max(estimate),
    outside_95 = sum(outside_95), outside_998 = sum(outside_998), clusters = n(),
    largest = loo_label(loo_delta_pp[which.max(abs(loo_delta_pp))], cluster[which.max(abs(loo_delta_pp))]),
    .by = c(sample, antigen)
  ) |>
  transmute(
    `Country, phase, method` = sample_label(sample),
    Outcome = antigen,
    `Branch estimate` = str_c(format_number(100 * branch_estimate), "%"),
    `Cluster coverage, median (range)` = str_c(format_number(100 * median), "% (", format_number(100 * low), "-", format_number(100 * high), "%)"),
    `Clusters outside the 95% / 99.8% binomial limits` = str_c(outside_95, " / ", outside_998, " of ", clusters),
    `Largest leave-one-out change` = largest
  )

save_table(table_G03, "output/tables/table_G03.csv")

## Annex Table G4 ----
# Comparisons between groups by dimension and method, pooled over countries
# and phases: those made, those distinguishable from zero and, of these,
# those with better coverage in the urban, female or more educated group.
# count the comparisons
dimension_labels <- c(residence = "Residence", child_sex = "Child sex", education = "Caregiver education", settlement = "Settlement")
g4_counts <- \(data) summarise(data, comparisons = n(), better = sum(distinguishable & better), distinguishable = sum(distinguishable), .by = c(variable, method))
table_G04 <-
  bind_rows(
    g4_counts(comparisons),
    g4_counts(mutate(comparisons, method = "all"))
  ) |>
  bind_rows(summarise(comparisons, variable = "all", method = "all", comparisons = n(), distinguishable = sum(distinguishable))) |>
  arrange(match(variable, c(names(dimension_labels), "all")), match(method, c(method_order, "all"))) |>
  filter(!(variable == "settlement" & method == "all")) |>
  transmute(
    Dimension = c(dimension_labels, all = "All dimensions")[variable],
    Method = c(method_labels, all = "All methods")[method],
    Comparisons = comparisons,
    `Distinguishable from zero` = distinguishable,
    `Distinguishable, better coverage in the urban, female, or more educated group` = if_else(variable %in% c("settlement", "all"), "-", as.character(better))
  ) |>
  mutate(Dimension = if_else(Dimension == lag(Dimension, default = ""), "", Dimension))

save_table(table_G04, "output/tables/table_G04.csv")

## Annex Tables G5 to G8 ----
# Percentage (95% CI) for each group, with the children in it and, for HH and
# CBE, the clusters it spans. A cautioned estimate is in parentheses with its
# interval in square brackets; a dagger marks an interval widened to the
# count-based one and a double dagger an estimate of 0% or 100%. Rows whose
# four estimates are all suppressed are left out.
# format one estimate
format_group <- function(estimate, lower, upper, caution, widened, boundary) {
  mark <- if_else(boundary, "‡", if_else(widened, "†", ""))
  value <- if_else(
    caution,
    sprintf("(%s [%s-%s])%s", format_number(100 * estimate), format_number(100 * lower), format_number(100 * upper), mark),
    str_c(format_ci_hyphen(estimate, lower, upper), mark)
  )
  if_else(is.na(estimate), "-", value)
}

# one row per sample and group, the four outcomes as columns
group_rows <-
  subgroups |>
  mutate(value = format_group(estimate, ci_lower, ci_upper, caution, widened, boundary)) |>
  summarise(
    shown = any(!suppressed), children = max(n_children), clusters = max(n_clusters),
    Penta3 = value[antigen == "Penta3"], OPV3 = value[antigen == "OPV3"], MCV1 = value[antigen == "MCV1"], `Zero-dose` = value[antigen == "Zero-dose"],
    .by = c(sample, country, phase, method, variable, level)
  ) |>
  filter(shown) |>
  mutate(
    `Country, phase` = str_c(country_labels[country], ", ", phase),
    Method = method_labels[method],
    Children = children,
    Clusters = if_else(is.na(clusters), "-", as.character(clusters)),
    order = match(sample, study_samples)
  ) |>
  arrange(order)

# blank the repeated country and method labels
blank_repeats <- \(table) {
  table |>
    mutate(
      Method = if_else(`Country, phase` == lag(`Country, phase`, default = "") & Method == lag(Method, default = ""), "", Method),
      `Country, phase` = if_else(`Country, phase` == lag(`Country, phase`, default = ""), "", `Country, phase`)
    )
}
outcome_columns <- c("Penta3", "OPV3", "MCV1", "Zero-dose")

# residence and child sex
table_G05 <-
  group_rows |>
  filter(variable %in% c("residence", "child_sex")) |>
  arrange(order, match(variable, c("residence", "child_sex")), if_else(variable == "residence", match(level, group_levels$residence), match(level, group_levels$child_sex))) |>
  mutate(Characteristic = str_c(c(residence = "Residence: ", child_sex = "Child sex: ")[variable], level)) |>
  select(`Country, phase`, Method, Characteristic, Children, Clusters, all_of(outcome_columns)) |>
  blank_repeats()

save_table(table_G05, "output/tables/table_G05.csv")

# caregiver education
table_G06 <-
  group_rows |>
  filter(variable == "education") |>
  arrange(order, match(level, group_levels$education)) |>
  mutate(`Caregiver education` = level) |>
  select(`Country, phase`, Method, `Caregiver education`, Children, Clusters, all_of(outcome_columns)) |>
  blank_repeats()

save_table(table_G06, "output/tables/table_G06.csv")

# place, largest first; Somalia Phase 2 CBE is in Table G8
table_G07 <-
  group_rows |>
  filter(variable == "place", sample != "som_phase2_cbe") |>
  arrange(order, -children, level) |>
  mutate(Place = level) |>
  select(`Country, phase`, Method, Place, Children, Clusters, all_of(outcome_columns)) |>
  blank_repeats()

save_table(table_G07, "output/tables/table_G07.csv")

# settlement groups of Somalia Phase 2 CBE, then Afmadow town without cluster 5
table_G08 <-
  group_rows |>
  filter(variable == "settlement" | (variable == "settlement_excl_cluster5" & level == "Afmadow town")) |>
  mutate(
    Method = if_else(row_number() == 1, "CBE", ""),
    Place = if_else(variable == "settlement", level, "Afmadow town (excluding cluster 5)")
  ) |>
  select(Method, Place, Children, Clusters, all_of(outcome_columns))

save_table(table_G08, "output/tables/table_G08.csv")

## Annex Table G9 ----
# Difference in dropout on each step, and in the share with MCV1 but no
# Penta1 (pp, first method minus second, 95% CI); a dagger marks a comparison
# with a count-based interval.
# format the differences, one column per step
table_G09 <-
  dropout_comparisons |>
  mutate(
    `Country, phase` = str_c(country_labels[country], ", ", phase),
    Comparison = str_c(method_labels[method_a], " − ", method_labels[method_b]),
    value = str_c(format_difference_range(diff_pp, diff_lower_pp, diff_upper_pp), if_else(boundary, "†", "")),
    column = factor(step_labels[step], levels = step_labels)
  ) |>
  select(`Country, phase`, Comparison, column, value) |>
  pivot_wider(names_from = column, values_from = value)

save_table(table_G09, "output/tables/table_G09.csv")

## Annex Table G10 ----
# Mean absolute difference in dropout and in end-dose coverage between two
# methods, with the steps not distinguishable and those with a count-based
# interval.
# format the summary of each comparison
table_G10 <-
  dropout_summary |>
  transmute(
    `Country, phase` = str_c(country_labels[country], ", ", phase),
    Comparison = str_c(method_labels[method_a], " and ", method_labels[method_b]),
    `Mean absolute difference in dropout (pp)` = format_number(mean_dropout),
    `Mean absolute difference in end-dose coverage (pp)` = format_number(mean_coverage),
    `Steps not distinguishable (of 3)` = not_different,
    `Steps with a count-based interval` = boundary
  )

save_table(table_G10, "output/tables/table_G10.csv")

# Figures ----

## Figure 11 ----
# Each method's estimate and 95% CI against the household survey's estimate
# and 95% CI.
# create forest plot
figure_11 <- coverage_forest(phase1, phase = 1)
save_figure(figure_11, "output/figures/figure_11.png", height_cm = 11)

## Annex Figure G1 ----
# Coverage of the six antigens as a heat map, one row per method within each
# country and phase; labels are white on the darker half of the scale.
# arrange the estimates by country and phase
six_antigens <-
  estimates |>
  filter(antigen != "Zero-dose") |>
  arrange(country, phase, method, match(antigen, c("OPV1", "OPV3", "Penta1", "Penta3", "MCV1", "MCV2"))) |>
  mutate(
    method = factor(str_to_lower(method), levels = rev(method_order)),
    panel = factor(
      str_c(country_labels[country], ", ", phase_labels[as.character(phase)]),
      levels = as.vector(t(outer(country_labels, phase_labels, paste, sep = ", ")))
    ),
    antigen = fct_drop(antigen),
    pct = 100 * estimate
  )

# draw the heat map
figure_G1 <-
  ggplot(six_antigens, aes(antigen, method, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f", pct), colour = pct > 55), size = 2.6, family = report_font, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink[["primary"]])) +
  scale_fill_gradient(low = "#EEF3F8", high = "#1F3F66", limits = c(0, 100), name = "Coverage (%)") +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = method_labels) +
  labs(x = NULL, y = NULL) +
  theme_report() +
  theme(
    strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
    panel.border = element_blank(), panel.grid.major = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", colour = ink[["primary"]], size = 8.5)
  )

save_figure(figure_G1, "output/figures/figure_G1.png", height_cm = 15)

## Annex Figure G2 ----
# Share of each dose given on time, given late and without a recorded date,
# one bar per sample (country, phase, method) and one panel per dose, with a
# gap between country-phase groups.
# the three segments of each bar
dose_labels <- c(penta_1 = "Penta1", penta_2 = "Penta2", penta_3 = "Penta3", opv_1 = "OPV1", opv_2 = "OPV2", opv_3 = "OPV3", mcv_1 = "MCV1", mcv_2 = "MCV2")
timeliness_colours <- c(Timely = "#2F6FA3", Late = accent_colour, Undated = "#D5D5DB")
timeliness_bars <-
  bind_rows(
    mutate(timeliness, segment = "Timely", share = timely),
    mutate(timeliness, segment = "Late", share = late),
    mutate(timeliness, segment = "Undated", share = undated)
  ) |>
  mutate(
    segment = factor(segment, levels = rev(names(timeliness_colours))),
    dose = factor(dose_labels[dose], levels = dose_labels),
    group = str_c(country, " ", phase),
    row = str_c(country_labels[country], " P", phase, " ", method_labels[method])
  ) |>
  arrange(country, phase, match(method, method_order))

# rows top to bottom, with a blank row between groups
groups <- unique(timeliness_bars$group)
row_levels <- rev(unlist(map(seq_along(groups), \(i) {
  rows <- unique(timeliness_bars$row[timeliness_bars$group == groups[[i]]])
  if (i < length(groups)) c(rows, strrep(" ", i)) else rows
})))
timeliness_bars$row <- factor(timeliness_bars$row, levels = row_levels)

# draw the bars
figure_G2 <-
  ggplot(timeliness_bars, aes(share, row, fill = segment)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.5) +
  facet_wrap(vars(dose), ncol = 3) +
  scale_fill_manual(values = timeliness_colours, breaks = names(timeliness_colours), name = NULL) +
  scale_x_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
  scale_y_discrete(limits = row_levels, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_report(base_size = 8) +
  theme(
    strip.text = element_text(hjust = 0, face = "bold", size = 8.5), panel.spacing.x = unit(1.4, "lines"),
    panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7),
    plot.margin = margin(6, 14, 4, 6)
  )

save_figure(figure_G2, "output/figures/figure_G2.png", height_cm = 17)

## Annex Figure G3 ----
# Age at each dose in 15-day bins, CBE in Ethiopia: doses dated from the card
# above the axis and from recall below, with the scheduled age as a line;
# each row's axis is symmetric about zero.
# mirror the recall bins below the axis
offsets_plotted <-
  source_offsets |>
  arrange(match(sample, c("eth_phase1_cbe", "eth_phase2_cbe")), match(dose, names(dose_labels)), arm, bin) |>
  mutate(
    dose = factor(dose_labels[dose], levels = dose_labels),
    branch = factor(sample, levels = c("eth_phase1_cbe", "eth_phase2_cbe"), labels = c("Ethiopia Phase 1 CBE", "Ethiopia Phase 2 CBE")),
    signed = if_else(arm == "recall", -pct, pct),
    arm = factor(arm, levels = c("card", "recall"))
  )
targets <- doses |> transmute(dose = factor(dose_labels[dose], levels = dose_labels), target_days = target * 30.4375)
reach <- offsets_plotted |> summarise(reach = max(pct), .by = dose)
reach <- bind_rows(transmute(reach, dose, signed = reach), transmute(reach, dose, signed = -reach))

# draw the mirrored histograms
figure_G3 <-
  ggplot(offsets_plotted, aes(bin, signed, fill = arm)) +
  geom_col(position = "identity", width = 13) +
  geom_blank(data = reach, aes(x = 0, y = signed), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, colour = ink[["muted"]], linewidth = 0.3) +
  geom_vline(data = targets, aes(xintercept = target_days), colour = method_colours[["hh"]], linewidth = 0.4) +
  facet_grid(dose ~ branch, scales = "free_y", switch = "y") +
  # a few doses are recorded at implausible ages; the axis shows the plausible range
  coord_cartesian(xlim = c(0, 720)) +
  scale_fill_manual(values = c(card = "#2F6FA3", recall = accent_colour), labels = c(card = "Card (above the axis)", recall = "Recall (below)"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 720, 180), expand = c(0, 0)) +
  scale_y_continuous(labels = \(x) scales::label_percent(accuracy = 1)(abs(x)), expand = expansion(mult = c(0.04, 0.04))) +
  labs(x = "Age at dose (days, 15-day bins)", y = "Share of doses") +
  theme_report(base_size = 8) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1), axis.text.y = element_text(size = 6.5))

save_figure(figure_G3, "output/figures/figure_G3.png", height_cm = 20)

## Annex Figure G4 ----
# Coverage in each cluster, one row per sample: points sized by the children
# with the outcome recorded and hollow under 10, with the sample estimate and
# 95% CI over each row. Points are spread vertically by a fixed sequence,
# larger clusters first.
# place the cluster points
cluster_rows <- c(
  eth_phase1_hh = "Ethiopia Phase 1 HH", eth_phase1_cbe = "Ethiopia Phase 1 CBE", eth_phase2_cbe = "Ethiopia Phase 2 CBE",
  som_phase1_hh = "Somalia Phase 1 HH", som_phase1_cbe = "Somalia Phase 1 CBE", som_phase2_cbe = "Somalia Phase 2 CBE"
)
clusters_plotted <-
  cluster_coverage |>
  mutate(antigen = factor(antigen, levels = c("OPV3", "Penta3", "MCV1", "Zero-dose"))) |>
  arrange(match(sample, c("eth_phase1_hh", "eth_phase1_cbe", "som_phase1_hh", "som_phase1_cbe", "eth_phase2_cbe", "som_phase2_cbe")), antigen) |>
  filter(!is.na(estimate)) |>
  mutate(
    branch = factor(cluster_rows[sample], levels = rev(cluster_rows)),
    antigen = factor(antigen, levels = primary_outcomes)
  ) |>
  arrange(-n_observed) |>
  mutate(y = as.numeric(branch) + ((row_number() * 0.6180339887) %% 1 - 0.5) * 0.56)
branch_marks <- distinct(clusters_plotted, sample, antigen, .keep_all = TRUE)

# draw the clusters and the sample estimates
figure_G4 <-
  ggplot() +
  geom_point(data = clusters_plotted, aes(100 * estimate, y, size = n_observed, colour = method, shape = small_n, alpha = small_n), stroke = 0.5) +
  geom_linerange(data = branch_marks, aes(y = as.numeric(branch), xmin = 100 * branch_ci_lower, xmax = 100 * branch_ci_upper), linewidth = 1, colour = ink[["primary"]]) +
  geom_point(data = branch_marks, aes(100 * branch_estimate, as.numeric(branch)), shape = 23, fill = "white", colour = ink[["primary"]], size = 2.4, stroke = 0.8) +
  facet_wrap(vars(antigen), ncol = 2) +
  scale_colour_manual(values = method_colours[c("hh", "cbe")], breaks = c("hh", "cbe"), labels = method_labels[c("hh", "cbe")], name = "Method") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = "10 or more children", `TRUE` = "Fewer than 10 children"), name = "Cluster") +
  scale_alpha_manual(values = c(`FALSE` = 0.5, `TRUE` = 0.9), guide = "none") +
  scale_size_area(max_size = 3, name = "Children observed") +
  scale_y_continuous(breaks = seq_along(cluster_rows), labels = rev(cluster_rows), expand = c(0, 0), limits = c(0.5, length(cluster_rows) + 0.5)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 2.6, alpha = 0.9, shape = 16)),
    shape = guide_legend(order = 3, override.aes = list(size = 2.6, colour = ink[["secondary"]], alpha = 0.9)),
    size = guide_legend(order = 2, override.aes = list(colour = ink[["secondary"]], alpha = 0.7))
  ) +
  labs(x = "Cluster coverage (%)", y = NULL) +
  theme_report(base_size = 8.5) +
  theme(
    legend.box = "horizontal", legend.title.position = "top", legend.spacing.x = unit(0.5, "cm"),
    panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7.5)
  )

save_figure(figure_G4, "output/figures/figure_G4.png", height_cm = 14)

## Annex Figure G5 ----
# Dropout on each step with its 95% CI, methods side by side within each step;
# rows are countries and columns phases, and an open point marks a count-based
# interval.
# place each method's point within its step
dropout_points <-
  dropout |>
  filter(metric == "retention") |>
  # rounded, so a retention of 1 that the survey package returns as 1 + 2e-16 plots at 0
  mutate(
    dropout = round(100 * (1 - estimate), 10), lower = round(100 * (1 - ci_upper_used), 10), upper = round(100 * (1 - ci_lower_used), 10),
    fill_colour = if_else(boundary, "white", method_colours[method]),
    step_index = match(step, dropout_steps$step)
  ) |>
  group_by(country, phase) |>
  mutate(
    method = factor(method, levels = method_order),
    offset = (as.integer(fct_drop(method)) - (n_distinct(method) + 1) / 2) * 0.17,
    x = step_index + offset
  ) |>
  ungroup() |>
  mutate(
    country_f = factor(country_labels[country], levels = country_labels),
    panel_phase = factor(phase_labels[as.character(phase)], levels = phase_labels)
  ) |>
  arrange(country, phase, method, step_index)

# draw the dropout points, Penta1 to Penta3 on top
figure_G5 <-
  ggplot(dropout_points, aes(x = x, colour = method)) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.55, show.legend = FALSE) +
  geom_point(aes(y = dropout, fill = fill_colour, shape = factor(boundary, levels = c(FALSE, TRUE))), size = 2.4, stroke = 0.8) +
  facet_grid(country_f ~ panel_phase, switch = "y") +
  scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 21), breaks = "TRUE", labels = "Count-based arm (open symbol)", name = NULL,
    guide = guide_legend(order = 2, override.aes = list(fill = "white", colour = ink[["secondary"]], size = 2.4, stroke = 0.8))
  ) +
  method_scale(method_order) +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = 1:3, labels = c("Penta1 to\nPenta3", "OPV1 to\nOPV3", "Penta1 to\nMCV1"),
    limits = c(3.5, 0.5), trans = "reverse", expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0.03, 0.03)) +
  coord_flip() +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 21, fill = unname(method_colours[method_order]), size = 2.4))) +
  labs(x = NULL, y = "Dropout (%)") +
  theme_report() +
  theme(panel.grid.major.y = element_blank(), panel.spacing.y = unit(2, "lines"))

save_figure(country_divider(figure_G5, 12), "output/figures/figure_G5.png", height_cm = 12)
