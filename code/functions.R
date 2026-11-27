# functions.R
# Helpers shared by the cleaning and analysis scripts

# Study constants ----

# the fourteen samples (country, phase and method), in report order
study_samples <- c(
  "eth_phase1_hh", "eth_phase1_cbe", "eth_phase1_rds", "eth_phase1_vbs",
  "som_phase1_hh", "som_phase1_cbe", "som_phase1_rds", "som_phase1_vbs",
  "eth_phase2_cbe", "eth_phase2_rds", "eth_phase2_vbs",
  "som_phase2_cbe", "som_phase2_rds", "som_phase2_vbs"
)

# children aged 12-23 months in each RDS study area, from the protocol (HRSL
# under-fives / 5 in Ethiopia, WorldPop 2025 under-ones in Somalia)
population_size <- c(eth_phase1_rds = 5973, som_phase1_rds = 12093, eth_phase2_rds = 5699, som_phase2_rds = 11274)

study_seed <- 20260902

# outcome fields, and the three primary antigens and zero-dose in display order
antigen_outcomes <- c(
  "Penta1" = "penta_1p_clean", "Penta3" = "penta_3p_clean",
  "OPV1" = "opv_1p_clean", "OPV3" = "opv_3p_clean",
  "MCV1" = "measles_1p_clean", "MCV2" = "measles_2p_clean",
  "Zero-dose" = "zero_dose"
)
primary_outcomes <- c("Penta3", "OPV3", "MCV1", "Zero-dose")

antigen_fields <- c("opv_1p_clean", "opv_3p_clean", "penta_1p_clean", "penta_3p_clean", "measles_1p_clean", "measles_2p_clean")

# Labels ----

## Participation status ----
# A respondent who agreed and reported at least one eligible child (own or in
# the network module) is enrolled. A missing consent answer counts as enrolled
# only when the interview went on past consent.
classify_consent <- function(main, progressed) {
  consent <- str_to_lower(str_trim(as.character(main$consent)))
  own_children <- coalesce(as.numeric(as.character(main$hh_child_num)), 0)
  network_children <- coalesce(as.numeric(as.character(main$nsum_child_total)), 0)
  answered <-
    main |>
    select(any_of(c("hh_child_num", "nsum_child_total", "respondent_age", "respondent_sex", "relation_child"))) |>
    map(\(x) !is.na(x)) |>
    reduce(`|`)

  case_when(
    consent == "yes" & own_children == 0 & network_children == 0 ~ "Ineligible",
    consent == "yes" ~ "Enrolled",
    consent == "no" ~ "Refused",
    is.na(consent) & (answered | progressed) ~ "Enrolled",
    is.na(consent) ~ "Ineligible"
  )
}

## Caregiver relation ----
classify_relation <- function(relation) {
  relation <- str_trim(as.character(relation))
  case_when(
    relation == "Primary female caregiver" ~ "Female caregiver",
    relation == "Primary male caregiver (e.g. widower, single father, or other male relative)" ~ "Male caregiver",
    relation == "Other caregiver (e.g. grandmother, aunt, sibling, or non-parental guardian)" ~ "Other caregiver",
    relation == "Head of household" ~ "Head of household"
  ) |>
    factor(levels = c("Female caregiver", "Head of household", "Male caregiver", "Other caregiver"), ordered = TRUE)
}

## Dose counts ----
# Number of doses as asked on the form, and the labels used in the analysis.
dose_answers <- list(opv = c("1", "2", "3", "4 (or more)"), penta = c("1", "2", "3"), measles = c("1", "2 (or more)"))
dose_labels <- list(
  opv = c("0", "1", "2", "3", "4+", "Doses Unk", "Vacc Unk"),
  penta = c("0", "1", "2", "3", "Doses Unk", "Vacc Unk"),
  measles = c("0", "1", "2+", "Doses Unk", "Vacc Unk")
)

classify_dose_count <- function(ever, n, family) {
  ever <- as.character(ever)
  n <- as.character(n)
  answered <- match(n, dose_answers[[family]])
  case_when(
    ever == "No" | is.na(ever) ~ "0",
    !is.na(answered) ~ dose_labels[[family]][answered + 1],
    n == "Do not know" ~ "Doses Unk",
    ever == "Do not know" ~ "Vacc Unk"
  ) |>
    factor(levels = dose_labels[[family]], ordered = TRUE)
}

# Dose-count labels and the threshold indicators built on them, for children
# reported directly (child) or through the network module (network).
classify_doses <- function(data, type = "child") {
  if (type == "network") {
    data$opv_ever <- data$nsum_opv_ever
    data$penta_ever <- data$nsum_penta_ever
    data$measles_ever <- data$nsum_measles_ever
    data$opv_n_raw <- data$nsum_opv_n
    data$penta_n_raw <- data$nsum_penta_n
    data$measles_n_raw <- data$nsum_measles_n
  } else {
    data$opv_n_raw <- data$opv_n
    data$penta_n_raw <- data$penta_n
    data$measles_n_raw <- data$measles_n
  }
  data$opv_n <- classify_dose_count(data$opv_ever, data$opv_n_raw, "opv")
  data$penta_n <- classify_dose_count(data$penta_ever, data$penta_n_raw, "penta")
  data$measles_n <- classify_dose_count(data$measles_ever, data$measles_n_raw, "measles")

  data |>
    mutate(
      opv_1p = as.numeric(opv_n %in% c("1", "2", "3", "4+")),
      opv_2p = as.numeric(opv_n %in% c("2", "3", "4+")),
      opv_3p = as.numeric(opv_n %in% c("3", "4+")),
      opv_4p = as.numeric(opv_n %in% "4+"),
      opv_n_dnk = as.numeric(opv_n %in% "Doses Unk"),
      opv_dnk = as.numeric(opv_n %in% "Vacc Unk"),
      penta_1p = as.numeric(penta_n %in% c("1", "2", "3")),
      penta_2p = as.numeric(penta_n %in% c("2", "3")),
      penta_3p = as.numeric(penta_n %in% "3"),
      penta_n_dnk = as.numeric(penta_n %in% "Doses Unk"),
      penta_dnk = as.numeric(penta_n %in% "Vacc Unk"),
      measles_1p = as.numeric(measles_n %in% c("1", "2+")),
      measles_2p = as.numeric(measles_n %in% "2+"),
      measles_n_dnk = as.numeric(measles_n %in% "Doses Unk"),
      measles_dnk = as.numeric(measles_n %in% "Vacc Unk")
    )
}

# Antigens ----

## Antigen indicators ----
# First doses come from the "ever received" answer; later doses also need the
# number of doses (any OPV dose counts, the birth dose included).
first_dose <- function(ever) {
  ever <- as.character(ever)
  case_when(ever == "Yes" ~ 1, ever == "No" ~ 0, .default = NA_real_)
}

later_dose <- function(ever, n, target, below) {
  ever <- as.character(ever)
  n <- as.character(n)
  case_when(
    ever == "No" ~ 0,
    ever == "Yes" & n %in% target ~ 1,
    ever == "Yes" & n %in% below ~ 0,
    .default = NA_real_
  )
}

derive_antigens <- function(data) {
  data |>
    mutate(
      opv_1p_clean = first_dose(opv_ever),
      opv_3p_clean = later_dose(opv_ever, opv_n, target = c("3", "4+"), below = c("1", "2")),
      penta_1p_clean = first_dose(penta_ever),
      penta_3p_clean = later_dose(penta_ever, penta_n, target = "3", below = c("1", "2")),
      measles_1p_clean = first_dose(measles_ever),
      measles_2p_clean = later_dose(measles_ever, measles_n, target = "2+", below = "1")
    )
}

## Never-vaccinated children ----
# On the household form, a caregiver without a card who said the child had
# never been vaccinated skipped the dose questions; the child counts as
# unvaccinated for every antigen.
complete_never_vaccinated <- function(data) {
  never <- as.character(data$reason_no_card) %in% "Child has never been vaccinated"
  for (field in antigen_fields) {
    data[[field]][never] <- 0
  }
  data
}

## Source of each answer ----
# The source question (card or recall) was asked only when a dose was reported,
# and on the household form only when the caregiver had a card. A reported
# non-receipt takes the card as its source for card holders and recall
# otherwise; on the household form, any answer from a caregiver without a card
# is recall.
code_answer_source <- function(data, holder_field, household = FALSE) {
  holder <- as.character(data[[holder_field]])
  no_card <- !is.na(holder) & holder == "No"
  reason_never <- household & as.character(data$reason_no_card) %in% "Child has never been vaccinated"

  for (family in c("penta", "opv", "measles")) {
    ever <- as.character(data[[paste0(family, "_ever")]])
    field <- paste0(family, "_source")
    source <- as.character(data[[field]])
    not_received <- !is.na(ever) & ever == "No"
    source[not_received & !is.na(holder) & holder == "Yes"] <- "Vaccination card"
    source[not_received & no_card] <- "Respondent recall"
    if (household) {
      source[!is.na(ever) & ever == "Yes" & no_card & is.na(source)] <- "Respondent recall"
      source[is.na(ever) & no_card & reason_never & is.na(source)] <- "Respondent recall"
    }
    data[[field]] <- factor(source, levels = union(levels(data[[field]]), c("Vaccination card", "Respondent recall")))
  }
  data
}

# Venue-based sampling ----

## Interview venue ----
classify_venue <- function(venue, country) {
  venue <- str_trim(as.character(venue))
  if (country == "eth") {
    case_when(
      venue == "Worabe referral hospital" ~ "Worabe hospital",
      venue == "Worabe bus station" ~ "Worabe bus",
      venue %in% c("Inseno bus station", "In Inseno") ~ "Inseno bus",
      venue == "Inseno market" ~ "Inseno market",
      .default = "Other"
    )
  } else {
    case_when(
      venue == "Galkayo North MCH" ~ "G North MCH",
      venue == "Adado NIRA Center" ~ "Adado NIRA",
      venue == "Gaalkayo North market" ~ "G North market",
      .default = "Other"
    )
  }
}

## Earlier ESTIMATE interviews ----
# Where a Phase 1 VBS respondent had been interviewed before: at a venue, by
# phone (RDS) or at home (HH or CBE).
classify_previous_survey <- function(main, country) {
  venues <- if (country == "eth") {
    c("where_previous_survey_healthfac", "where_previous_survey_bus", "where_previous_survey_inseno", "where_previous_survey_market")
  } else {
    c("where_previous_survey_healthfac", "where_previous_survey_adado_nid", "where_previous_survey_market")
  }
  as_number <- \(x) as.numeric(as.character(x))
  main$previous_survey_vbs <- main[venues] |> map(as_number) |> as_tibble() |> rowSums(na.rm = TRUE) > 0
  main$previous_survey_rds <- as_number(main$where_previous_survey_phone) %in% 1
  main$previous_survey_hh <- as_number(main$where_previous_survey_home) %in% 1
  main
}

## Participation status ----
classify_participant_status <- function(hh_child_num) {
  n <- as.numeric(as.character(hh_child_num))
  case_when(
    n == 0 ~ "Enrolled without Eligible Child",
    n >= 1 ~ "Enrolled with Eligible Child",
    is.na(n) ~ "Not Enrolled"
  )
}

## Mobility weights ----
# Inverse of the number of visits to the venue town in the past 4 weeks (bucket
# midpoints, "0" scored as "1-4" since the respondent was met in town), scaled
# to a mean of 1. Respondents who preferred not to answer keep a weight of 1.
mobility_weight <- function(travel, multiplier = 1) {
  travel <- as.character(travel)
  score <- case_when(travel %in% c("0", "1-4") ~ 2.5, travel == "5-9" ~ 7, travel == "10+" ~ 10) * multiplier
  weight <- (1 / score) / mean(1 / score, na.rm = TRUE)
  replace_na(weight, 1)
}

# the travel answer for the town where each respondent was met
vbs_travel <- function(cohort, sample) {
  switch(
    sample,
    eth_phase1_vbs = if_else(str_detect(cohort$int_location, "Worabe"), as.character(cohort$worabe_travel), as.character(cohort$inseno_travel)),
    som_phase1_vbs = if_else(str_detect(cohort$int_location, "G North"), as.character(cohort$ngalkayo_travel), as.character(cohort$adado_travel)),
    eth_phase2_vbs = as.character(cohort$bahir_dar_travel),
    som_phase2_vbs = as.character(cohort$kismayo_travel)
  )
}

# Somalia Phase 1 asked about Galkayo North per week and Adado per 4 weeks, so
# Galkayo North answers are put on a 4-week basis
travel_multiplier <- function(cohort, sample) {
  if (sample == "som_phase1_vbs") if_else(str_detect(cohort$int_location, "G North"), 4, 1) else 1
}

vbs_weight <- function(cohort, sample) mobility_weight(vbs_travel(cohort, sample), travel_multiplier(cohort, sample))

# Cohorts ----

## Directly reported children ----
# Enrolled respondents' children after the analytic exclusions, with zero-dose
# as no Penta1; every CBE record is a child, and VBS keeps the caregivers of an
# eligible child.
direct_cohort <- function(frame, sample) {
  cohort <- frame |> filter(analytic) |> mutate(zero_dose = 1 - penta_1p_clean)
  if (str_detect(sample, "hh|rds")) cohort <- filter(cohort, consent == "Enrolled")
  if (str_detect(sample, "vbs")) cohort <- filter(cohort, consent == "Enrolled", participant_status == "Enrolled with Eligible Child")
  cohort
}

## Directly reporting respondents ----
# One row per respondent in the cohort (a CBE record is its own respondent),
# with the respondent key in `key`.
analytic_respondents <- function(frame, sample) {
  cohort <- direct_cohort(frame, sample)
  cohort$key <- if (str_detect(sample, "cbe")) cohort$record_id else cohort$record_ix
  distinct(cohort, key, .keep_all = TRUE)
}

## Sample labels ----
# Country, phase and method of a sample name such as "eth_phase1_hh".
sample_parts <- function(sample) {
  tibble(
    sample = sample,
    country = str_sub(sample, 1, 3),
    phase = as.integer(str_extract(sample, "(?<=phase)\\d")),
    method = str_extract(sample, "[a-z]+$")
  )
}

# Respondent characteristics ----

## Education and household wealth ----
# Education grouped into four stages, and household wealth as "High asset" when
# the respondent owns more than half of the assets scored.
education_stages <- c(
  "None" = "No formal education", "Read and write (no formal education)" = "No formal education",
  "Religious only" = "No formal education", "Primary" = "Primary", "Some primary" = "Primary",
  "Completed primary" = "Primary", "Secondary" = "Secondary", "Some secondary" = "Secondary",
  "Completed secondary" = "Secondary", "More than secondary" = "Higher", "Higher education" = "Higher"
)

# household assets asked about in at least one sample
asset_candidates <- c("owns_radio", "owns_tv", "owns_fridge", "owns_computer", "owns_aircon", "owns_mobile", "mattress")

# score household assets on the items given
wealth_group <- function(d, assets) {
  if (length(assets) == 0) return(rep(NA_character_, nrow(d)))
  scores <- map(assets, \(a) case_when(d[[a]] == "Yes" ~ 1, d[[a]] == "No" ~ 0)) |> reduce(cbind) |> matrix(nrow = nrow(d))
  answered <- rowSums(!is.na(scores))
  case_when(answered == 0 ~ NA, rowSums(scores, na.rm = TRUE) > length(assets) / 2 ~ "High asset", .default = "Low asset")
}

# the assets asked about in every sample of a sample's country and phase
common_assets <- function(cohorts, sample) {
  group <- str_sub(sample, 1, 10)
  columns <- reduce(map(cohorts[str_starts(names(cohorts), group)], names), intersect)
  intersect(asset_candidates, columns)
}

# Estimation ----

## Seeds ----
# Each sample and outcome draws from its own seed, derived from the study seed
# by hashing.
derived_seed <- function(seed, key) {
  hash <- digest::digest(paste(seed, key, sep = "::"), algo = "sha256", serialize = FALSE)
  as.integer(strtoi(substr(hash, 1, 7), base = 16)) + 1L
}

sample_seed <- function(sample) derived_seed(study_seed, str_replace(sample, "_phase(\\d)_", "/phase_\\1/"))

## Design-based proportion ----
# Weighted mean with a Wilson interval at the design's effective sample size
# and degrees of freedom; a zero standard error gives a point interval.
survey_proportion <- function(design, field) {
  formula <- as.formula(paste0("~", field))
  mean <- svymean(formula, design, na.rm = TRUE)
  estimate <- unname(coef(mean)[[1]])
  standard_error <- unname(SE(mean)[[1]])
  if (standard_error == 0) {
    return(tibble(estimate, standard_error, ci_lower = estimate, ci_upper = estimate))
  }
  interval <- attr(svyciprop(formula, design, method = "wilson", na.rm = TRUE, df = degf(design)), "ci")
  tibble(estimate, standard_error, ci_lower = unname(interval[[1]]), ci_upper = unname(interval[[2]]))
}

## Wilson interval from an estimate and its standard error ----
wilson_interval <- function(estimate, standard_error, z = 1.96) {
  if (standard_error == 0) return(c(estimate, estimate))
  n_eff <- estimate * (1 - estimate) / standard_error^2
  centre <- estimate + z^2 / (2 * n_eff)
  half <- z / (2 * n_eff) * sqrt(4 * n_eff * estimate * (1 - estimate) + z^2)
  c(centre - half, centre + half) / (1 + z^2 / n_eff)
}

## RDS recruitment network ----
# One node per enrolled respondent, with the recruiter (0 for a seed), the
# degree (eligible children known, at least the number recruited and at least
# 1) and the mean of each outcome over the respondent's children. With a
# tracking sheet, a referral code that matches no interviewed respondent is
# added as a node without data when the sheet records that person, and treated
# as the start of a chain otherwise.
rds_network <- function(cohort, outcomes, tracking = NULL) {
  nodes <-
    cohort |>
    distinct(date, id_participant_calc, id_referring, seed_participant, nsum_child_total) |>
    mutate(id_participant_calc = as.character(id_participant_calc), id_referring = as.character(id_referring))
  nodes$id_referring[as.character(nodes$seed_participant) == "Yes"] <- "0"
  nodes$degree_raw <- suppressWarnings(as.integer(nodes$nsum_child_total))

  for (field in outcomes) {
    nodes[[field]] <- node_mean(cohort, nodes$id_participant_calc, field)
  }

  if (!is.null(tracking)) {
    tracking <- mutate(tracking, across(c(id_participant, id_referring), as.character))
    for (pass in 1:10) {
      orphans <- setdiff(na.omit(nodes$id_referring), c(nodes$id_participant_calc, "0"))
      if (length(orphans) == 0) break
      tracked <- filter(tracking, id_participant %in% orphans)
      if (nrow(tracked) > 0) {
        nodes <- bind_rows(nodes, tibble(
          date = tracked$date,
          id_participant_calc = tracked$id_participant,
          id_referring = coalesce(tracked$id_referring, "0"),
          seed_participant = "No"
        ))
      }
      nodes$id_referring[nodes$id_referring %in% setdiff(orphans, tracked$id_participant)] <- "0"
    }
  }

  recruits <- table(nodes$id_referring[!is.na(nodes$id_referring) & !nodes$id_referring %in% c("", "NA")])
  nodes$n_recruited <- coalesce(as.integer(recruits[nodes$id_participant_calc]), 0L)
  nodes$degree <- pmax(1L, nodes$n_recruited, nodes$degree_raw, na.rm = TRUE)
  nodes
}

# mean outcome over each node's children, and the number of children with an
# observed outcome
node_mean <- function(cohort, node_ids, field) {
  means <- tapply(cohort[[field]], as.character(cohort$id_participant_calc), mean, na.rm = TRUE)
  out <- as.vector(means[node_ids])
  out[is.nan(out)] <- NA_real_
  out
}

node_children <- function(cohort, node_ids, field) {
  observed <- table(as.character(cohort$id_participant_calc)[!is.na(cohort[[field]])])
  coalesce(as.integer(observed[node_ids]), 0L)
}

## RDS successive-sampling estimate ----
# Gile's successive-sampling weights for the respondent network, applied to
# each respondent's children: the estimate is the weighted share of children
# vaccinated, with a second-order empirical-likelihood standard error (summed
# over all pairs of respondents) and a Wilson interval.
rds_child_mean <- function(weight, outcome, population) {
  keep <- !is.na(weight) & !is.na(outcome)
  w <- weight[keep]
  y <- outcome[keep]
  estimate <- sum(y * w) / sum(w)
  wij <- outer(w, w)
  wij[wij > 10 * median(wij)] <- 10 * median(wij)
  off_diagonal <- row(wij) != col(wij)
  yij <- 0.5 * outer(y, y, "+")
  centre <- sum(yij * wij * off_diagonal) / sum(wij * off_diagonal)
  gi <- rowSums(wij^2 * (yij - centre) * off_diagonal)
  vi <- rowSums(wij^2 * off_diagonal)
  variance <- 4 * sum(gi^2) / sum(vi)^2
  list(estimate = estimate, standard_error = sqrt((population - length(outcome)) * variance / (population - 1)))
}

rds_proportion <- function(design, cohort, field, population, seed) {
  fit <- withr::with_seed(seed, RDS.SS.estimates(design, outcome.variable = field, N = population))
  children <- node_children(cohort, as.character(design$id_participant_calc), field)
  result <- rds_child_mean(as.numeric(fit$weights) * children, as.numeric(design[[field]]), population)
  interval <- wilson_interval(result$estimate, result$standard_error)
  tibble(
    estimate = result$estimate, standard_error = result$standard_error,
    ci_lower = interval[[1]], ci_upper = interval[[2]], n_observed = sum(children)
  )
}

## Estimates for any outcome column ----
# Each sample's own estimator (as in 02_estimates.R) applied to derived 0/1/NA
# columns of its cohort, such as an outcome within a subgroup or a later dose
# among children with the first. A child whose value is missing is left out
# inside the design, so the estimate is a domain estimate. RDS draws a seed for
# each column. An estimate the estimator cannot compute (a handful of children)
# is returned as missing.
estimate_fields <- function(cohort, sample, fields, tracking = NULL) {
  method <- str_extract(sample, "[a-z]+$")
  safe <- \(f) tryCatch(f(), error = \(e) tibble(estimate = NA_real_, standard_error = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  if (method == "rds") {
    nodes <- rds_network(cohort, fields, tracking = tracking)
    design <- as.rds.data.frame(
      nodes,
      id = "id_participant_calc", recruiter.id = "id_referring",
      network.size = "degree", population.size = population_size[[sample]]
    )
    return(map(fields, \(field) {
      seed <- derived_seed(sample_seed(sample), paste0("rds_outcome/", field))
      safe(\() rds_proportion(design, cohort, field, population_size[[sample]], seed)) |>
        mutate(n_observed = sum(!is.na(cohort[[field]])))
    }) |>
      list_rbind() |>
      mutate(field = fields, .before = 1))
  }
  if (method == "vbs") {
    cohort$weight <- vbs_weight(cohort, sample)
    design <- svydesign(id = ~record_ix, weights = ~weight, data = cohort)
  } else {
    if (method == "cbe") cohort <- mutate(cohort, cluster_id = cluster, record_ix = record_id)
    design <- svydesign(id = ~ cluster_id + record_ix, weights = ~weight, nest = TRUE, data = mutate(cohort, weight = 1))
  }
  map(fields, \(field) safe(\() survey_proportion(design, field)) |> mutate(n_observed = sum(!is.na(cohort[[field]])))) |>
    list_rbind() |>
    mutate(field = fields, .before = 1)
}

## RDS child weights ----
# The successive-sampling weight of each child's respondent for one outcome,
# drawn with the outcome's own seed as in the primary estimate.
rds_child_weights <- function(design, cohort, field, population, seed) {
  fit <- withr::with_seed(seed, RDS.SS.estimates(design, outcome.variable = field, N = population))
  stopifnot(length(fit$weights) == nrow(design))
  as.numeric(fit$weights)[match(as.character(cohort$id_participant_calc), as.character(design$id_participant_calc))]
}

## RDS recruitment chains ----
# The first node of each node's chain, following referrals back to a seed or
# to a recruiter who is not in the network.
chain_root <- function(nodes) {
  ids <- as.character(nodes$id_participant_calc)
  parent <- as.character(nodes$id_referring)
  map_chr(ids, \(id) {
    seen <- character(0)
    repeat {
      if (id %in% seen) break
      seen <- c(seen, id)
      up <- parent[match(id, ids)]
      if (is.na(up) || up == "0" || !up %in% ids) break
      id <- up
    }
    id
  })
}

# Comparisons ----

## Unweighted counts ----
# Children with an observed outcome (n_domain) and with the outcome (n_events)
# in each sample's cohort, for the count-based intervals below.
event_counts <- function(frames) {
  samples <- setdiff(names(frames), "som_phase2_rds_all_enumerators")
  map(samples, \(sample) {
    cohort <- direct_cohort(frames[[sample]], sample)
    map(names(antigen_outcomes), \(antigen) {
      values <- cohort[[antigen_outcomes[[antigen]]]]
      tibble(sample = sample, antigen = antigen, n_domain = sum(!is.na(values)), n_events = sum(values == 1, na.rm = TRUE))
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}

## Newcombe hybrid-score interval ----
# Interval for the difference of two independent proportions, built from each
# proportion's own interval.
newcombe_interval <- function(p1, l1, u1, p2, l2, u2) {
  difference <- p1 - p2
  list(
    lower = difference - sqrt((p1 - l1)^2 + (u2 - p2)^2),
    upper = difference + sqrt((u1 - p1)^2 + (p2 - l2)^2)
  )
}

## Count-based Wilson interval ----
count_wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  half <- (z / denominator) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(max(0, centre - half), min(1, centre + half))
}

## Boundary arms ----
# An estimate whose standard error is zero, that sits at 0 or 1, or that rests
# on fewer than five children with or without the outcome enters a comparison
# with the count-based Wilson interval on its unweighted counts (widened to
# contain a weighted estimate); every other estimate keeps its own interval.
boundary_arm <- function(estimate, standard_error, ci_lower, ci_upper, n_domain, n_events) {
  flag <- standard_error < 1e-9 | estimate < 1e-9 | estimate > 1 - 1e-9 | n_events < 5 | (n_domain - n_events) < 5
  for (i in which(flag)) {
    interval <- count_wilson(n_events[[i]], n_domain[[i]])
    ci_lower[[i]] <- max(0, min(interval[[1]], estimate[[i]]))
    ci_upper[[i]] <- min(1, max(interval[[2]], estimate[[i]]))
  }
  list(flag = flag, ci_lower = ci_lower, ci_upper = ci_upper)
}

## Phase 1 comparison with the household survey ----
# One row per country, alternative method and outcome: the difference from the
# household survey, the width of the method's interval, and whether Newcombe's
# interval for the difference contains zero (on the count-based interval for a
# boundary arm). `counts` gives each sample's unweighted counts; the household
# rows of `estimates` are the reference.
phase1_performance <- function(estimates, counts) {
  phase1 <-
    estimates |>
    filter(phase == 1, antigen %in% primary_outcomes) |>
    mutate(antigen = as.character(antigen), method = as.character(method)) |>
    left_join(counts, by = c("sample", "antigen"))
  reference <-
    phase1 |>
    filter(method == "HH") |>
    select(
      country, antigen, ref_estimate = estimate, ref_se = standard_error, ref_lower = ci_lower, ref_upper = ci_upper,
      ref_n_domain = n_domain, ref_n_events = n_events
    )

  performance <-
    phase1 |>
    filter(method != "HH") |>
    left_join(reference, by = c("country", "antigen"))
  arm <- boundary_arm(performance$estimate, performance$standard_error, performance$ci_lower, performance$ci_upper, performance$n_domain, performance$n_events)
  ref <- boundary_arm(performance$ref_estimate, performance$ref_se, performance$ref_lower, performance$ref_upper, performance$ref_n_domain, performance$ref_n_events)
  interval <- newcombe_interval(performance$estimate, arm$ci_lower, arm$ci_upper, performance$ref_estimate, ref$ci_lower, ref$ci_upper)

  performance$difference_pp <- 100 * (performance$estimate - performance$ref_estimate)
  performance$ci_width_pp <- 100 * (performance$ci_upper - performance$ci_lower)
  performance$boundary <- arm$flag | ref$flag
  performance$not_different <- interval$lower <= 0 & interval$upper >= 0
  performance |>
    mutate(method = factor(method, levels = c("CBE", "RDS", "VBS")), antigen = factor(antigen, levels = primary_outcomes)) |>
    arrange(country, method, antigen)
}

## Phase 1 summary ----
# Per country and pooled over the two countries: RMSE, the noise floor (the
# RMSE expected from the sampling error of both surveys alone), mean interval
# width, the antigens with an interval at most 10 pp wide, and the antigens not
# distinguishable from the household survey, for the three primary antigens
# and for zero-dose separately.
phase1_summary <- function(performance) {
  scored <- performance |> mutate(scope = if_else(antigen == "Zero-dose", "zero_dose", "primary_antigens"))
  summarise_cells <- \(data) {
    data |>
      summarise(
        n_cells = n(),
        rmse_pp = sqrt(mean(difference_pp^2)),
        noise_floor_pp = sqrt(mean((100 * standard_error)^2 + (100 * ref_se)^2)),
        mean_ci_width_pp = mean(ci_width_pp),
        n_ci_width_le_10 = sum(ci_width_pp <= 10 + 1e-9),
        n_not_different = sum(not_different),
        n_not_different_boundary = sum(not_different & boundary),
        .groups = "drop"
      )
  }
  bind_rows(
    scored |> group_by(country, method, scope) |> summarise_cells(),
    scored |> group_by(method, scope) |> summarise_cells() |> mutate(country = "pooled")
  ) |>
    mutate(country = factor(country, levels = c("eth", "som", "pooled"))) |>
    arrange(country, method, scope)
}

## Phase 2 pairwise comparisons ----
# Every pair of methods within a country, on the four outcomes: the difference
# (first method minus second) with Newcombe's interval, count-based for a
# boundary arm.
pairwise_tests <- function(estimates, counts) {
  phase2 <-
    estimates |>
    filter(phase == 2, antigen %in% primary_outcomes) |>
    mutate(antigen = as.character(antigen), method = as.character(method)) |>
    left_join(counts, by = c("sample", "antigen"))
  pairs <- list(c("CBE", "RDS"), c("CBE", "VBS"), c("RDS", "VBS"))
  map(pairs, \(pair) {
    a <- phase2 |> filter(method == pair[[1]])
    b <- phase2 |> filter(method == pair[[2]])
    joined <- inner_join(a, b, by = c("country", "antigen"), suffix = c("_a", "_b"))
    arm_a <- boundary_arm(joined$estimate_a, joined$standard_error_a, joined$ci_lower_a, joined$ci_upper_a, joined$n_domain_a, joined$n_events_a)
    arm_b <- boundary_arm(joined$estimate_b, joined$standard_error_b, joined$ci_lower_b, joined$ci_upper_b, joined$n_domain_b, joined$n_events_b)
    interval <- newcombe_interval(joined$estimate_a, arm_a$ci_lower, arm_a$ci_upper, joined$estimate_b, arm_b$ci_lower, arm_b$ci_upper)
    tibble(
      country = joined$country, antigen = joined$antigen, method_a = pair[[1]], method_b = pair[[2]],
      diff_pp = 100 * (joined$estimate_a - joined$estimate_b),
      diff_lower_pp = 100 * interval$lower, diff_upper_pp = 100 * interval$upper,
      distinguishable = !(interval$lower <= 0 & interval$upper >= 0),
      boundary_a = arm_a$flag, boundary_b = arm_b$flag, boundary = arm_a$flag | arm_b$flag
    )
  }) |>
    list_rbind() |>
    mutate(antigen = factor(antigen, levels = primary_outcomes)) |>
    arrange(country, antigen, method_a, method_b)
}

## Inter-method range ----
# Highest minus lowest of the Phase 2 method estimates, per country and outcome.
intermethod_range <- function(estimates) {
  estimates |>
    filter(phase == 2, antigen %in% primary_outcomes) |>
    group_by(country, antigen) |>
    summarise(low = 100 * min(estimate), high = 100 * max(estimate), range_pp = high - low, .groups = "drop")
}

# Stopping point ----

## Percentile over enrolment orders ----
# A percentile of the stopping points of the random enrolment orders. An order
# that never meets the criterion counts at the full sample, and a percentile
# that falls among those orders is missing.
order_percentile <- function(x, prob, n_final) {
  rank <- max(1, ceiling(prob * length(x) - 1e-9))
  if (rank > sum(!is.na(x))) return(NA_real_)
  sort(replace_na(x, n_final))[rank]
}

# Field effort ----

## Team hours ----
# Total team hours including travel per sample, split into enumerator time on
# site, enumerator travel and supervisor time (time-motion segments plus the
# supervisors' tracker hours). CBE sheets record no supervisor time.
team_hours <- function(hours) {
  hours |>
    filter(source != "coordinator timesheet") |>
    mutate(component = case_when(
      role == "supervisor" ~ "supervisor",
      activity == "travel" ~ "enumerator_travel",
      .default = "enumerator_on_site"
    )) |>
    summarise(hours = sum(hours, na.rm = TRUE), .by = c(phase, country, method, component)) |>
    pivot_wider(names_from = component, values_from = hours, values_fill = 0) |>
    mutate(total_hours = enumerator_on_site + enumerator_travel + supervisor)
}

## Core interview duration ----
# Minutes from the end of consent to the start of the network module (the
# socio-demographic and direct modules), per directly reporting respondent. A
# module that took longer than 30 minutes (a paused form) or a negative time
# leaves the respondent without a duration.
core_interview_minutes <- function(respondents) {
  module_minutes <- function(from, to) {
    minutes <- as.numeric(difftime(to, from, units = "mins"))
    if_else(minutes < 0 | minutes > 30, NA_real_, minutes)
  }
  module_minutes(respondents$t_end_consent, respondents$t_end_socdem) +
    module_minutes(respondents$t_end_socdem, respondents$t_start_nsum)
}

## Field staff ----
# Enumerators and their person-days in the time-motion sheets, and the
# supervisors in the sheets or, in Somalia RDS and VBS, those who logged field
# hours in the weekly tracker. A person working in two methods is counted in
# each.
staff_counts <- function(hours) {
  field <- filter(hours, source != "coordinator timesheet")
  enumerators <-
    field |>
    filter(role == "enumerator") |>
    summarise(n_enumerators = n_distinct(person), n_person_days = n_distinct(person, date), .by = c(phase, country, method))
  supervisors <-
    field |>
    filter(role == "supervisor", source == "time-motion sheet" | hours > 0) |>
    summarise(n_supervisors = n_distinct(person), .by = c(phase, country, method))
  enumerators |>
    left_join(supervisors, by = c("phase", "country", "method")) |>
    mutate(n_supervisors = coalesce(n_supervisors, 0L))
}

## Contacts ----
# The people each method approached: households approached for HH and CBE
# (absent households included where a tally recorded them), persons screened
# for RDS and travellers screened for VBS. HH counts the records of the
# KoboToolbox export before cleaning; CBE forms hold the enrolled households
# only, the others coming from the tally sheets. The Somalia Phase 1 CBE sheets
# of 8 of the 30 clusters were not received, and their counts are imputed from
# the mean of the other clusters.
sample_contacts <- function(frames, tallies, record_counts, samples = study_samples) {
  map(samples, \(sample) {
    parts <- sample_parts(sample)
    records <- switch(
      parts$method,
      hh = filter(record_counts, country == parts$country, phase == parts$phase, method == "hh")$n_records,
      cbe = nrow(frames[[sample]]),
      n_distinct(frames[[sample]]$record_ix)
    )
    sheets <- filter(tallies, country == parts$country, phase == parts$phase, method == parts$method)
    if (nrow(sheets) == 0) return(tibble(sample, n_contacts = records))
    n_missing <- if (parts$method == "cbe") n_distinct(frames[[sample]]$cluster) - n_distinct(sheets$cluster) else 0
    tally <-
      sheets |>
      summarise(across(c(absent, ineligible, refused), \(x) sum(x, na.rm = TRUE))) |>
      mutate(across(everything(), \(x) x + round(x / n_distinct(sheets$cluster) * n_missing)))
    tibble(sample, n_contacts = records + tally$absent + tally$ineligible + tally$refused)
  }) |>
    list_rbind()
}

# Tables ----

## Number formats ----
# Half-up rounding as printed in the report. A negative value carries a minus
# sign (U+2212), a percentage with its 95% CI reads "78.4 (73.8–82.3)", and a
# missing value is an en dash.
round_half_up <- function(x, digits = 1) sign(x) * floor(abs(x) * 10^digits + 0.5 + 1e-9) / 10^digits

format_number <- function(x, digits = 1) {
  if_else(is.na(x), "–", formatC(round_half_up(x, digits), format = "f", digits = digits, big.mark = ","))
}

format_count <- function(x) if_else(is.na(x), "–", formatC(round_half_up(x, 0), format = "d", big.mark = ","))

format_dollars <- function(x) if_else(is.na(x), "–", paste0("$", format_count(x)))

format_signed <- function(x, digits = 1) {
  r <- round_half_up(x, digits)
  text <- sprintf(paste0("%.", digits, "f"), abs(r))
  if_else(is.na(x), "–", if_else(r < 0, paste0("−", text), text))
}

# a percentage with its 95% CI, from proportions
format_ci <- function(estimate, lower, upper) {
  sprintf("%s (%s–%s)", format_number(100 * estimate), format_number(100 * lower), format_number(100 * upper))
}

# a difference in pp with its interval
format_difference <- function(estimate, lower, upper) {
  sprintf("%s (%s to %s)", format_signed(estimate), format_signed(lower), format_signed(upper))
}

# a percentage with its 95% CI, the bounds joined by a hyphen
format_ci_hyphen <- function(estimate, lower, upper) {
  sprintf("%s (%s-%s)", format_number(100 * estimate), format_number(100 * lower), format_number(100 * upper))
}

# a difference in pp with its interval, the bounds joined by an en dash when
# both are positive and by "to" otherwise
format_difference_range <- function(estimate, lower, upper) {
  positive <- round_half_up(lower) >= 0 & round_half_up(upper) >= 0
  if_else(
    positive,
    sprintf("%s (%s–%s)", format_signed(estimate), format_signed(lower), format_signed(upper)),
    sprintf("%s (%s to %s)", format_signed(estimate), format_signed(lower), format_signed(upper))
  )
}

## Labels ----
# Outcome and verdict wording, and a "phase, country" label, as printed in the
# annex tables.
outcome_label <- function(antigen) if_else(antigen == "Zero-dose", "Zero-dose prevalence", as.character(antigen))

verdict_label <- function(distinguishable) if_else(distinguishable, "Distinguishable", "Not distinguishable")

phase_country <- function(phase, country) str_c(phase, ", ", country_labels[country])

## Save ----
# One UTF-8 CSV per report table.
save_table <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_csv(table, path)
}

# Figures ----

## Style ----
# Method colours (HH is the charcoal reference), text colours and labels used
# in every figure.
method_order <- c("hh", "cbe", "rds", "vbs")
method_labels <- c(hh = "HH", cbe = "CBE", rds = "RDS", vbs = "VBS")
method_colours <- c(hh = "#4D4D4D", cbe = "#1B9E77", rds = "#7570B3", vbs = "#D95F02")
accent_colour <- "#B1005B"
ink <- c(primary = "#1E1E22", secondary = "#4A4A52", muted = "#7A7A82", grid = "#E6E6EA")
country_labels <- c(eth = "Ethiopia", som = "Somalia")
phase_labels <- c(`1` = "Phase 1", `2` = "Phase 2")
figure_width_cm <- 16

## Font ----
# Aptos, the report font, ships with Microsoft Office; figures fall back to
# Arial where it is not installed.
register_report_font <- function() {
  faces <- c("Aptos.ttf", "Aptos-Bold.ttf", "Aptos-Italic.ttf", "Aptos-Bold-Italic.ttf")
  folders <- Sys.glob("/Applications/Microsoft */Contents/Resources/DFonts")
  folders <- keep(folders, \(folder) all(file.exists(file.path(folder, faces))))
  if (length(folders) == 0) return("Arial")
  paths <- file.path(folders[[1]], faces)
  systemfonts::register_font("Aptos", plain = paths[1], bold = paths[2], italic = paths[3], bolditalic = paths[4])
  "Aptos"
}
report_font <- register_report_font()

## Theme ----
# White surface, hairline grid, no ticks, legend at the bottom.
theme_report <- function(base_size = 8.5) {
  theme_minimal(base_size = base_size, base_family = report_font) +
    theme(
      text = element_text(colour = ink[["primary"]]),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = ink[["grid"]], linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      strip.text = element_text(face = "bold", size = base_size + 0.5, colour = ink[["primary"]], margin = margin(2, 3, 4, 3)),
      strip.placement = "outside",
      axis.text = element_text(colour = ink[["secondary"]], size = base_size - 1),
      axis.title = element_text(colour = ink[["primary"]], size = base_size),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 0.5, colour = ink[["primary"]]),
      legend.key.size = unit(0.4, "cm"),
      legend.margin = margin(t = 0),
      plot.margin = margin(6, 8, 4, 6)
    )
}

## Method scale ----
method_scale <- function(methods = method_order, aesthetic = c("colour", "fill")) {
  aesthetic <- match.arg(aesthetic)
  methods <- method_order[method_order %in% methods]
  scale_discrete_manual(aesthetic, values = method_colours[methods], breaks = methods, labels = method_labels[methods], name = "Method")
}

## Save ----
# PNG at 300 dpi, 16 cm wide (the report text column).
save_figure <- function(plot, path, height_cm, width_cm = figure_width_cm) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, device = ragg::agg_png, width = width_cm, height = height_cm, units = "cm", dpi = 300, bg = "white")
}

## Country divider ----
# A hairline across the full width between two equal-height country rows.
country_divider <- function(plot, height_cm) {
  # lay the plot out on a device of the size it is saved at
  ragg::agg_png(tempfile(fileext = ".png"), width = figure_width_cm, height = height_cm, units = "cm", res = 300)
  on.exit(dev.off())
  gt <- if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot) else ggplotGrob(plot)
  area <- which(gt$layout$name == "panel-area")
  if (length(area) == 1) {
    top <- gt$layout$t[[area]]
    bottom <- gt$layout$b[[area]]
  } else {
    panels <- gt$layout[str_detect(gt$layout$name, "^panel"), ]
    top <- min(panels$t)
    bottom <- max(panels$b)
  }
  rule <- grid::segmentsGrob(
    x0 = unit(0.15, "cm"), x1 = unit(1, "npc") - unit(0.15, "cm"), y0 = unit(0.5, "npc"), y1 = unit(0.5, "npc"),
    gp = grid::gpar(col = ink[["muted"]], lwd = 0.3 * .pt)
  )
  gt <- gtable::gtable_add_grob(gt, rule, t = top, l = 1, b = bottom, r = ncol(gt), z = Inf, clip = "off", name = "country-divider")
  patchwork::wrap_elements(full = gt)
}

## Coverage forest ----
# Methods CBE, RDS and VBS on x, coverage on y, countries as rows and
# outcomes as columns. The three antigens share a 0-100% axis and zero-dose
# has its own. In Phase 1 the household survey is drawn as a line at its
# estimate over a band at its 95% CI; in Phase 2 a bracket spans the lowest
# to the highest estimate and is labelled with the inter-method range.
forest_block <- function(alt, hh, range_pp, ylim, ybreaks, first) {
  methods <- c("cbe", "rds", "vbs")
  p <- ggplot()
  if (!is.null(hh)) {
    p <- p +
      geom_rect(data = hh, aes(xmin = -Inf, xmax = Inf, ymin = 100 * ci_lower, ymax = 100 * ci_upper), fill = method_colours[["hh"]], alpha = 0.13) +
      geom_hline(data = hh, aes(yintercept = 100 * estimate, linetype = "hh"), colour = method_colours[["hh"]], linewidth = 0.55)
  }
  if (!is.null(range_pp)) {
    bx <- 3.42
    p <- p +
      geom_segment(data = range_pp, aes(x = bx, xend = bx, y = low, yend = high), colour = ink[["muted"]], linewidth = 0.4) +
      geom_segment(data = range_pp, aes(x = bx - 0.07, xend = bx + 0.07, y = low, yend = low), colour = ink[["muted"]], linewidth = 0.4) +
      geom_segment(data = range_pp, aes(x = bx - 0.07, xend = bx + 0.07, y = high, yend = high), colour = ink[["muted"]], linewidth = 0.4) +
      geom_text(
        data = range_pp, aes(x = bx + 0.14, y = (low + high) / 2, label = label), hjust = 0, size = 2.4,
        colour = ink[["secondary"]], family = report_font
      )
  }
  p <- p +
    geom_errorbar(data = alt, aes(x = method_x, ymin = 100 * ci_lower, ymax = 100 * ci_upper, colour = method_c), width = 0.18, linewidth = 0.55) +
    geom_point(data = alt, aes(x = method_x, y = 100 * estimate, fill = method_c), shape = 21, colour = "white", size = 2.6, stroke = 0.6) +
    facet_grid(country_f ~ outcome_f, switch = "y") +
    method_scale(methods, "fill") +
    scale_colour_manual(values = method_colours[methods], guide = "none") +
    (if (is.null(hh)) NULL else scale_linetype_manual(values = c(hh = "solid"), name = NULL, labels = c(hh = "HH estimate and 95% CI"))) +
    scale_x_discrete(expand = expansion(add = c(0.6, if (is.null(range_pp)) 0.6 else 2.05))) +
    scale_y_continuous(limits = ylim, breaks = ybreaks, labels = \(x) paste0(x, "%"), expand = c(0.02, 0)) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 2.6, colour = "white")),
      linetype = if (is.null(hh)) "none" else guide_legend(order = 2, override.aes = list(linewidth = 0.55))
    ) +
    labs(x = NULL, y = NULL) +
    theme_report() +
    theme(panel.grid.major.x = element_blank(), panel.spacing.y = unit(2, "lines"))
  if (!first) p <- p + theme(strip.text.y.left = element_blank(), strip.background.y = element_blank())
  p
}

# `estimates` has country, phase, method (HH, CBE, RDS, VBS), antigen,
# estimate, ci_lower and ci_upper; `ranges` (Phase 2) has country, antigen,
# low, high and range_pp.
coverage_forest <- function(estimates, phase, ranges = NULL, height_cm = 11) {
  labelled <- \(d) mutate(
    d,
    country_f = factor(country_labels[country], levels = country_labels),
    outcome_f = factor(as.character(antigen), levels = primary_outcomes)
  )
  d <-
    estimates |>
    filter(phase == .env$phase, antigen %in% primary_outcomes) |>
    mutate(method = str_to_lower(as.character(method))) |>
    labelled()
  alt <-
    d |>
    filter(method != "hh") |>
    mutate(
      method_x = factor(method_labels[method], levels = method_labels[c("cbe", "rds", "vbs")]),
      method_c = factor(method, levels = c("cbe", "rds", "vbs"))
    )
  hh <- if (phase == 1) filter(d, method == "hh") else NULL
  range_pp <- if (phase == 2) ranges |> labelled() |> mutate(label = sprintf("%.1f pp", range_pp)) else NULL
  is_antigen <- alt$antigen != "Zero-dose"
  zero_top <- ceiling(max(100 * alt$ci_upper[!is_antigen]) / 10) * 10
  part <- function(antigens) {
    keep <- \(x) if (is.null(x)) NULL else filter(x, (antigen != "Zero-dose") == antigens)
    forest_block(
      keep(alt), keep(hh), keep(range_pp),
      ylim = if (antigens) c(0, 100) else c(-2, zero_top),
      ybreaks = if (antigens) seq(0, 100, 20) else seq(0, zero_top, 10),
      first = antigens
    )
  }
  p <- patchwork::wrap_plots(part(TRUE), part(FALSE), nrow = 1, widths = c(3, 1), guides = "collect") &
    theme(legend.position = "bottom")
  country_divider(p, height_cm)
}

## Scorecards ----
# Tile grids of the verdict of each method on each measure and outcome. Each
# tile is filled by its verdict and prints its value; a key below the panels
# explains the fills and the tile text.
verdict_levels <- c("meets", "fails", "not_readable", "no_threshold", "beyond_ns", "not_read")
verdict_fill <- c(meets = "#CBE8D6", fails = "#F5CDBB", not_readable = "white", no_threshold = "#EEEEF0", beyond_ns = "#FBE7B5", not_read = "#EEEEF0")
verdict_edge <- c(meets = "white", fails = "white", not_readable = ink[["muted"]], no_threshold = "white", beyond_ns = "white", not_read = "white")
verdict_linetype <- c(meets = "solid", fails = "solid", not_readable = "dashed", no_threshold = "solid", beyond_ns = "solid", not_read = "solid")
measure_standards <- c(
  Accuracy = "within ±10 pp (±5 pp zero-dose)", Precision = "95% CI width ≤ 10 pp",
  Calibration = "read beyond the accuracy tolerance", Agreement = "range < 10 pp; pairwise tests"
)

# width of strings in cm, measured in the figure font
text_width_cm <- function(x, size_pt) {
  ragg::agg_png(tempfile(fileext = ".png"))
  on.exit(grDevices::dev.off(), add = TRUE)
  systemfonts::string_widths_dev(x, family = report_font, size = size_pt, unit = "cm")
}

# one key item: a swatch in a verdict's style, optionally carrying tile text,
# then its explanation; a text-only note when verdict is NULL
key_item <- function(text, verdict = NULL, tile = NULL) list(text = text, verdict = verdict, tile = tile)

# lay the key rows out centred across the figure width, drawn in centimetres
scorecard_key <- function(rows, text_pt = 7, tile_pt = 6.5, swatch_h = 0.36, row_h = 0.44, item_gap = 0.45, text_gap = 0.14, swatch_min_w = 0.5) {
  width <- figure_width_cm - 14 / 72.27 * 2.54
  swatches <- list()
  texts <- list()
  for (i in seq_along(rows)) {
    items <- rows[[i]]
    text_w <- text_width_cm(map_chr(items, \(it) it$text), text_pt)
    tile_w <- map_dbl(items, \(it) if (is.null(it$tile)) 0 else text_width_cm(it$tile, tile_pt))
    swatch_w <- map_dbl(seq_along(items), \(k) if (is.null(items[[k]]$verdict)) 0 else max(swatch_min_w, tile_w[k] + 0.24))
    total <- sum(if_else(swatch_w > 0, swatch_w + text_gap, 0) + text_w) + item_gap * (length(items) - 1)
    x <- (width - total) / 2
    y <- (length(rows) - i + 0.5) * row_h
    for (k in seq_along(items)) {
      it <- items[[k]]
      if (!is.null(it$verdict)) {
        v <- it$verdict
        swatches[[length(swatches) + 1]] <- tibble(
          xmin = x, xmax = x + swatch_w[k], ymin = y - swatch_h / 2, ymax = y + swatch_h / 2,
          fill = verdict_fill[[v]], edge = verdict_edge[[v]], linetype = verdict_linetype[[v]]
        )
        if (!is.null(it$tile)) {
          texts[[length(texts) + 1]] <- tibble(
            x = x + swatch_w[k] / 2, y = y, label = it$tile, hjust = 0.5, size = tile_pt,
            face = if (v == "not_readable") "italic" else "plain", colour = if (v == "not_readable") ink[["muted"]] else ink[["primary"]]
          )
        }
        x <- x + swatch_w[k] + text_gap
      }
      texts[[length(texts) + 1]] <- tibble(x = x, y = y, label = it$text, hjust = 0, size = text_pt, face = "plain", colour = ink[["primary"]])
      x <- x + text_w[k] + item_gap
    }
  }
  swatches <- bind_rows(swatches)
  texts <- bind_rows(texts)
  plot <- ggplot() +
    geom_rect(data = swatches, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill), colour = I(edge), linetype = I(linetype)), linewidth = 0.4) +
    geom_text(
      data = texts, aes(x = x, y = y, label = label, hjust = hjust, size = I(size / .pt), fontface = face, colour = I(colour)),
      family = report_font
    ) +
    scale_x_continuous(limits = c(0, width), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, length(rows) * row_h), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(0, 8, 0, 6))
  list(plot = plot, height_cm = length(rows) * row_h)
}

# draw the tiles, faceted by country (rows) and measure (columns), with the key below
scorecard_plot <- function(tiles, verdicts, key) {
  strip_label <- \(measure) {
    sprintf("<b>%s</b><br><span style='color:%s;font-size:7pt'>%s</span>", measure, ink[["muted"]], measure_standards[as.character(measure)])
  }
  readable <- tiles |> filter(verdict != "not_readable")
  count_based <- tiles |> filter(verdict == "not_readable")
  main <- ggplot(tiles, aes(outcome, row)) +
    geom_tile(aes(fill = verdict, colour = verdict, linetype = verdict), width = 0.94, height = 0.9, linewidth = 0.4) +
    geom_text(
      data = readable, aes(label = label, fontface = if_else(verdict == "fails", "bold", "plain")),
      size = 2.6, colour = ink[["primary"]], family = report_font, lineheight = 0.9
    ) +
    geom_text(data = count_based, aes(label = label), size = 2.3, colour = ink[["muted"]], family = report_font, fontface = "italic", lineheight = 0.9) +
    facet_grid(country ~ measure, scales = "free_y", space = "free_y", switch = "y", labeller = labeller(measure = strip_label, country = label_value)) +
    scale_fill_manual(values = verdict_fill, limits = verdicts, drop = FALSE) +
    scale_colour_manual(values = verdict_edge, limits = verdicts, drop = FALSE) +
    scale_linetype_manual(values = verdict_linetype, limits = verdicts, drop = FALSE) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(x = NULL, y = NULL) +
    theme_report() +
    theme(
      strip.text.x = ggtext::element_markdown(family = report_font, size = 9, colour = ink[["primary"]], face = "plain", margin = margin(2, 3, 4, 3)),
      panel.grid.major = element_blank(), panel.spacing = unit(0.35, "lines"),
      axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 8, lineheight = 0.9),
      legend.position = "none"
    )
  patchwork::wrap_plots(main, patchwork::wrap_elements(full = key$plot), ncol = 1, heights = grid::unit.c(unit(1, "null"), unit(key$height_cm, "cm"))) +
    patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
}
