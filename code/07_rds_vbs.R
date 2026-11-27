# 07_rds_vbs.R
# Method diagnostics: VBS travel frequency, mobility weights and venues; design
# effects; RDS recruitment, homophily and convergence
# Reproduces: Annex Tables E1-E5, Annex Figures E1-E4

# Setup ----

pacman::p_load(survey, RDS, igraph, ggh4x, patchwork, tidyverse)

source("code/functions.R")

# Data ----

frames <- read_rds("data/clean/frames.rds")
estimates <- read_rds("data/clean/estimates.rds")
tracking <- read_rds("data/clean/rds_tracking.rds")
stopping_rule <- read_rds("data/clean/stopping.rds")

# the samples of each method, in report order
samples <- study_samples
vbs_samples <- str_subset(samples, "vbs")
rds_samples <- str_subset(samples, "rds")

# each sample's directly reported children
cohorts <- imap(frames[samples], direct_cohort)

# the four outcomes of the diagnostics and their fields
outcome_fields <- antigen_outcomes[primary_outcomes]

# "Ethiopia, 1" row labels
country_phase <- function(country, phase) str_c(country_labels[country], ", ", phase)

# Analysis ----

## VBS travel frequency and mobility weights ----
# Each directly reporting caregiver's answer on how often they travelled to the
# town of the venue where they were met, with the score and relative weight the
# primary VBS estimator gives that answer (Galkayo North answers were per week
# and are put on a 4-week basis), and the Kish design effect of the weights
# over child records.
# town and the period each town's question asked about
vbs_town <- function(cohort, sample) {
  switch(
    sample,
    eth_phase1_vbs = if_else(str_detect(cohort$int_location, "Worabe"), "Worabe", "Inseno"),
    som_phase1_vbs = if_else(str_detect(cohort$int_location, "G North"), "Galkayo North", "Adado"),
    eth_phase2_vbs = rep("Bahir Dar", nrow(cohort)),
    som_phase2_vbs = rep("Kismayo", nrow(cohort))
  )
}
towns <- c("Worabe", "Inseno", "Bahir Dar", "Galkayo North", "Adado", "Kismayo")
period_asked <- c(Worabe = "past 4 weeks", Inseno = "past 4 weeks", `Bahir Dar` = "past 4 weeks", `Galkayo North` = "per week", Adado = "past 4 weeks", Kismayo = "per week")
travel_levels <- c("0", "1-4", "5-9", "10+", "Prefer not to answer")

# each child record with its caregiver's town, answer, score and weight
travel <- map(vbs_samples, \(sample) {
  cohort <- cohorts[[sample]]
  answer <- vbs_travel(cohort, sample)
  score <- unname(c(`0` = 2.5, `1-4` = 2.5, `5-9` = 7, `10+` = 10)[answer] * travel_multiplier(cohort, sample))
  tibble(
    record_ix = cohort$record_ix,
    town = vbs_town(cohort, sample),
    answer = factor(answer, levels = travel_levels),
    score = score,
    weight = vbs_weight(cohort, sample)
  ) |>
    mutate(sample = .env$sample, .before = 1)
}) |>
  list_rbind() |>
  left_join(sample_parts(vbs_samples), by = "sample") |>
  mutate(town = factor(town, levels = towns))

# caregivers by town and answer
travel_caregivers <-
  travel |>
  distinct(sample, country, phase, town, record_ix, answer) |>
  count(country, phase, town, answer, .drop = FALSE, name = "n") |>
  mutate(pct = 100 * n / sum(n), .by = c(country, phase, town)) |>
  filter(sum(n) > 0, .by = c(country, phase, town))

# score and weight of each answer at each town
travel_weights <-
  travel |>
  summarise(score = first(score), weight = first(weight), .by = c(country, phase, town, answer))

# Kish design effect of the weights in each sample
kish <- travel |> summarise(kish_deff = n() * sum(weight^2) / sum(weight)^2, .by = c(country, phase))

## VBS coverage by venue ----
# Each outcome within each venue, on the sample's weighted design. A venue with
# fewer than 10 children observed has no estimate; an estimate of 0 or 100%
# takes the count-based Wilson interval.
# estimate each outcome in each venue of each sample
venue_estimates <- map(vbs_samples, \(sample) {
  cohort <- cohorts[[sample]]
  cohort$weight <- vbs_weight(cohort, sample)
  design <- svydesign(id = ~record_ix, weights = ~weight, data = cohort)
  map(str_sort(unique(cohort$int_location), locale = "en"), \(venue) {
    in_venue <- cohort$int_location == venue
    map(names(outcome_fields), \(antigen) {
      field <- outcome_fields[[antigen]]
      n_observed <- sum(!is.na(cohort[[field]][in_venue]))
      result <- tibble(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
      count_based <- FALSE
      if (n_observed >= 10) {
        result <- survey_proportion(subset(design, int_location == venue), field)
        if (abs(result$estimate) < 1e-9 || abs(result$estimate - 1) < 1e-9) {
          interval <- count_wilson(if (result$estimate < 0.5) 0 else n_observed, n_observed)
          result$ci_lower <- interval[[1]]
          result$ci_upper <- interval[[2]]
          count_based <- TRUE
        }
      }
      tibble(
        sample, venue, antigen, estimate = result$estimate, ci_lower = result$ci_lower, ci_upper = result$ci_upper,
        n_children = sum(in_venue), small_n = n_observed < 10, count_based
      )
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  left_join(sample_parts(vbs_samples), by = "sample")

## Design effects ----
# The variance of each primary estimate divided by that of a simple random
# sample of the children observed, and, for the cluster designs, the
# intraclass correlation from the design effect of the cluster totals and the
# mean children per cluster.
# design effect of every sample and outcome
counts <- event_counts(frames)
design_effects <-
  estimates |>
  filter(antigen %in% primary_outcomes) |>
  mutate(antigen = as.character(antigen), method = str_to_lower(as.character(method))) |>
  left_join(counts, by = c("sample", "antigen")) |>
  mutate(
    deff = standard_error^2 / (estimate * (1 - estimate) / n_observed),
    sparse = n_events < 5 | (n_domain - n_events) < 5
  )

# intraclass correlation of one outcome in a cluster sample
cluster_icc <- function(y, cluster) {
  keep <- !is.na(y)
  y <- y[keep]
  cluster <- cluster[keep]
  p <- mean(y)
  residual <- tapply(y, cluster, sum) - p * tapply(y, cluster, length)
  k <- length(residual)
  deff <- (k / (k - 1)) * sum(residual^2) / length(y)^2 / (p * (1 - p) / length(y))
  (deff - 1) / (length(y) / k - 1)
}

# ICC of every HH and CBE sample and outcome
iccs <- map(str_subset(samples, "hh|cbe"), \(sample) {
  cohort <- cohorts[[sample]]
  cluster <- if (str_detect(sample, "cbe")) cohort$cluster else cohort$cluster_id
  tibble(sample, antigen = names(outcome_fields), icc = map_dbl(outcome_fields, \(field) cluster_icc(cohort[[field]], cluster)))
}) |>
  list_rbind()

## RDS recruitment networks ----
# The recruitment network of each RDS sample as estimated: one node per
# respondent, plus in Somalia Phase 2 the recruiters named on the tracking
# sheet who were not interviewed. Each node's wave is its distance from the
# start of its chain.
# build each network, with waves and chains
networks <- map(set_names(rds_samples), \(sample) {
  cohort <- cohorts[[sample]]
  nodes <- rds_network(cohort, outcome_fields, tracking = if (sample == "som_phase2_rds") tracking) |> as.data.frame()
  design <- as.rds.data.frame(nodes, id = "id_participant_calc", recruiter.id = "id_referring", network.size = "degree", population.size = population_size[[sample]])
  nodes$wave <- as.integer(get.wave(design))
  nodes$interviewed <- nodes$id_participant_calc %in% as.character(cohort$id_participant_calc)
  nodes$chain <- chain_root(nodes)
  chain_size <- table(nodes$chain)
  chain_order <- names(chain_size)[order(-as.integer(chain_size), as.numeric(names(chain_size)))]
  nodes$chain_rank <- match(nodes$chain, chain_order)
  nodes$chain_size <- as.integer(chain_size[nodes$chain])
  list(cohort = cohort, nodes = nodes, design = design)
})

## Recruitment structure ----
# Seeds (and those who recruited at least one respondent), respondents, chains
# (and those of 10 or more nodes), waves of the respondents, the largest
# chain's share of nodes, and the mean recruits per recruiter. Recruiters named
# but not interviewed count in the chains and recruits, not in respondents or
# waves.
# summarise each network
recruitment <- imap(networks, \(network, sample) {
  nodes <- network$nodes
  respondents <- distinct(network$cohort, id_participant_calc, .keep_all = TRUE)
  seeds <- as.character(respondents$id_participant_calc[respondents$seed_participant %in% "Yes"])
  waves <- nodes$wave[nodes$interviewed]
  chain_size <- table(nodes$chain)
  tibble(
    sample,
    n_seeds = length(seeds),
    n_seeds_productive = sum(seeds %in% nodes$id_referring[nodes$interviewed]),
    n_respondents = nrow(respondents),
    n_chains = length(chain_size),
    n_chains_10 = sum(chain_size >= 10),
    max_wave = max(waves),
    median_wave = median(waves),
    p25_wave = quantile(waves, 0.25, names = FALSE),
    p75_wave = quantile(waves, 0.75, names = FALSE),
    largest_chain_pct = 100 * max(chain_size) / nrow(nodes),
    recruits_per_recruiter = mean(table(nodes$id_referring[nodes$id_referring != "0"]))
  )
}) |>
  list_rbind() |>
  left_join(sample_parts(rds_samples), by = "sample")

## Recruitment homophily ----
# For each variable, the recruiter and recruit pairs in the same category
# divided by the number expected if respondents recruited at random with
# respect to the variable (1: no homophily), over the pairs in which both have
# a value. A respondent's outcome is their children's, left out where their
# children differ. Education is in four stages; wealth is the asset group.
# homophily of one variable, with the share of same-category pairs and
# whether any expected cell count of the pair table is below 5
homophily <- function(design, values, population) {
  design$hom_variable <- factor(values)
  fit <- homophily.estimates(design, outcome.variable = "hom_variable", recruitment = TRUE, N = population)
  pairs <- as.matrix(unclass(fit@estimate))
  expected <- suppressWarnings(chisq.test(pairs))$expected
  tibble(homophily = as.numeric(attr(fit@estimate, "homophily")), same_pct = 100 * sum(diag(pairs)) / sum(pairs), sparse = min(expected) < 5)
}

# variables in display order
homophily_variables <- c(
  Penta3 = "penta_3p_clean", OPV3 = "opv_3p_clean", MCV1 = "measles_1p_clean", `Zero-dose status` = "zero_dose",
  Residence = "residence", Education = "education", `Household wealth` = "wealth"
)

# homophily of each variable in each sample
homophily_results <- imap(networks, \(network, sample) {
  nodes <- network$nodes
  respondents <- distinct(network$cohort, id_participant_calc, .keep_all = TRUE)
  at <- match(nodes$id_participant_calc, as.character(respondents$id_participant_calc))
  values <- c(
    map(homophily_variables[1:4], \(field) if_else(nodes[[field]] %in% c(0, 1), as.character(nodes[[field]]), NA_character_)),
    list(
      residence = as.character(respondents$residence)[at],
      education = unname(education_stages[as.character(respondents$education)])[at],
      wealth = wealth_group(respondents, common_assets(cohorts, sample))[at]
    )
  )
  names(values) <- names(homophily_variables)
  imap(values, \(v, variable) homophily(network$design, v, population_size[[sample]]) |> mutate(sample, variable)) |>
    list_rbind()
}) |>
  list_rbind() |>
  left_join(sample_parts(rds_samples), by = "sample")

## Cumulative estimates within chains ----
# Penta3 among respondents in interview order (date, then wave), within each
# chain and over the whole sample, each respondent weighted by the inverse of
# their network size. A descriptive trajectory, not the RDS estimator.
# cumulative inverse-degree weighted Penta3 per chain and overall
chain_trajectories <- imap(networks, \(network, sample) {
  d <-
    network$nodes |>
    filter(interviewed, !is.na(penta_3p_clean)) |>
    arrange(date, wave, as.numeric(id_participant_calc)) |>
    mutate(weight = 1 / pmax(degree, 1))
  cumulative <- \(x, chain) tibble(chain, n = seq_len(nrow(x)), estimate = cumsum(x$weight * x$penta_3p_clean) / cumsum(x$weight))
  bind_rows(
    map(sort(unique(d$chain_rank)), \(rank) cumulative(filter(d, chain_rank == rank), sprintf("chain_%02d", rank))) |> list_rbind(),
    cumulative(d, "all")
  ) |>
    mutate(sample)
}) |>
  list_rbind() |>
  left_join(sample_parts(rds_samples), by = "sample")

# Tables ----

## Annex Table E1 ----
# Caregivers by number of trips to the venue town, by site.
# format count and percentage, one row per site
table_E01 <-
  travel_caregivers |>
  arrange(country, phase, town) |>
  mutate(cell = str_c(format_count(n), " (", format_number(pct), "%)"), caregivers = sum(n), .by = c(country, phase, town)) |>
  select(country, phase, town, answer, cell, caregivers) |>
  pivot_wider(names_from = answer, values_from = cell) |>
  transmute(
    `Country, phase` = country_phase(country, phase), Site = as.character(town), `Period asked` = period_asked[as.character(town)],
    across(all_of(travel_levels)), Caregivers = format_count(caregivers)
  )

save_table(table_E01, "output/tables/table_E01.csv")

## Annex Table E2 ----
# Score and relative weight of each answer as applied, by site (a dash where no
# caregiver gave the answer), and the Kish design effect of the weights.
# one cell per site listing the four trip categories
trip_cells <- function(values, answers, digits) {
  map_chr(travel_levels[1:4], \(level) if (level %in% answers) format_number(values[answers == level], digits) else "–") |> str_c(collapse = " / ")
}

table_E02 <-
  travel_weights |>
  arrange(country, phase, town) |>
  summarise(
    score_cell = trip_cells(score, answer, 1),
    weight_cell = trip_cells(weight, answer, 3),
    declined = if ("Prefer not to answer" %in% answer) format_number(weight[answer == "Prefer not to answer"], 3) else "–",
    .by = c(country, phase, town)
  ) |>
  left_join(kish, by = c("country", "phase")) |>
  transmute(
    `Country, phase` = country_phase(country, phase), Site = as.character(town),
    `Score: 0 / 1-4 / 5-9 / 10+ trips` = score_cell,
    `Relative weight: 0 / 1-4 / 5-9 / 10+ trips` = weight_cell,
    `Relative weight, prefer not to answer` = declined,
    `Kish design effect of the weights` = format_number(kish_deff, 3)
  )

save_table(table_E02, "output/tables/table_E02.csv")

## Annex Table E3 ----
# Design effect and ICC by sample and outcome; a dagger marks fewer than five
# children with or without the outcome.
# one row per sample, DEFF then ICC columns
table_E03 <-
  design_effects |>
  left_join(iccs, by = c("sample", "antigen")) |>
  mutate(
    deff = str_c(format_number(deff, 2), if_else(sparse, "†", "")),
    icc = format_number(icc, 3),
    antigen = if_else(antigen == "Zero-dose", "zero-dose", antigen),
    method = factor(method, levels = method_order)
  ) |>
  arrange(phase, country, method) |>
  select(phase, country, method, antigen, deff, icc) |>
  pivot_wider(names_from = antigen, values_from = c(deff, icc), names_glue = "{toupper(.value)} {antigen}") |>
  transmute(`Country, phase` = country_phase(country, phase), Method = method_labels[as.character(method)], across(starts_with(c("DEFF", "ICC"))))

save_table(table_E03, "output/tables/table_E03.csv")

## Annex Table E4 ----
# Seeds, respondents, chains, waves and recruits per recruiter, RDS samples.
# format one row per sample
table_E04 <-
  recruitment |>
  arrange(phase, country) |>
  transmute(
    `Country, phase` = country_phase(country, phase),
    `Seeds enrolled (productive)` = str_glue("{n_seeds} ({n_seeds_productive})"),
    Respondents = format_count(n_respondents),
    `Chains (10 or more respondents)` = str_glue("{n_chains} ({n_chains_10})"),
    `Maximum wave` = as.character(max_wave),
    `Median wave (IQR)` = str_glue("{median_wave} ({p25_wave}–{p75_wave})"),
    `Largest chain (% of sample)` = format_number(largest_chain_pct),
    `Recruits per recruiter` = format_number(recruits_per_recruiter)
  )

save_table(table_E04, "output/tables/table_E04.csv")

## Annex Table E5 ----
# Homophily with the percentage of same-category pairs; a dagger marks fewer
# than five pairs expected in a cell of the recruiter-by-recruit table.
# one row per variable, one column per sample
table_E05 <-
  homophily_results |>
  arrange(phase, country) |>
  mutate(
    cell = str_c(format_number(homophily, 2), " (", format_number(same_pct), ")", if_else(sparse, "†", "")),
    column = country_phase(country, phase),
    variable = factor(variable, levels = names(homophily_variables))
  ) |>
  select(variable, column, cell) |>
  pivot_wider(names_from = column, values_from = cell) |>
  arrange(variable) |>
  rename(Variable = variable) |>
  mutate(Variable = as.character(Variable))

save_table(table_E05, "output/tables/table_E05.csv")

# Figures ----

## Annex Figure E1 ----
# Each venue's estimate and 95% CI per outcome, venues grouped by country and
# phase, with the children per venue in the right-hand column; in Phase 1 the
# household survey's estimate and 95% CI.
# venue names as the annex text gives them: short Phase 1 names and the
# unnamed Somalia Phase 2 bus station spelt out, "Hospital" capitalised in a
# proper name
venue_display <- c(
  "Inseno bus" = "Inseno bus station", "Worabe bus" = "Worabe bus station",
  "G North market" = "Galkayo North market", "G North MCH" = "Galkayo North MCH",
  "Felege Hiwot Comprehensive hospital" = "Felege Hiwot Comprehensive Hospital", "Bus station" = "Kismayo bus station"
)

# venue labels without a leading "At (the)", in country and phase order
venues <-
  venue_estimates |>
  mutate(
    venue_label = str_remove(venue, "^At (the )?") |> (\(x) str_c(str_to_upper(str_sub(x, 1, 1)), str_sub(x, 2)))(),
    venue_label = coalesce(venue_display[venue_label], venue_label),
    country_f = factor(country_labels[country], levels = country_labels),
    phase_f = factor(phase_labels[as.character(phase)], levels = phase_labels),
    outcome_f = factor(antigen, levels = c(primary_outcomes, "Children (n)")),
    style = if_else(count_based, "count", "estimate")
  ) |>
  arrange(country, phase)
venues$venue_f <- factor(venues$venue_label, levels = rev(unique(venues$venue_label)))

# children per venue
venue_n <-
  venues |>
  distinct(country_f, phase_f, venue_f, n_children, small_n) |>
  summarise(n_children = first(n_children), small_n = all(small_n), .by = c(country_f, phase_f, venue_f)) |>
  mutate(outcome_f = factor("Children (n)", levels = levels(venues$outcome_f)), label = if_else(small_n, str_c(n_children, ", too few"), as.character(n_children)))

# household survey estimates, Phase 1
venue_hh <-
  estimates |>
  filter(phase == 1, method == "HH", antigen %in% primary_outcomes) |>
  transmute(
    country_f = factor(country_labels[country], levels = country_labels), phase_f = factor(phase_labels[["1"]], levels = phase_labels),
    outcome_f = factor(as.character(antigen), levels = levels(venues$outcome_f)),
    estimate = 100 * estimate, lower = 100 * ci_lower, upper = 100 * ci_upper
  )

# shown estimates, and the x scale of each column
shown <- filter(venues, !is.na(estimate))
column_scales <- c(
  rep(list(scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100), labels = c("0", "50", "100"), expand = expansion(add = 4))), 4),
  list(scale_x_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0)))
)

# create venue forest plot
figure_E1 <-
  ggplot() +
  geom_rect(data = venue_hh, aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf), fill = method_colours[["hh"]], alpha = 0.13) +
  geom_vline(data = venue_hh, aes(xintercept = estimate, linetype = "HH estimate, 95% CI"), colour = method_colours[["hh"]], linewidth = 0.55) +
  geom_linerange(data = shown, aes(y = venue_f, xmin = 100 * ci_lower, xmax = 100 * ci_upper), colour = method_colours[["vbs"]], linewidth = 0.55, orientation = "y") +
  geom_point(data = shown, aes(100 * estimate, venue_f, fill = style), shape = 21, colour = method_colours[["vbs"]], size = 2.5, stroke = 0.8) +
  geom_text(data = venue_n, aes(x = 0.02, y = venue_f, label = label), hjust = 0, size = 2.4, colour = ink[["secondary"]], family = report_font) +
  facet_nested(
    rows = vars(country_f, phase_f), cols = vars(outcome_f), scales = "free", space = "free_y", switch = "y",
    strip = strip_nested(
      text_y = elem_list_text(angle = c(90, 90, 0, 0, 0, 0), face = "bold", hjust = c(0.5, 0.5, 1, 1, 1, 1)),
      by_layer_y = FALSE, clip = "off"
    )
  ) +
  facetted_pos_scales(x = column_scales) +
  force_panelsizes(cols = unit(c(1, 1, 1, 1, 0.8), "null")) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(
    values = c(estimate = method_colours[["vbs"]], count = "white"), breaks = c("estimate", "count"),
    labels = c("VBS estimate, 95% CI", "Count-based interval"), name = NULL,
    guide = guide_legend(order = 1, override.aes = list(linetype = 0, size = 2.5, shape = 21, colour = method_colours[["vbs"]]))
  ) +
  scale_linetype_manual(
    values = c(`HH estimate, 95% CI` = "solid"), name = NULL,
    guide = guide_legend(order = 2, override.aes = list(colour = method_colours[["hh"]], linewidth = 0.55))
  ) +
  labs(x = "Coverage (%), children aged 12 to 23 months", y = NULL) +
  scale_y_discrete(labels = scales::label_wrap(24)) +
  theme_report() +
  theme(
    panel.grid.major.y = element_blank(), strip.placement = "outside",
    axis.text.y = element_text(colour = ink[["primary"]], size = 7.5, margin = margin(l = 6, r = 2)),
    axis.text.x = element_text(size = 7)
  )

save_figure(figure_E1, "output/figures/figure_E1.png", height_cm = 12.5)

## Annex Figure E2 ----
# The RDS estimate and 95% CI after each respondent in the realised recruitment
# order, Phase 1, with the final estimate and the household survey's estimate
# and 95% CI; the zero-dose panels have their own scale.
# legend entries and their colours
series_names <- c(rds = "RDS estimate and 95% CI", final = "RDS final estimate", hh = "HH estimate and 95% CI")
series_colours <- set_names(c(method_colours[["rds"]], ink[["muted"]], method_colours[["hh"]]), series_names)

# realised paths in percent, with the final estimate
convergence <-
  stopping_rule$realised |>
  filter(sample %in% rds_samples) |>
  left_join(sample_parts(rds_samples), by = "sample") |>
  mutate(
    antigen = names(outcome_fields)[match(outcome, outcome_fields)],
    final = 100 * estimate[n == max(n)],
    .by = c(sample, outcome)
  ) |>
  transmute(
    country, antigen, n, estimate = 100 * estimate, lower = 100 * ci_lower, upper = 100 * ci_upper, final,
    country_f = factor(country_labels[country], levels = country_labels), outcome_f = factor(antigen, levels = primary_outcomes),
    series = series_names[["rds"]], series_final = series_names[["final"]]
  )

# household survey reference
convergence_hh <-
  estimates |>
  filter(phase == 1, method == "HH", antigen %in% primary_outcomes) |>
  transmute(
    country_f = factor(country_labels[country], levels = country_labels), outcome_f = factor(as.character(antigen), levels = primary_outcomes),
    estimate = 100 * estimate, lower = 100 * ci_lower, upper = 100 * ci_upper, series = series_names[["hh"]]
  )

# one block for the antigens (0-100%) and one for zero-dose (0-50%)
convergence_block <- function(zero_dose, limit, step, first) {
  pick <- \(x) filter(x, (outcome_f == "Zero-dose") == zero_dose)
  ggplot() +
    geom_rect(data = pick(convergence_hh), aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper, fill = series), alpha = 0.13) +
    geom_hline(data = pick(convergence_hh), aes(yintercept = estimate, colour = series), linewidth = 0.55) +
    geom_ribbon(data = pick(convergence), aes(n, ymin = lower, ymax = upper, fill = series), alpha = 0.15) +
    geom_line(data = pick(convergence), aes(n, final, colour = series_final), linewidth = 0.4) +
    geom_line(data = pick(convergence), aes(n, estimate, colour = series), linewidth = 0.6) +
    facet_grid(country_f ~ outcome_f, switch = "y") +
    scale_discrete_manual(c("colour", "fill"), values = series_colours, breaks = unname(series_names), name = NULL, guide = guide_legend(nrow = 1)) +
    scale_x_continuous(limits = c(0, max(convergence$n) * 1.02), breaks = c(0, 500, 1000), labels = scales::label_comma(), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, limit), breaks = seq(0, limit, step), labels = \(x) paste0(x, "%"), expand = c(0, 0), oob = scales::oob_keep) +
    labs(x = if (first) "RDS respondents recruited, in recruitment order" else NULL, y = NULL) +
    theme_report() +
    theme(panel.grid.major.x = element_blank(), strip.text.y.left = if (first) element_text(angle = 90) else element_blank())
}

# create convergence plot
figure_E2 <-
  wrap_plots(convergence_block(FALSE, 100, 25, TRUE), convergence_block(TRUE, 50, 10, FALSE), nrow = 1, widths = c(3, 1.05), guides = "collect") &
  theme(legend.position = "bottom")

save_figure(figure_E2, "output/figures/figure_E2.png", height_cm = 9.5)

## Annex Figure E3 ----
# The recruitment chains of each RDS sample, one panel per country and phase:
# chains left to right by size, each node at its wave, seeds as squares, and
# respondents filled when their child had received Penta3.
# panel labels and the grey of the no-value marks
panel_label <- \(country, phase) str_c(country_labels[country], " · Phase ", phase)
no_value_grey <- "#D6D6DC"

# lay out each forest with chains placed left to right, scaled to 0-1
tree_nodes <- imap(networks, \(network, sample) {
  x <- network$nodes
  parts <- sample_parts(sample)
  has_parent <- x$id_referring %in% x$id_participant_calc
  graph <- graph_from_data_frame(
    data.frame(from = x$id_referring[has_parent], to = x$id_participant_calc[has_parent]),
    vertices = data.frame(name = x$id_participant_calc), directed = TRUE
  )
  roots <- x$id_participant_calc[!has_parent][order(x$chain_rank[!has_parent])]
  layout <- layout_as_tree(graph, root = match(roots, V(graph)$name), mode = "out")
  x$lx <- layout[, 1]
  x$ly <- -x$wave
  offset <- 0
  for (root in roots) {
    keep <- x$chain == x$chain[match(root, x$id_participant_calc)]
    width <- diff(range(x$lx[keep]))
    x$lx[keep] <- x$lx[keep] - min(x$lx[keep]) + offset
    offset <- offset + width + 2.2
  }
  x$lx <- x$lx / max(x$lx)
  x$panel <- panel_label(parts$country, parts$phase)
  x$state <- case_when(!x$interviewed | is.na(x$penta_3p_clean) ~ "none", x$penta_3p_clean == 1 ~ "received", x$penta_3p_clean == 0 ~ "not_received", .default = "none")
  x$key <- str_c(sample, x$id_participant_calc)
  x$parent_key <- str_c(sample, x$id_referring)
  x
}) |>
  list_rbind() |>
  mutate(panel = factor(panel, levels = unique(panel)))

# lines from each recruiter to their recruits
tree_edges <-
  tree_nodes |>
  filter(parent_key %in% key) |>
  mutate(px = tree_nodes$lx[match(parent_key, tree_nodes$key)], py = tree_nodes$ly[match(parent_key, tree_nodes$key)])

# key drawn beneath the panels
tree_key <- function() {
  label <- \(x, text) annotate("text", x = x, y = 0, label = text, hjust = 0, family = report_font, size = 2.4, colour = ink[["primary"]])
  ggplot() +
    annotate("point", x = 0, y = 0, shape = 22, fill = ink[["primary"]], colour = "white", size = 2.2, stroke = 0.3) +
    label(0.25, "Seed") +
    annotate("point", x = 2, y = 0, shape = 21, fill = method_colours[["rds"]], colour = method_colours[["rds"]], size = 1.6) +
    label(2.25, "Penta3 received") +
    annotate("point", x = 5, y = 0, shape = 21, fill = "white", colour = method_colours[["rds"]], size = 1.6, stroke = 0.5) +
    label(5.25, "Penta3 not received") +
    annotate("point", x = 8.4, y = 0, shape = 21, fill = no_value_grey, colour = no_value_grey, size = 1.4) +
    label(8.65, "No Penta3 record, or recruiter outside the analysis") +
    scale_x_continuous(limits = c(-0.3, 16)) +
    theme_void()
}

# create recruitment forest
nodes_in <- \(state) filter(tree_nodes, state == .env$state)
figure_E3_forest <-
  ggplot() +
  geom_segment(data = tree_edges, aes(x = px, y = py, xend = lx, yend = ly), colour = "#C9C7D6", linewidth = 0.18) +
  geom_point(data = nodes_in("not_received"), aes(lx, ly), shape = 21, fill = "white", colour = method_colours[["rds"]], size = 0.75, stroke = 0.3) +
  geom_point(data = nodes_in("received"), aes(lx, ly), shape = 21, fill = method_colours[["rds"]], colour = method_colours[["rds"]], size = 0.75, stroke = 0.3) +
  geom_point(data = nodes_in("none"), aes(lx, ly), shape = 21, fill = no_value_grey, colour = no_value_grey, size = 0.6) +
  geom_point(data = filter(tree_nodes, seed_participant %in% "Yes"), aes(lx, ly), shape = 22, fill = ink[["primary"]], colour = "white", size = 1.5, stroke = 0.3) +
  facet_wrap(~panel, ncol = 1, scales = "free") +
  scale_y_continuous(breaks = -seq(0, 30, 5), labels = seq(0, 30, 5), expand = expansion(add = 0.8)) +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.01)) +
  labs(x = NULL, y = "Wave") +
  theme_report() +
  theme(strip.text = element_text(hjust = 0, face = "bold"), panel.grid.major.y = element_line(colour = "#F0F0F3", linewidth = 0.3))

figure_E3 <- wrap_plots(figure_E3_forest, tree_key(), ncol = 1, heights = c(40, 1))

save_figure(figure_E3, "output/figures/figure_E3.png", height_cm = 17.5)

## Annex Figure E4 ----
# Cumulative Penta3 within chains in interview order: chains of 20 or more
# respondents in colour with a point at their final value, chains of 5 to 19
# in grey, the whole sample in black, and the reported RDS estimate dashed.
# chains by length, and the panel of each sample
trajectories <-
  chain_trajectories |>
  mutate(
    panel = factor(panel_label(country, phase), levels = panel_label(c("eth", "som", "eth", "som"), c(1, 1, 2, 2))),
    group = str_c(sample, chain),
    length = max(n),
    .by = c(sample, chain)
  )
long_chains <- filter(trajectories, chain != "all", length >= 20)
short_chains <- filter(trajectories, chain != "all", length >= 5, length < 20)
chain_ends <- filter(long_chains, n == length)
whole_sample <- filter(trajectories, chain == "all")

# reported RDS estimates
rds_reference <-
  estimates |>
  filter(method == "RDS", antigen == "Penta3") |>
  mutate(panel = factor(panel_label(country, phase), levels = levels(trajectories$panel)))

# create chain trajectory plot
figure_E4 <-
  ggplot() +
  geom_hline(data = rds_reference, aes(yintercept = 100 * estimate), colour = ink[["primary"]], linewidth = 0.4, linetype = "22") +
  geom_line(data = short_chains, aes(n, 100 * estimate, group = group), colour = "#D9D8E3", linewidth = 0.3) +
  geom_line(data = long_chains, aes(n, 100 * estimate, group = group), colour = method_colours[["rds"]], linewidth = 0.45, alpha = 0.85) +
  geom_point(data = chain_ends, aes(n, 100 * estimate), colour = method_colours[["rds"]], size = 0.9) +
  geom_line(data = whole_sample, aes(n, 100 * estimate, group = group), colour = ink[["primary"]], linewidth = 0.7) +
  facet_wrap(~panel, ncol = 2) +
  scale_y_continuous(limits = c(0, 100), labels = \(x) paste0(x, "%")) +
  scale_x_continuous(labels = scales::label_comma()) +
  labs(x = "Respondents interviewed (cumulative, in interview order)", y = "Penta3 coverage, cumulative estimate") +
  theme_report() +
  theme(strip.text = element_text(hjust = 0, face = "bold"))

save_figure(figure_E4, "output/figures/figure_E4.png", height_cm = 11)
