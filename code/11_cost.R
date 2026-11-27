# 11_cost.R
# Cost and resource requirements of each method, the cost at the Penta3
# stopping point, and the cost detail of Annex L
# Reproduces: Table 14, Figures 17, 18 and 21, Annex Tables L1-L9, Annex
# Figure L1

# Setup ----

pacman::p_load(ggh4x, ggtext, tidyverse)

source("code/functions.R")

# Data ----

cost_line_items <- read_rds("data/clean/cost_line_items.rds")
fx_cpi <- read_rds("data/clean/fx_cpi.rds")
frames <- read_rds("data/clean/frames.rds")
hours <- read_rds("data/clean/hours.rds")
stopping_rule <- read_rds("data/clean/stopping.rds")
estimates <- read_rds("data/clean/estimates.rds")
tallies <- read_rds("data/clean/tallies.rds")
record_counts <- read_rds("data/clean/kobo_record_counts.rds")
cost_shared_items <- read_rds("data/clean/cost_shared_items.rds")
research_team_cost <- read_rds("data/clean/research_team_cost.rds")

# the fourteen samples, in report order
samples <- study_samples

# Analysis ----

## Constant 2026 US dollars ----
# Birr are converted at the mean rate over the phase's fieldwork in Ethiopia,
# Swiss francs at the mean rate of the calendar year they were budgeted in,
# and Phase 1 (2025) amounts are taken to 2026 prices with the US CPI-U.
# birr per US dollar by phase, US dollars per Swiss franc by calendar year,
# and the 2025-to-2026 price index ratio
rates <- function(quantity, country) {
  rows <- filter(fx_cpi, quantity == .env$quantity, country == .env$country)
  set_names(rows$value, rows$phase_or_period)
}
etb_per_usd <- rates("etb_per_usd", "eth")
usd_per_chf <- rates("usd_per_chf", "all")
cpi_inflator <- filter(fx_cpi, quantity == "inflator_2025_to_2026")$value

# convert every line
cost_line_items <-
  cost_line_items |>
  mutate(
    to_usd = case_when(
      currency == "ETB" ~ 1 / etb_per_usd[as.character(study_phase)],
      currency == "CHF" ~ usd_per_chf[fx_period],
      .default = 1
    ),
    amount_usd = amount * unname(to_usd) * if_else(study_phase == 1, cpi_inflator, 1)
  )

## Cost per sample ----
# The programmatic cost of each sample (the cost of running the method; the
# research team's own lines are left out), split into a fixed part and a
# marginal part that grows with the number of respondents: the variable lines
# other than training. The tablets are a fixed cost.
# respondents who reported directly on their own children
respondents <- map_int(set_names(samples), \(sample) nrow(analytic_respondents(frames[[sample]], sample)))

# programmatic, variable and marginal cost per sample
cost <-
  cost_line_items |>
  filter(cost_scope == "programmatic") |>
  summarise(
    total_cost = sum(amount_usd),
    variable_cost = sum(amount_usd[cost_variability == "variable"]),
    variable_training = sum(amount_usd[cost_variability == "variable" & cost_category == "training"]),
    .by = c(country, study_phase, method)
  ) |>
  transmute(
    sample = paste0(country, "_phase", study_phase, "_", method), country, phase = study_phase, method, total_cost,
    marginal_cost = variable_cost - variable_training,
    fixed_cost = total_cost - marginal_cost,
    n_respondents = respondents[sample],
    cost_per_respondent = total_cost / n_respondents,
    marginal_cost_per_respondent = marginal_cost / n_respondents
  )

# ratio of total cost to the household survey of the same country (Phase 1)
cost <-
  cost |>
  mutate(total_cost_ratio = total_cost / first(total_cost[method == "hh"]), .by = c(country, phase))

## Field effort ----
# Total team hours including travel: enumerator time on site and travelling,
# and supervisor time. CBE recorded no supervisor time.
# team hours and their ratio to the household survey (Phase 1)
effort <-
  team_hours(hours) |>
  mutate(team_hours_ratio = total_hours / first(total_hours[method == "hh"]), .by = c(country, phase))

## Interview duration ----
# Median and interquartile range of the core interview, from the end of
# consent to the start of the network module, among directly reporting
# respondents with a valid duration. CBE was recorded on paper, without
# timestamps.
# summarise the core interview of each timestamped sample
duration <-
  str_subset(samples, "cbe", negate = TRUE) |>
  set_names() |>
  map(\(sample) {
    minutes <- core_interview_minutes(analytic_respondents(frames[[sample]], sample))
    tibble(
      n_valid = sum(!is.na(minutes)),
      median = median(minutes, na.rm = TRUE),
      p25 = quantile(minutes, 0.25, na.rm = TRUE, names = FALSE),
      p75 = quantile(minutes, 0.75, na.rm = TRUE, names = FALSE)
    )
  }) |>
  list_rbind(names_to = "sample")

## Cost at the stopping point ----
# Fixed cost plus the marginal cost per respondent times the respondents at
# the Penta3 stopping point (and at the 10th and 90th percentiles over
# enrolment orders), and its ratio to the household survey of the same
# country.
# cost for n respondents
cost_at <- function(n, fixed, marginal) fixed + marginal * n

# cost at the stopping point of each Phase 1 sample (median and 10th and 90th
# percentiles), and its ratio to the household survey
stopping <-
  stopping_rule$stopping |>
  left_join(select(cost, sample, fixed_cost, marginal_cost_per_respondent), by = "sample") |>
  mutate(
    cost_at_stopping = cost_at(stopping_n, fixed_cost, marginal_cost_per_respondent),
    cost_at_stopping_q10 = cost_at(stopping_n_q10, fixed_cost, marginal_cost_per_respondent),
    cost_at_stopping_q90 = cost_at(stopping_n_q90, fixed_cost, marginal_cost_per_respondent),
    sampling_unit = c(hh = "clusters", cbe = "clusters", rds = "respondents", vbs = "venue-days")[method]
  ) |>
  mutate(cost_ratio = cost_at_stopping / first(cost_at_stopping[method == "hh"]), .by = country)

## Cost decomposition ----
# Training is the variable training cost (left out of marginal cost). The
# additive total counts each shared cost once: the lines a partner booked once
# for a country and phase are divided across its methods by their share of
# directly reporting respondents, and the tablet setting across the methods
# using tablets (HH, RDS and VBS) the same way.
# shared lines and the tablet setting of each country and phase
shared_cost <-
  cost_line_items |>
  semi_join(cost_shared_items, by = c("country", "study_phase", "organisation", "item", "amount")) |>
  distinct(country, study_phase, organisation, item, amount, amount_usd) |>
  summarise(shared = sum(amount_usd), .by = c(country, study_phase))
tablet_setting <-
  cost_line_items |>
  filter(item == "Tablets (capital equipment)") |>
  summarise(tablet = first(amount_usd), .by = c(country, study_phase))

# training cost and the additive total of each sample
decomposition <-
  cost_line_items |>
  filter(cost_scope == "programmatic") |>
  summarise(training = sum(amount_usd[cost_variability == "variable" & cost_category == "training"]), .by = c(country, study_phase, method)) |>
  right_join(cost, by = c("country", "study_phase" = "phase", "method")) |>
  left_join(shared_cost, by = c("country", "study_phase")) |>
  left_join(tablet_setting, by = c("country", "study_phase")) |>
  mutate(
    device = method != "cbe",
    share = n_respondents / sum(n_respondents),
    device_share = if_else(device, n_respondents / sum(n_respondents[device]), 0),
    .by = c(country, study_phase)
  ) |>
  mutate(
    additive_total = total_cost - shared * (1 - share) - if_else(device, tablet * (1 - device_share), 0),
    tablet_share = tablet * device_share
  ) |>
  rename(phase = study_phase)

## Research-team time ----
# The UNIGE research team's Phase 1 hours at their day rates, with the 13%
# overhead, in 2026 dollars. Hours logged for both countries are split equally
# between them, and hours logged for no single method equally across the four.
# cost of each country and method logged, apportioned to country and method
research_time <-
  research_team_cost |>
  mutate(
    cost = day_rate_cost_usd * 1.13 * cpi_inflator,
    share = if_else(country == "both", 0.5, 1) * if_else(method == "cross-cutting", 0.25, 1)
  ) |>
  left_join(tibble(country = "both", to_country = c("eth", "som")), by = "country", relationship = "many-to-many") |>
  left_join(tibble(method = "cross-cutting", to_method = method_order), by = "method", relationship = "many-to-many") |>
  mutate(country = coalesce(to_country, country), method = coalesce(to_method, method)) |>
  summarise(research_time = sum(cost * share), .by = c(country, method))

## Study cost ----
# Each cost counted once: the partner budgets (additive totals without the
# tablet setting, whose actual cost is UNIGE's tablet line), UNIGE's non-staff
# lines, divided equally between the countries, and the research-team time.
# UNIGE lines in report order, with their labels
unige_labels <- c(
  Consultants = "Consultants", `Admin assistant` = "Administrative assistant", Avion = "International flights",
  `Hébergement et repas` = "Accommodation and meals", Tablette = "Tablets", Assurance = "Insurance",
  `Impressions/fournitures` = "Printing and supplies", Visa = "Visas",
  `Internet & mobile phone` = "Communications (internet, mobile telephone and SIM cards)",
  `Mobile SIM cards` = "Communications (internet, mobile telephone and SIM cards)", `Location véhicule` = "Vehicle hire"
)

# each component by country, the UNIGE lines indented as printed
unige_lines <-
  cost_line_items |>
  filter(organisation == "UNIGE") |>
  summarise(amount = sum(amount_usd), .by = c(country, item)) |>
  mutate(component = unige_labels[item]) |>
  summarise(amount = sum(amount), .by = c(country, component)) |>
  arrange(country, desc(amount))
study_cost <- bind_rows(
  decomposition |> filter(phase == 1) |> summarise(amount = sum(additive_total - tablet_share), .by = country) |>
    mutate(component = "Partner budget (AAU in Ethiopia, SORDI in Somalia), Phase 1"),
  decomposition |> filter(phase == 2) |> summarise(amount = sum(additive_total - tablet_share), .by = country) |>
    mutate(component = "Partner budget, Phase 2"),
  unige_lines |> summarise(amount = sum(amount), .by = country) |> mutate(component = "UNIGE non-staff budget lines, Phase 1"),
  unige_lines |> mutate(component = str_c("\u200b ", component)),
  research_time |> summarise(amount = sum(research_time), .by = country) |> mutate(component = "UNIGE research-team time, Phase 1")
)

## Unit costs ----
# Marginal and total cost per sampling unit (clusters; RDS seeds and recruits,
# the respondents other than seeds; venues with a directly reporting
# caregiver), per contact and per directly reporting respondent.
# sampling units and contacts of each sample
sampling_units <- map(samples, \(sample) {
  method <- str_extract(sample, "[a-z]+$")
  cohort <- direct_cohort(frames[[sample]], sample)
  switch(
    method,
    hh = tibble(units = n_distinct(cohort$cluster_id), recruits = NA_integer_),
    cbe = tibble(units = n_distinct(cohort$cluster), recruits = NA_integer_),
    rds = distinct(cohort, record_ix, .keep_all = TRUE) |> summarise(units = sum(seed_participant %in% "Yes"), recruits = sum(seed_participant %in% "No")),
    vbs = tibble(units = n_distinct(cohort$int_location), recruits = NA_integer_)
  )
}) |>
  set_names(samples) |>
  list_rbind(names_to = "sample") |>
  left_join(sample_contacts(frames, tallies, record_counts), by = "sample")

## Cost and precision ----
# Kish's ratio of the cost of a unit of precision in the household survey to
# that in each method, Phase 1: marginal cost per respondent times the
# variance per respondent (design effect times p(1 - p) at the household
# survey's estimate, which cancels), and with each method's own variance per
# respondent (children observed times the squared standard error).
# design effect and variance per respondent of the three primary antigens
precision <-
  estimates |>
  filter(phase == 1, antigen %in% c("Penta3", "OPV3", "MCV1")) |>
  transmute(
    sample, country, method = str_to_lower(as.character(method)), antigen = as.character(antigen),
    deff = standard_error^2 / (estimate * (1 - estimate) / n_observed),
    own_variance = n_observed * standard_error^2
  ) |>
  left_join(select(cost, sample, marginal_cost_per_respondent), by = "sample") |>
  mutate(
    ratio = first((marginal_cost_per_respondent * deff)[method == "hh"]) / (marginal_cost_per_respondent * deff),
    ratio_own = first((marginal_cost_per_respondent * own_variance)[method == "hh"]) / (marginal_cost_per_respondent * own_variance),
    .by = c(country, antigen)
  )

## Stopping conditions ----
# The sample size at which each condition was first met on its own and both
# together, for Penta3, and both for all three primary antigens (the latest
# of the three in each enrolment order), as medians and 10th to 90th
# percentiles over orders; an order that never meets a condition counts at
# the full sample. The cost ratio is at the all-three stopping point.
# stopping point for all three antigens in each order
all_three <-
  stopping_rule$order_stopping |>
  summarise(stopping_all = if (anyNA(stopping_n)) NA_real_ else max(stopping_n), .by = c(sample, order))

# median and range of each condition, and the cost at the all-three point
conditions <-
  stopping_rule$order_stopping |>
  filter(outcome == "penta_3p_clean") |>
  left_join(all_three, by = c("sample", "order")) |>
  left_join(select(stopping, sample, n_final), by = "sample") |>
  summarise(
    stability_median = order_percentile(stability_n, 0.5, n_final[1]),
    stability_q10 = order_percentile(stability_n, 0.1, n_final[1]),
    stability_q90 = order_percentile(stability_n, 0.9, n_final[1]),
    precision_median = order_percentile(precision_n, 0.5, n_final[1]),
    precision_q10 = order_percentile(precision_n, 0.1, n_final[1]),
    precision_q90 = order_percentile(precision_n, 0.9, n_final[1]),
    all_median = order_percentile(stopping_all, 0.5, n_final[1]),
    all_q10 = order_percentile(stopping_all, 0.1, n_final[1]),
    all_q90 = order_percentile(stopping_all, 0.9, n_final[1]),
    .by = sample
  ) |>
  left_join(select(cost, sample, country, method, fixed_cost, marginal_cost_per_respondent), by = "sample") |>
  mutate(
    cost_all = cost_at(all_median, fixed_cost, marginal_cost_per_respondent),
    cost_ratio_all = cost_all / first(cost_all[method == "hh"]),
    .by = country
  )

## Tablet cost annualised ----
# The flat tablet setting replaced by the equivalent annual cost of the
# devices over a useful life of 3 years at a 3% discount rate. One pool of
# devices per country (Ethiopia 16 bought at CHF 99 each, Somalia 15 at USD
# 125) served HH, RDS and VBS in both phases; each method is charged for its
# device-days, the enumerators and supervisors using tablets summed over the
# calendar days of its fieldwork, with each day's devices shared in proportion
# to users where methods together exceeded the pool. CBE used paper.
# replacement value of a device in 2026 dollars, and its equivalent annual cost
devices <- tibble(
  country = c("eth", "som"), n_devices = c(16, 15),
  unit_price = c(99 * usd_per_chf[["cal_2025"]], 125) * cpi_inflator
)
useful_life <- 3
discount_rate <- 0.03
annuity_factor <- (1 - (1 + discount_rate)^-useful_life) / discount_rate
devices$annual_cost <- devices$unit_price / annuity_factor

# fieldwork days and tablet users of each device method
device_use <-
  map(str_subset(samples, "cbe", negate = TRUE), \(sample) {
    dates <- as.Date(distinct(frames[[sample]], record_ix, .keep_all = TRUE)$date)
    tibble(sample, start = min(dates, na.rm = TRUE), end = max(dates, na.rm = TRUE))
  }) |>
  list_rbind() |>
  left_join(sample_parts(samples), by = "sample") |>
  left_join(staff_counts(hours), by = c("country", "phase", "method")) |>
  mutate(users = n_enumerators + n_supervisors) |>
  left_join(select(devices, country, n_devices), by = "country")

# device-days, sharing each day's pool among the methods in the field
device_days <-
  device_use |>
  mutate(day = map2(start, end, \(a, b) seq(a, b, by = "day"))) |>
  unnest(day) |>
  mutate(devices = users * pmin(1, n_devices / sum(users)), .by = c(country, phase, day)) |>
  summarise(device_days = sum(devices), .by = c(country, phase, method))

# annualised charge, adjusted totals and costs at the stopping point
annualised <-
  decomposition |>
  left_join(device_days, by = c("country", "phase", "method")) |>
  left_join(select(devices, country, annual_cost), by = "country") |>
  left_join(select(stopping, sample, stopping_basis, cost_at_stopping, cost_ratio), by = "sample") |>
  mutate(
    device_days = coalesce(device_days, 0),
    flat = if_else(device, tablet, 0),
    charge = annual_cost / 365 * device_days,
    total_adjusted = total_cost - flat + charge,
    stopping_adjusted = cost_at_stopping - flat + charge
  ) |>
  mutate(ratio_adjusted = stopping_adjusted / first(stopping_adjusted[method == "hh"]), .by = c(country, phase))

# cost verdict of a ratio to the household survey
cost_verdict <- function(ratio) {
  case_when(is.na(ratio) ~ NA_character_, ratio <= 1.10 + 1e-9 ~ "Meets", ratio <= 1.25 + 1e-9 ~ "Borderline", .default = "Fails")
}

## Planning calculation ----
# Respondents required for a 95% CI of a given half-width at a prevalence of
# 50%, using each Phase 1 sample's own design effect (the mean over Penta3,
# OPV3 and MCV1) and children per respondent, and their fixed-plus-marginal
# cost.
# target margins of error (pp)
target_moe <- c(15, 10, 7.5, 5, 3)

# respondents and cost at each margin
planning <-
  estimates |>
  filter(phase == 1, antigen %in% c("Penta3", "OPV3", "MCV1")) |>
  summarise(deff = mean(standard_error^2 / (estimate * (1 - estimate) / n_observed)), .by = sample) |>
  left_join(cost, by = "sample") |>
  mutate(children_per_respondent = map_int(sample, \(s) nrow(direct_cohort(frames[[s]], s))) / n_respondents) |>
  cross_join(tibble(moe = target_moe)) |>
  mutate(
    n_children = ceiling(1.96^2 * deff * 0.25 / (moe / 100)^2),
    n_required = ceiling(n_children / children_per_respondent),
    total = fixed_cost + marginal_cost_per_respondent * n_required
  )

# Tables ----

## Table 14 ----
# Cost, cost per directly reporting respondent, team hours, core interview
# duration and the stopping point, by phase, country and method. Ratios are to
# the household survey in Phase 1. A double dagger marks the two costs per
# respondent based on fewer than 100 directly reporting respondents, which
# the table footnote gives.
# assemble one row per sample in report order
table_14_rows <-
  cost |>
  left_join(select(effort, country, phase, method, total_hours, team_hours_ratio), by = c("country", "phase", "method")) |>
  left_join(duration, by = "sample") |>
  left_join(
    select(stopping, sample, stopping_basis, stopping_n, stopping_n_q10, stopping_n_q90, n_units_at_stopping, sampling_unit,
           cost_at_stopping, cost_at_stopping_q10, cost_at_stopping_q90, cost_ratio),
    by = "sample"
  ) |>
  mutate(method = factor(method, levels = method_order)) |>
  arrange(phase, country, method) |>
  mutate(
    Country = if_else(row_number() == 1, country_labels[country], ""),
    .by = c(phase, country)
  )

# format the cells as printed
table_14 <-
  table_14_rows |>
  transmute(
    phase,
    Country,
    Method = method_labels[as.character(method)],
    `Total cost ($)` = format_dollars(total_cost),
    `Ratio to HH` = format_number(total_cost_ratio, 2),
    `Cost per directly reporting respondent ($)` = str_c(format_dollars(cost_per_respondent), if_else(n_respondents < 100, "‡", "")),
    `Total team hours` = format_count(total_hours),
    `Ratio to HH ` = format_number(team_hours_ratio, 2),
    `Median core interview, min (IQR)` = if_else(is.na(median), "–", str_glue("{format_number(median)} ({format_number(p25)}–{format_number(p75)})")),
    `Stopping N, median (10th–90th)` = case_when(
      stopping_basis == "observed" ~ str_glue("{format_count(stopping_n)} ({format_count(stopping_n_q10)}–{format_count(stopping_n_q90)})"),
      stopping_basis == "modelled" ~ str_glue("{format_count(stopping_n)} ({n_units_at_stopping} {sampling_unit})"),
      .default = "–"
    ),
    `Stopping basis` = replace_na(str_to_sentence(str_replace(stopping_basis, "_", " ")), "–"),
    `Cost at the stopping point ($), median (10th–90th)` = case_when(
      stopping_basis == "observed" ~ str_glue("{format_dollars(cost_at_stopping)} ({format_dollars(cost_at_stopping_q10)}–{format_dollars(cost_at_stopping_q90)})"),
      .default = format_dollars(cost_at_stopping)
    ),
    `Ratio to HH  ` = format_number(cost_ratio, 2)
  )

# add a heading row for each phase
table_14 <-
  table_14 |>
  group_by(phase) |>
  group_modify(\(rows, key) bind_rows(tibble(Country = str_c("Phase ", key$phase)), rows)) |>
  ungroup() |>
  select(-phase) |>
  mutate(across(everything(), \(x) replace_na(as.character(x), "")))

# the three ratio columns share the report's heading
names(table_14) <- str_trim(names(table_14))

save_table(table_14, "output/tables/table_14.csv")

## Annex Table L1 ----
# Total, fixed, marginal and training cost, the additive total and the UNIGE
# research-team time by sample (Phase 1 only).
# "1, Ethiopia" row labels
# one row per sample in report order
table_L01 <-
  decomposition |>
  left_join(research_time, by = c("country", "method")) |>
  mutate(method = factor(method, levels = method_order)) |>
  arrange(phase, country, method) |>
  transmute(
    `Phase, country` = phase_country(phase, country),
    Method = method_labels[as.character(method)],
    Total = format_count(total_cost),
    Fixed = format_count(fixed_cost),
    Marginal = format_count(marginal_cost),
    Training = format_count(training),
    `Additive total` = format_count(additive_total),
    `UNIGE research-team time` = if_else(phase == 1, format_count(research_time), "not built")
  )

save_table(table_L01, "output/tables/table_L01.csv")

## Annex Table L2 ----
# Sample size at which each stopping condition was first met (median and 10th
# to 90th percentile over enrolment orders), Phase 1, and the cost ratio at the
# all-three-antigen stopping point.
# a median with its range, or "Not met"
with_range <- \(median, q10, q90) if_else(is.na(median), "Not met", str_glue("{format_count(median)} ({format_count(q10)}–{format_count(q90)})"))

# one row per sample, the country on its first row
table_L02 <-
  conditions |>
  left_join(select(stopping, sample, n_final, stopping_basis, stopping_n, stopping_n_q10, stopping_n_q90), by = "sample") |>
  mutate(method = factor(method, levels = method_order)) |>
  arrange(country, method) |>
  mutate(Country = if_else(row_number() == 1, country_labels[country], ""), .by = country) |>
  transmute(
    Country,
    Method = method_labels[as.character(method)],
    `Final sample` = format_count(n_final),
    `Stability alone, Penta3` = with_range(stability_median, stability_q10, stability_q90),
    `Precision alone, Penta3` = with_range(precision_median, precision_q10, precision_q90),
    `Stopping point, Penta3` = if_else(stopping_basis == "modelled", str_c(format_count(stopping_n), " (modelled)"), with_range(stopping_n, stopping_n_q10, stopping_n_q90)),
    `Stopping point, all three primary antigens` = with_range(all_median, all_q10, all_q90),
    `Cost ratio to HH, all three primary antigens` = if_else(is.na(cost_ratio_all), "Not met", format_number(cost_ratio_all, 2))
  )

save_table(table_L02, "output/tables/table_L02.csv")

## Annex Table L3 ----
# Cost of the study by component and country, each cost counted once.
# components in report order (UNIGE lines indented), then the total
study_cost_wide <-
  study_cost |>
  pivot_wider(names_from = country, values_from = amount) |>
  mutate(both = eth + som)
study_total <-
  study_cost_wide |>
  filter(!str_starts(component, "\u200b")) |>
  summarise(across(c(eth, som, both), sum)) |>
  mutate(component = "Total")
table_L03 <-
  bind_rows(study_cost_wide, study_total) |>
  transmute(Component = component, Ethiopia = format_count(eth), Somalia = format_count(som), `Both countries` = format_count(both))

save_table(table_L03, "output/tables/table_L03.csv")

## Annex Table L4 ----
# Unit costs at the realised sample: marginal and total cost per sampling
# unit (RDS total per recruit only), per contact and per directly reporting
# respondent; a dagger marks fewer than 100 respondents.
# one row per sample in report order
table_L04 <-
  cost |>
  left_join(sampling_units, by = "sample") |>
  mutate(method = factor(method, levels = method_order), rds = method == "rds") |>
  arrange(phase, country, method) |>
  transmute(
    `Phase, country` = phase_country(phase, country),
    Method = method_labels[as.character(method)],
    `Sampling unit` = c(hh = "Cluster", cbe = "Cluster", rds = "Seed; recruit", vbs = "Venue")[as.character(method)],
    Units = if_else(rds, str_c(format_count(units), "; ", format_count(recruits)), format_count(units)),
    `Marginal cost per unit` = if_else(rds, str_c(format_count(marginal_cost / units), "; ", format_count(marginal_cost / recruits)), format_count(marginal_cost / units)),
    `Total cost per unit` = if_else(rds, str_c("–; ", format_count(total_cost / recruits)), format_count(total_cost / units)),
    Contacts = format_count(n_contacts),
    `Marginal cost per contact` = format_number(marginal_cost / n_contacts, 2),
    `Total cost per contact` = format_number(total_cost / n_contacts, 2),
    `Directly reporting respondents` = format_count(n_respondents),
    `Marginal cost per respondent` = str_c(format_count(marginal_cost_per_respondent), if_else(n_respondents < 100, "†", ""))
  )

save_table(table_L04, "output/tables/table_L04.csv")

## Annex Table L5 ----
# Cost at the Penta3 stopping point per respondent, per sampling unit (cluster
# or venue-day; none for RDS, whose seeds were all recruited at the start) and
# per contact, with units and contacts scaled from the realised sample; an
# asterisk marks a modelled stopping point.
# one row per Phase 1 sample
table_L05 <-
  stopping |>
  left_join(select(cost, sample, n_respondents), by = "sample") |>
  left_join(select(sampling_units, sample, n_contacts), by = "sample") |>
  mutate(method = factor(method, levels = method_order)) |>
  arrange(country, method) |>
  transmute(
    Country = country_labels[country],
    Method = method_labels[as.character(method)],
    `Stopping point (respondents)` = str_c(format_count(stopping_n), if_else(stopping_basis == "modelled", "*", "")),
    `Cost at the stopping point` = format_count(cost_at_stopping),
    `Per respondent` = format_count(cost_at_stopping / stopping_n),
    `Per sampling unit` = if_else(method == "rds", "–", format_count(cost_at_stopping / n_units_at_stopping)),
    `Per contact` = format_number(cost_at_stopping / (n_contacts / n_respondents * stopping_n), 2)
  )

save_table(table_L05, "output/tables/table_L05.csv")

## Annex Table L6 ----
# Marginal cost per respondent, design effect and cost per effective
# respondent for Penta3, and the Kish ratio for each primary antigen and with
# each method's own Penta3 variance, Phase 1.
# one row per sample
table_L06 <-
  precision |>
  select(sample, country, method, antigen, marginal_cost_per_respondent, deff, ratio, ratio_own) |>
  pivot_wider(names_from = antigen, values_from = c(deff, ratio, ratio_own)) |>
  mutate(method = factor(method, levels = method_order)) |>
  arrange(country, method) |>
  transmute(
    Country = country_labels[country],
    Method = method_labels[as.character(method)],
    `Marginal cost per respondent (USD)` = format_number(marginal_cost_per_respondent, 2),
    `Design effect, Penta3` = format_number(deff_Penta3, 2),
    `Cost per effective respondent, Penta3 (USD)` = format_number(marginal_cost_per_respondent * deff_Penta3, 2),
    `Ratio, Penta3` = format_number(ratio_Penta3, 2),
    `Ratio, OPV3` = format_number(ratio_OPV3, 2),
    `Ratio, MCV1` = format_number(ratio_MCV1, 2),
    `Ratio with the method's own variance, Penta3` = format_number(ratio_own_Penta3, 2)
  )

save_table(table_L06, "output/tables/table_L06.csv")

## Annex Table L7 ----
# Exchange rates (means of daily mid-market rates) and the CPI-U price
# adjustment used to express costs in constant 2026 US dollars.
# a date window in words
window_text <- function(start, end) {
  start <- as.Date(start)
  end <- as.Date(end)
  if (format(start, "%m-%d") == "01-01" && format(end, "%m-%d") == "12-31") return(str_c("Calendar year ", format(end, "%Y")))
  str_c(as.integer(format(start, "%d")), " ", format(start, "%B"), " to ", as.integer(format(end, "%d")), " ", format(end, "%B %Y"))
}
rate_window <- \(quantity, period) filter(fx_cpi, quantity == .env$quantity, phase_or_period == period)

# one row per conversion
inflator_text <- format_number(cpi_inflator, 4)
table_L07 <- tibble(
  Partner = c("AAU", "AAU", "UNIGE", "UNIGE", "SORDI"),
  Currency = c("Ethiopian birr", "Ethiopian birr", "Swiss franc", "Swiss franc", "US dollar"),
  `Costs converted` = c("Ethiopia Phase 1", "Ethiopia Phase 2", "2025 budget lines, both countries", "2026 budget lines, both countries", "Somalia Phases 1 and 2"),
  `Period of the mean` = c(
    str_c("Fieldwork, ", with(rate_window("etb_per_usd", "1"), window_text(window_start, window_end))),
    str_c("Fieldwork, ", with(rate_window("etb_per_usd", "2"), window_text(window_start, window_end))),
    with(rate_window("usd_per_chf", "cal_2025"), window_text(window_start, window_end)),
    with(rate_window("usd_per_chf", "cal_2026"), window_text(window_start, window_end)),
    "–"
  ),
  Rate = c(
    str_c(format_number(etb_per_usd[["1"]], 2), " birr per USD"), str_c(format_number(etb_per_usd[["2"]], 2), " birr per USD"),
    str_c(format_number(usd_per_chf[["cal_2025"]], 4), " USD per CHF"), str_c(format_number(usd_per_chf[["cal_2026"]], 4), " USD per CHF"), "–"
  ),
  `Price adjustment to 2026` = c(inflator_text, "None", inflator_text, inflator_text, str_c(inflator_text, " (Phase 1); none (Phase 2)"))
)

save_table(table_L07, "output/tables/table_L07.csv")

## Annex Table L8 ----
# Cost at the Penta3 stopping point with the flat tablet setting and with the
# tablet cost annualised, Phase 1, and the verdict against the 1.10 and 1.25
# ratios to HH; an asterisk marks a modelled stopping point.
# one row per Phase 1 sample
table_L08 <-
  annualised |>
  filter(phase == 1) |>
  mutate(method = factor(method, levels = method_order), hh = method == "hh") |>
  arrange(country, method) |>
  transmute(
    Country = country_labels[country],
    Method = method_labels[as.character(method)],
    `Device-days` = format_count(device_days),
    `Annualised tablet charge` = format_count(charge),
    `Flat setting: ratio to HH` = format_number(cost_ratio, 2),
    `Flat setting: verdict` = if_else(hh, "–", cost_verdict(cost_ratio)),
    `Annualised: cost at the stopping point` = str_c(format_count(stopping_adjusted), if_else(stopping_basis == "modelled", "*", "")),
    `Annualised: ratio to HH` = format_number(ratio_adjusted, 2),
    `Annualised: verdict` = if_else(hh, "–", cost_verdict(ratio_adjusted))
  )

save_table(table_L08, "output/tables/table_L08.csv")

## Annex Table L9 ----
# Total cost and cost per directly reporting respondent with the flat tablet
# setting and annualised; a dagger marks fewer than 100 respondents.
# one row per sample in report order
table_L09 <-
  annualised |>
  mutate(method = factor(method, levels = method_order), few = if_else(n_respondents < 100, "†", "")) |>
  arrange(phase, country, method) |>
  transmute(
    `Phase, country` = phase_country(phase, country),
    Method = method_labels[as.character(method)],
    `Flat setting: total cost` = format_count(total_cost),
    `Annualised: total cost` = format_count(total_adjusted),
    `Flat setting: cost per respondent` = str_c(format_count(total_cost / n_respondents), few),
    `Annualised: cost per respondent` = str_c(format_count(total_adjusted / n_respondents), few)
  )

save_table(table_L09, "output/tables/table_L09.csv")

# Figures ----

## Figure 17 ----
# Programmatic cost of each sample divided by cost category, labelled with
# the total and the cost per directly reporting respondent.
# cost categories in stacking order, darkest first
cost_categories <- c(
  coordination_field = "Coordination and supervision",
  direct_field = "Field operations and equipment",
  enumerator = "Enumerators",
  training = "Training",
  incentive = "Participant incentives"
)

# cost by category, with country and phase as nested row strips; segments are
# drawn in stacking order so each white edge falls the same way
segments <-
  cost_line_items |>
  filter(cost_scope == "programmatic") |>
  summarise(amount_usd = sum(amount_usd), .by = c(country, study_phase, method, cost_category)) |>
  mutate(
    category = factor(cost_category, levels = names(cost_categories), labels = cost_categories),
    method = factor(method, levels = rev(method_order)),
    country_f = factor(country_labels[country], levels = country_labels),
    phase_f = factor(phase_labels[as.character(study_phase)], levels = phase_labels)
  ) |>
  arrange(category)

# total and cost per respondent of each sample, for the bar labels
totals <-
  cost |>
  mutate(
    method = factor(method, levels = rev(method_order)),
    country_f = factor(country_labels[country], levels = country_labels),
    phase_f = factor(phase_labels[as.character(phase)], levels = phase_labels),
    label = str_glue("{format_dollars(total_cost)} · {format_dollars(cost_per_respondent)} per respondent")
  )

# create stacked bar chart
figure_17 <-
  ggplot() +
  geom_col(
    data = segments, aes(amount_usd, method, fill = method, alpha = category),
    orientation = "y", width = 0.68, colour = "white", linewidth = 0.5, position = position_stack(reverse = TRUE)
  ) +
  geom_text(data = totals, aes(total_cost, method, label = label), hjust = 0, nudge_x = 500, size = 2.4, colour = ink[["primary"]], family = report_font) +
  facet_nested(
    rows = vars(country_f, phase_f), scales = "free_y", space = "free_y", switch = "y",
    strip = strip_nested(
      text_y = elem_list_text(angle = c(90, 90, 0, 0, 0, 0), face = "bold", hjust = c(0.5, 0.5, 1, 1, 1, 1)),
      by_layer_y = FALSE, clip = "off"
    )
  ) +
  scale_fill_manual(values = method_colours, guide = "none") +
  scale_alpha_manual(
    values = set_names(c(1, 0.72, 0.5, 0.33, 0.2), cost_categories), name = NULL,
    guide = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(fill = ink[["secondary"]], colour = "white"))
  ) +
  scale_y_discrete(labels = method_labels) +
  scale_x_continuous(labels = scales::label_dollar(), breaks = seq(0, 30000, 10000), limits = c(0, 41000), expand = c(0, 0)) +
  labs(x = "Programmatic cost (USD)", y = NULL) +
  theme_report() +
  theme(
    panel.grid.major.y = element_blank(), strip.placement = "outside",
    axis.text.y = element_text(colour = ink[["primary"]], size = 8, margin = margin(l = 8, r = 2))
  )

save_figure(figure_17, "output/figures/figure_17.png", height_cm = 12)

## Figure 18 ----
# Total cost, total team hours and cost at the Penta3 stopping point of CBE,
# RDS and VBS as ratios to the household survey (Phase 1), with the 1.10 and
# 1.25 standards. A modelled stopping point is an open, labelled square.
# measures and their marker shapes
resource_measures <- c(total_cost = "Total cost", team_hours = "Total team hours", stopping_cost = "Cost at the Penta3 stopping point")

# the three ratios of CBE, RDS and VBS, Phase 1
ratios <-
  bind_rows(
    cost |> filter(phase == 1) |> transmute(country, method, measure = "total_cost", ratio = total_cost_ratio, basis = "observed"),
    effort |> filter(phase == 1) |> transmute(country, method, measure = "team_hours", ratio = team_hours_ratio, basis = "observed"),
    stopping |> transmute(country, method, measure = "stopping_cost", ratio = cost_ratio, basis = stopping_basis)
  ) |>
  filter(method != "hh", !is.na(ratio)) |>
  mutate(
    measure = factor(measure, levels = names(resource_measures), labels = resource_measures),
    method_f = factor(method_labels[method], levels = method_labels[c("cbe", "rds", "vbs")]),
    method = factor(method, levels = c("vbs", "rds", "cbe")),
    country_f = factor(country_labels[country], levels = country_labels),
    modelled = basis == "modelled"
  )

# labels of the two standards, at the top of the CBE panels
standards <- crossing(
  tibble(method_f = factor("CBE", levels = levels(ratios$method_f)), x = c(1.10, 1.25), label = c("1.10", "1.25")),
  country_f = factor(country_labels, levels = country_labels)
)

# create dot plot
figure_18 <-
  ggplot(ratios, aes(ratio, measure)) +
  geom_vline(xintercept = c(1.10, 1.25), linetype = "dashed", colour = ink[["muted"]], linewidth = 0.4) +
  geom_vline(xintercept = 1, colour = method_colours[["hh"]], linewidth = 0.55) +
  geom_text(data = standards, aes(x = x, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.25, vjust = 1.4, size = 2.3, colour = ink[["muted"]], family = report_font) +
  geom_point(data = filter(ratios, !modelled), aes(shape = measure, fill = method), colour = "white", size = 3, stroke = 0.6) +
  geom_point(data = filter(ratios, modelled), aes(shape = measure, colour = method), fill = "white", size = 3, stroke = 0.9, show.legend = FALSE) +
  geom_text(data = filter(ratios, modelled), aes(label = "modelled"), vjust = -1.1, size = 2.3, family = report_font, colour = ink[["muted"]], fontface = "italic") +
  facet_grid(method_f ~ country_f, switch = "y") +
  scale_fill_manual(
    values = method_colours[c("cbe", "rds", "vbs")], breaks = c("cbe", "rds", "vbs"), labels = method_labels[c("cbe", "rds", "vbs")], name = "Method",
    guide = guide_legend(override.aes = list(shape = 21, size = 2.8, colour = "white"))
  ) +
  scale_colour_manual(values = method_colours[c("cbe", "rds", "vbs")], guide = "none") +
  scale_shape_manual(values = set_names(c(21, 24, 22), resource_measures), guide = "none") +
  scale_y_discrete(limits = rev, labels = scales::label_wrap(19)) +
  scale_x_continuous(limits = c(0.5, 1.5), breaks = seq(0.6, 1.4, 0.2), expand = c(0, 0)) +
  labs(x = "Ratio to the household survey", y = NULL) +
  theme_report() +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7, lineheight = 0.9))

save_figure(figure_18, "output/figures/figure_18.png", height_cm = 10)

## Figure 21 ----
# Penta3 along the realised enrolment order: the width of the 95% CI against
# the 10-point threshold, and the estimate minus the final estimate within a
# 5-point band; below, the stopping point (median and 10th to 90th percentile
# over random orders, or an open square past the end of fieldwork where it is
# modelled).
# panel rows, with a muted note under each title
trajectory_measures <- c(precision = "Precision", stability = "Stability", stopping = "Stopping point")
trajectory_notes <- c(
  precision = "95% CI width (pp),<br>realised enrolment order",
  stability = "estimate minus final<br>estimate (pp), realised order",
  stopping = "median and 10th to 90th percentile<br>over orders; open square: modelled"
)
strip_labels <- set_names(
  str_glue("**{trajectory_measures}**<br><span style='color:{ink[['muted']]};font-size:6.5pt'>{trajectory_notes}</span>"),
  trajectory_measures
)

# precision and stability lines along the realised order
lines <-
  stopping_rule$realised |>
  filter(outcome == "penta_3p_clean") |>
  mutate(final = estimate[n == max(n)], .by = sample) |>
  left_join(sample_parts(unique(stopping_rule$realised$sample)), by = "sample") |>
  reframe(
    tibble(measure = trajectory_measures[c("precision", "stability")], value = c(ci_width_pp, 100 * (estimate - final))),
    .by = c(sample, country, method, n)
  ) |>
  mutate(
    measure = factor(measure, levels = trajectory_measures),
    method = factor(method, levels = method_order),
    country_f = factor(country_labels[country], levels = country_labels)
  )

# stopping points, one row per method (HH at the top)
markers <-
  stopping |>
  mutate(
    measure = factor(trajectory_measures[["stopping"]], levels = trajectory_measures),
    method = factor(method, levels = method_order),
    country_f = factor(country_labels[country], levels = country_labels),
    slot = length(method_order) + 1 - as.integer(method)
  )
observed_markers <- filter(markers, stopping_basis == "observed")
modelled_markers <- filter(markers, stopping_basis == "modelled")

# reference lines and band
precision_line <- tibble(measure = factor(trajectory_measures[["precision"]], levels = trajectory_measures), y = 10)
stability_zero <- tibble(measure = factor(trajectory_measures[["stability"]], levels = trajectory_measures), y = 0)
stability_band <- tibble(measure = factor(trajectory_measures[["stability"]], levels = trajectory_measures), ymin = -5, ymax = 5)
reference_labels <- tibble(
  measure = factor(trajectory_measures[c("precision", "stability")], levels = trajectory_measures),
  country_f = factor("Ethiopia", levels = country_labels),
  y = c(10.6, 5.6),
  label = c("10 pp", "±5 pp")
)

# create trajectory plot
figure_21 <-
  ggplot() +
  geom_rect(data = stability_band, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), fill = ink[["muted"]], alpha = 0.14, inherit.aes = FALSE) +
  geom_hline(data = precision_line, aes(yintercept = y), colour = ink[["muted"]], linewidth = 0.4, linetype = "22") +
  geom_hline(data = stability_zero, aes(yintercept = y), colour = ink[["muted"]], linewidth = 0.3) +
  geom_line(data = lines, aes(n, value, colour = method, group = method), linewidth = 0.6) +
  geom_linerange(data = observed_markers, aes(y = slot, xmin = stopping_n_q10, xmax = stopping_n_q90, colour = method), linewidth = 0.9, show.legend = FALSE) +
  geom_point(data = observed_markers, aes(stopping_n, slot, colour = method, fill = method), shape = 21, size = 2.8, stroke = 0.9, show.legend = FALSE) +
  geom_point(data = modelled_markers, aes(n_final, slot, colour = method), shape = 124, size = 3.2, show.legend = FALSE) +
  geom_point(data = modelled_markers, aes(stopping_n, slot, colour = method), shape = 22, fill = "white", size = 2.8, stroke = 0.9, show.legend = FALSE) +
  geom_text(
    data = modelled_markers, aes(x = n_final, y = slot, label = "modelled"), hjust = 1.15, vjust = -0.9, size = 2.4, fontface = "italic",
    colour = ink[["muted"]], family = report_font
  ) +
  geom_text(data = reference_labels, aes(x = Inf, y = y, label = label), hjust = 1.1, vjust = 0, size = 2.4, colour = ink[["muted"]], family = report_font) +
  facet_grid2(measure ~ country_f, scales = "free_y", switch = "y", labeller = labeller(measure = strip_labels)) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 32), breaks = c(0, 10, 20, 30), oob = scales::oob_keep, expand = c(0, 0)),
    scale_y_continuous(limits = c(-15, 10), breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_keep, expand = c(0, 0)),
    scale_y_continuous(limits = c(0.4, 4.6), breaks = 1:4, labels = unname(method_labels[rev(method_order)]), expand = c(0, 0))
  )) +
  force_panelsizes(rows = unit(c(3.4, 3.0, 2.0), "cm")) +
  scale_fill_manual(values = method_colours, guide = "none") +
  method_scale(aesthetic = "colour") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.9, shape = NA, fill = NA))) +
  scale_x_continuous(limits = c(0, 1120), breaks = seq(0, 1000, 250), labels = scales::label_comma(), expand = c(0, 0), oob = scales::oob_keep) +
  labs(x = "Directly reporting respondents enrolled", y = NULL) +
  theme_report() +
  theme(
    strip.text.y.left = element_markdown(angle = 0, hjust = 1, vjust = 0.5, lineheight = 1.15, margin = margin(0, 4, 0, 2)),
    strip.placement = "outside",
    panel.grid.major.x = element_blank(),
    panel.spacing.y = unit(1, "lines")
  )

save_figure(figure_21, "output/figures/figure_21.png", height_cm = 12)

## Annex Figure L1 ----
# Respondents required (logarithmic scale) and total cost by target margin of
# error, by country, with each method labelled at the end of its line and the
# report's 5-point criterion dashed.
# panel rows and the criterion's position on x
planning_measures <- c(respondents = "Respondents required", cost = "Total cost (US$, thousands)")
criterion_x <- match(5, target_moe)

# two measures stacked, x as the position of the margin
planning_lines <-
  bind_rows(
    planning |> mutate(measure = planning_measures[["respondents"]], value = n_required),
    planning |> mutate(measure = planning_measures[["cost"]], value = total / 1000)
  ) |>
  transmute(
    country, country_f = factor(country_labels[country], levels = country_labels), method = factor(method, levels = method_order),
    xpos = match(moe, target_moe), measure = factor(measure, levels = planning_measures), value
  ) |>
  arrange(country_f, measure, method, xpos)

# move labels of one panel apart by at least `gap`, keeping their order
spread_labels <- function(y, gap) {
  o <- order(y)
  z <- y[o]
  for (i in seq_along(z)[-1]) z[i] <- max(z[i], z[i - 1] + gap)
  y[o] <- z
  y
}

# method labels at the right end of each line
planning_ends <-
  planning_lines |>
  mutate(range_log = diff(range(log10(value))), range = diff(range(value)), .by = c(country_f, measure)) |>
  filter(xpos == length(target_moe)) |>
  mutate(
    label = method_labels[as.character(method)],
    label_y = if (first(measure) == planning_measures[["respondents"]]) 10^spread_labels(log10(value), 0.075 * first(range_log)) else spread_labels(value, 0.075 * first(range)),
    .by = c(country_f, measure)
  )

# criterion note above the first panel
criterion_note <- tibble(
  country_f = factor(country_labels[[1]], levels = country_labels), measure = factor(planning_measures[["respondents"]], levels = planning_measures),
  x = criterion_x - 0.08,
  y = 1.04 * max(filter(planning_lines, measure == planning_measures[["respondents"]], country_f == country_labels[[1]])$value)
)

# create planning plot
figure_L1 <-
  ggplot(planning_lines, aes(xpos, value, colour = I(method_colours[as.character(method)]), group = method)) +
  geom_vline(xintercept = criterion_x, colour = accent_colour, linewidth = 0.4, linetype = "22") +
  geom_text(
    data = criterion_note, aes(x = x, y = y, label = "±5 pp: the report’s\nprecision criterion"), inherit.aes = FALSE,
    hjust = 1, vjust = 1, family = report_font, size = 2.2, colour = accent_colour, lineheight = 0.95
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(fill = I(method_colours[as.character(method)])), shape = 21, colour = "white", size = 1.9, stroke = 0.45) +
  geom_text(data = planning_ends, aes(x = length(target_moe) + 0.18, y = label_y, label = label), hjust = 0, family = report_font, fontface = "bold", size = 2.4) +
  facet_grid2(measure ~ country_f, scales = "free_y", switch = "y") +
  facetted_pos_scales(y = list(
    measure == planning_measures[["respondents"]] ~ scale_y_log10(breaks = c(50, 100, 200, 500, 1000, 2000, 5000), labels = scales::label_comma()),
    measure != planning_measures[["respondents"]] ~ scale_y_continuous(limits = c(10, NA))
  )) +
  scale_x_continuous(breaks = seq_along(target_moe), labels = paste0("±", target_moe), expand = expansion(add = c(0.25, 0.75))) +
  labs(x = "Target margin of error (percentage points)", y = NULL) +
  theme_report() +
  theme(
    strip.placement = "outside", strip.text.y.left = element_text(angle = 90, face = "bold"),
    panel.grid.major.x = element_blank(), panel.spacing.y = unit(1.1, "lines")
  )

save_figure(figure_L1, "output/figures/figure_L1.png", height_cm = 11)
