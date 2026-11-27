# 06_phase2.R
# Phase 2: agreement between CBE, RDS and VBS
# Reproduces: Tables 10 and 12, Figures 13-15, Annex Tables J1 and J2

# Setup ----

pacman::p_load(RDS, ggtext, patchwork, tidyverse)

source("code/functions.R")

# Data ----

estimates <- read_rds("data/clean/estimates.rds")
frames <- read_rds("data/clean/frames.rds")
tracking <- read_rds("data/clean/rds_tracking.rds")

# Analysis ----

## Inter-method range ----
# Highest minus lowest of the three method estimates per country and outcome.
# calculate the range per country and outcome
ranges <- intermethod_range(estimates)

## Pairwise comparisons ----
# Each pair of methods within a country, with Newcombe's interval for the
# difference (count-based for an arm with fewer than five children with or
# without the outcome).
# count children with an observed outcome and with the outcome
counts <- event_counts(frames)

# compare every pair of methods
pairwise <- pairwise_tests(estimates, counts)

## Precision ----
# Width of each 95% CI, and the primary antigens whose interval is at most
# 10 pp wide.
# calculate interval widths
phase2 <-
  estimates |>
  filter(phase == 2, antigen %in% primary_outcomes) |>
  mutate(antigen = factor(antigen, levels = primary_outcomes), method = as.character(method), ci_width_pp = 100 * (ci_upper - ci_lower))

# count the actionable antigens per country and method
actionable <-
  phase2 |>
  filter(antigen != "Zero-dose") |>
  group_by(country, method) |>
  summarise(n_narrow = sum(ci_width_pp <= 10 + 1e-9), .groups = "drop")

## Predicted against observed differences ----
# In Somalia, the Phase 2 difference between CBE or VBS and RDS set against the
# difference predicted from Phase 1 (each method's difference from the
# household survey minus that of RDS).
# predict the gap from the Phase 1 differences
predicted <-
  phase1_performance(estimates, counts) |>
  filter(country == "som") |>
  select(method, antigen, difference_pp) |>
  pivot_wider(names_from = method, values_from = difference_pp) |>
  transmute(antigen, CBE = CBE - RDS, VBS = VBS - RDS) |>
  pivot_longer(c(CBE, VBS), names_to = "method", values_to = "predicted_pp")

# observed gap: the pairwise difference, signed as method minus RDS
observed <-
  pairwise |>
  filter(country == "som", method_a == "RDS" | method_b == "RDS") |>
  mutate(
    method = if_else(method_a == "RDS", method_b, method_a),
    sign = if_else(method_a == "RDS", -1, 1),
    observed_pp = sign * diff_pp,
    observed_lower_pp = pmin(sign * diff_lower_pp, sign * diff_upper_pp),
    observed_upper_pp = pmax(sign * diff_lower_pp, sign * diff_upper_pp)
  ) |>
  select(method, antigen, observed_pp, observed_lower_pp, observed_upper_pp)

# compare direction and size
transfer <-
  predicted |>
  left_join(observed, by = c("method", "antigen")) |>
  mutate(
    same_direction = sign(predicted_pp) == sign(observed_pp),
    within_3pp = abs(predicted_pp - observed_pp) <= 3
  ) |>
  arrange(method, antigen)

## Somalia RDS record sets ----
# The Somalia RDS estimates from four sets of records (the forms of the four
# enumerators retained after the data-quality review or of all six, with every
# record or only those collected to 29 May 2026), each set in turn taking the
# place of the primary RDS estimates in the range and the pairwise comparisons.
record_sets <- tribble(
  ~record_set, ~frame, ~to_29_may,
  "Four enumerators, all records (primary)", "som_phase2_rds", FALSE,
  "Four enumerators, records to 29 May", "som_phase2_rds", TRUE,
  "Six enumerators, all records", "som_phase2_rds_all_enumerators", FALSE,
  "Six enumerators, records to 29 May", "som_phase2_rds_all_enumerators", TRUE
)

# estimate the four outcomes from one set of records, with the children counted
record_set_estimates <- function(frame, to_29_may) {
  cohort <- direct_cohort(frames[[frame]], frame)
  if (to_29_may) cohort <- filter(cohort, as.Date(date) < as.Date("2026-05-30"))
  population <- population_size[["som_phase2_rds"]]
  nodes <- rds_network(cohort, antigen_outcomes, tracking = tracking)
  design <- as.rds.data.frame(nodes, id = "id_participant_calc", recruiter.id = "id_referring", network.size = "degree", population.size = population)
  map(primary_outcomes, \(antigen) {
    field <- antigen_outcomes[[antigen]]
    seed_field <- if (field == "zero_dose") "penta_1p_clean" else field
    seed <- derived_seed(sample_seed("som_phase2_rds"), paste0("rds_outcome/", seed_field))
    rds_proportion(design, cohort, field, population, seed) |>
      mutate(antigen = antigen, n_domain = sum(!is.na(cohort[[field]])), n_events = sum(cohort[[field]] == 1, na.rm = TRUE))
  }) |>
    list_rbind()
}

# swap each set into the primary estimates and compare the methods again
record_set_results <- pmap(record_sets, \(record_set, frame, to_29_may) {
  version <- record_set_estimates(frame, to_29_may)
  swapped_estimates <- rows_update(
    estimates,
    version |> transmute(sample = "som_phase2_rds", antigen = factor(antigen, levels = levels(estimates$antigen)), estimate, standard_error, ci_lower, ci_upper),
    by = c("sample", "antigen")
  )
  swapped_counts <- rows_update(counts, version |> transmute(sample = "som_phase2_rds", antigen, n_domain, n_events), by = c("sample", "antigen"))
  tests <- pairwise_tests(swapped_estimates, swapped_counts) |> filter(country == "som")
  version |>
    mutate(antigen = factor(antigen, levels = primary_outcomes), record_set = record_set) |>
    left_join(intermethod_range(swapped_estimates) |> filter(country == "som") |> select(antigen, range_pp), by = "antigen") |>
    left_join(tests |> filter(method_a == "CBE", method_b == "RDS") |> select(antigen, diff_pp, diff_lower_pp, diff_upper_pp, cbe_rds = distinguishable), by = "antigen") |>
    left_join(tests |> filter(method_a == "RDS", method_b == "VBS") |> select(antigen, rds_vbs = distinguishable), by = "antigen")
}) |>
  list_rbind()

# Tables ----

## Table 10 ----
# Percentage (95% CI) by country and outcome, one column per method, with the
# inter-method range and, per method, the primary antigens with a CI at most
# 10 pp wide. A dagger marks an estimate that entered a comparison with its
# count-based interval.
# estimate with CI, marked where the arm used a count-based interval
count_based <-
  bind_rows(
    pairwise |> filter(boundary_a) |> select(country, antigen, method = method_a),
    pairwise |> filter(boundary_b) |> select(country, antigen, method = method_b)
  ) |>
  distinct() |>
  mutate(dagger = "†")

# one row per country and outcome, methods across columns
estimate_rows <-
  phase2 |>
  left_join(count_based, by = c("country", "antigen", "method")) |>
  mutate(value = paste0(format_ci(estimate, ci_lower, ci_upper), replace_na(dagger, ""))) |>
  select(country, antigen, method, value) |>
  pivot_wider(names_from = method, values_from = value) |>
  left_join(ranges |> select(country, antigen, range_pp), by = c("country", "antigen")) |>
  transmute(country, Outcome = if_else(antigen == "Zero-dose", "Zero-dose prevalence", as.character(antigen)), CBE, RDS, VBS,
            `Inter-method range (pp)` = format_number(range_pp), order = as.integer(antigen))

# actionable row per country
actionable_rows <-
  actionable |>
  mutate(value = sprintf("%s (%d of 3)", if_else(n_narrow == 3, "Yes", "No"), n_narrow)) |>
  select(country, method, value) |>
  pivot_wider(names_from = method, values_from = value) |>
  mutate(Outcome = "Actionable", `Inter-method range (pp)` = "", order = 5L)

# stack the rows of each country
table_10 <-
  bind_rows(estimate_rows, actionable_rows) |>
  arrange(country, order) |>
  mutate(Country = if_else(order == 1, country_labels[country], "")) |>
  select(Country, Outcome, CBE, RDS, VBS, `Inter-method range (pp)`)

save_table(table_10, "output/tables/table_10.csv")

## Table 12 ----
# Somalia: the Phase 2 difference of CBE and VBS from RDS against the
# difference predicted from Phase 1.
# format the predicted and observed differences
table_12 <-
  transfer |>
  transmute(
    Pair = paste(method, "− RDS"),
    Outcome = if_else(antigen == "Zero-dose", "Zero-dose prevalence", as.character(antigen)),
    `Predicted from Phase 1` = format_signed(predicted_pp),
    `Observed in Phase 2 (95% CI)` = format_difference(observed_pp, observed_lower_pp, observed_upper_pp),
    `Same direction` = if_else(same_direction, "Yes", "No"),
    `Within 3 pp` = if_else(within_3pp, "Yes", "No")
  )

save_table(table_12, "output/tables/table_12.csv")

## Annex Table J1 ----
# Every pairwise comparison of the methods in Phase 2: the two estimates, the
# difference with its Newcombe interval, the verdict, and the arm that entered
# with a count-based interval.
# estimates of the two methods in each pair
pair_estimates <-
  phase2 |>
  transmute(country, antigen = as.character(antigen), method, estimate)

# one row per country, outcome and pair
table_J01 <-
  pairwise |>
  mutate(antigen = as.character(antigen)) |>
  left_join(pair_estimates |> rename(method_a = method, estimate_a = estimate), by = c("country", "antigen", "method_a")) |>
  left_join(pair_estimates |> rename(method_b = method, estimate_b = estimate), by = c("country", "antigen", "method_b")) |>
  transmute(
    Country = country_labels[country],
    Outcome = outcome_label(antigen),
    Pair = paste0(method_a, "–", method_b),
    `First method (%)` = format_number(100 * estimate_a),
    `Second method (%)` = format_number(100 * estimate_b),
    `Difference, first minus second, pp (95% CI)` = format_difference_range(diff_pp, diff_lower_pp, diff_upper_pp),
    Verdict = verdict_label(distinguishable),
    `Count-based interval` = case_when(
      boundary_a & boundary_b ~ "Both arms", boundary_a ~ paste(method_a, "arm"), boundary_b ~ paste(method_b, "arm"), .default = "No"
    )
  )

save_table(table_J01, "output/tables/table_J01.csv")

## Annex Table J2 ----
# Somalia: the RDS estimate, the inter-method range and the CBE-RDS and RDS-VBS
# verdicts under each of the four sets of RDS records.
# one row per record set and outcome
table_J02 <-
  record_set_results |>
  transmute(
    `Record set` = record_set,
    Outcome = outcome_label(antigen),
    `RDS estimate (%)` = format_number(100 * estimate),
    `RDS children (n)` = format_count(n_domain),
    `Inter-method range (pp)` = format_number(range_pp),
    `CBE minus RDS, pp (95% CI)` = format_difference_range(diff_pp, diff_lower_pp, diff_upper_pp),
    `CBE–RDS verdict` = verdict_label(cbe_rds),
    `RDS–VBS verdict` = verdict_label(rds_vbs)
  )

save_table(table_J02, "output/tables/table_J02.csv")

# Figures ----

## Figure 13 ----
# Scorecard of agreement and precision for every method and outcome, with the
# inter-method range as the top row of each country.
# agreement: a method meets when it is not distinguishable from either other
# method, and fails (with the number of its two comparisons that differ) when
# it is distinguishable from one or both
agreement <-
  phase2 |>
  select(country, method, antigen) |>
  mutate(
    n_distinguishable = pmap_int(list(country, method, antigen), \(c, m, a) {
      p <- pairwise |> filter(country == c, antigen == a, method_a == m | method_b == m)
      sum(p$distinguishable)
    }),
    any_boundary = pmap_lgl(list(country, method, antigen), \(c, m, a) {
      p <- pairwise |> filter(country == c, antigen == a, method_a == m | method_b == m)
      any(p$boundary)
    }),
    measure = "Agreement",
    verdict = case_when(n_distinguishable > 0 ~ "fails", any_boundary ~ "not_readable", .default = "meets"),
    label = case_when(n_distinguishable > 0 ~ sprintf("%d of 2", n_distinguishable), any_boundary ~ "count-\nbased", .default = "")
  )

# precision tiles, and the inter-method range as a top row
precision <-
  phase2 |>
  mutate(
    measure = "Precision", label = sprintf("%.1f", ci_width_pp),
    verdict = case_when(antigen == "Zero-dose" ~ "no_threshold", ci_width_pp <= 10 + 1e-9 ~ "meets", .default = "fails")
  )
range_tiles <-
  ranges |>
  mutate(method = "range", measure = "Agreement", label = sprintf("%.1f", range_pp), verdict = if_else(range_pp < 10, "meets", "fails"))

# stack the tiles and order the rows and columns
tiles <-
  bind_rows(agreement, precision, range_tiles) |>
  mutate(
    country = factor(country_labels[country], levels = c("Ethiopia", "Somalia")),
    measure = factor(measure, levels = c("Agreement", "Precision")),
    outcome = factor(antigen, levels = primary_outcomes),
    row = factor(method, levels = rev(c("range", "CBE", "RDS", "VBS")), labels = rev(c("Inter-method\nrange", "CBE", "RDS", "VBS"))),
    verdict = factor(verdict, levels = verdict_levels)
  )

# key rows: the general fills, the agreement tiles and the range row
key <- scorecard_key(list(
  list(key_item("Standard met", "meets"), key_item("Standard not met", "fails"), key_item("Zero-dose precision: no standard", "no_threshold")),
  list(
    key_item("Not distinguishable from either other method", "meets"),
    key_item("Distinguishable from one or both other methods", "fails", "1 of 2, 2 of 2")
  ),
  list(key_item("Count-based interval in one comparison", "not_readable", "count-based"), key_item("Top row of each country: inter-method range (pp)"))
))

# draw the scorecard
figure_13 <- scorecard_plot(tiles, c("meets", "fails", "not_readable", "no_threshold"), key)
save_figure(figure_13, "output/figures/figure_13.png", height_cm = 7.25 + key$height_cm)

## Figure 14 ----
# Each method's estimate and 95% CI, with a bracket over the inter-method
# range; primary antigens and zero-dose on separate scales.
# create forest plot
figure_14 <- coverage_forest(estimates, phase = 2, ranges = ranges)
save_figure(figure_14, "output/figures/figure_14.png", height_cm = 11)

## Figure 15 ----
# Pairwise differences between the methods with their Newcombe intervals, by
# country and outcome.
# accent colour for a distinguishable pair, grey otherwise, open point for a count-based arm
grey <- "#9A9AA2"
pairs <-
  pairwise |>
  mutate(
    country = factor(country_labels[country], levels = c("Ethiopia", "Somalia")),
    pair = factor(paste(method_a, "minus", method_b), levels = rev(c("CBE minus RDS", "CBE minus VBS", "RDS minus VBS"))),
    clear = distinguishable & !boundary,
    line_colour = if_else(clear, accent_colour, grey),
    fill_colour = if_else(boundary, "white", line_colour)
  )
status_labels <- c(clear = "Interval excludes zero", none = "Interval includes zero", count_based = "Count-based arm (not interpreted)")

# draw the differences with their intervals
figure_15 <-
  ggplot(pairs, aes(diff_pp, pair)) +
  geom_vline(xintercept = 0, colour = ink[["secondary"]], linewidth = 0.4) +
  geom_linerange(aes(xmin = diff_lower_pp, xmax = diff_upper_pp, colour = line_colour), linewidth = 0.55) +
  geom_point(aes(fill = fill_colour, colour = line_colour), shape = 21, size = 2.2, stroke = 0.7) +
  scale_fill_identity() +
  scale_colour_identity() +
  geom_point(
    data = tibble(status = factor(names(status_labels), levels = names(status_labels))), aes(x = Inf, y = Inf, shape = status),
    inherit.aes = FALSE, alpha = 0, na.rm = TRUE
  ) +
  scale_shape_manual(values = c(clear = 21, none = 21, count_based = 21), labels = status_labels, name = NULL) +
  guides(shape = guide_legend(override.aes = list(
    alpha = 1, size = 2.2, stroke = 0.7, colour = c(accent_colour, grey, grey), fill = c(accent_colour, grey, "white")
  ))) +
  facet_grid(antigen ~ country, switch = "y") +
  labs(x = "Difference between methods (percentage points)", y = NULL) +
  theme_report() +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1))

save_figure(figure_15, "output/figures/figure_15.png", height_cm = 13)
