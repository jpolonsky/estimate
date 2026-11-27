# 05_validity.R
# Phase 1 validity of CBE, RDS and VBS against the household survey
# Reproduces: Table 9, Figures 10 and 12, Annex Table G11

# Setup ----

pacman::p_load(RDS, ggrepel, ggtext, patchwork, tidyverse)

source("code/functions.R")

# Data ----

estimates <- read_rds("data/clean/estimates.rds")
frames <- read_rds("data/clean/frames.rds")

# Analysis ----

## Comparison with the household survey ----
# Each method's difference from the household survey on the three primary
# antigens and zero-dose, with Newcombe's interval for the difference.
# count children with an observed outcome and with the outcome
counts <- event_counts(frames)

# compare each method with the household survey
performance <- phase1_performance(estimates, counts)

## Summary by country and pooled ----
# RMSE, noise floor, interval width and the verdict counts per method.
# summarise per country and over both countries
summary <- phase1_summary(performance)

## Paired bootstrap of RMSE differences ----
# The difference in RMSE against the household survey between two methods,
# with a 95% percentile interval from 10,000 bootstrap replicates of each
# Phase 1 sample: clusters for HH and CBE, caregivers for VBS and recruitment
# chains for RDS, each with its own seed. Each sample's replicates are
# rescaled to the standard error of its primary estimate, and both methods of
# a pair are compared with the same household-survey replicates.
# children, outcomes, weights and resampling units of each Phase 1 sample
boot_outcomes <- c(opv_3 = "OPV3", penta_3 = "Penta3", measles_1 = "MCV1", zero_dose = "Zero-dose")
boot_samples <- str_subset(study_samples, "phase1")
boot_units <- map(set_names(boot_samples), \(sample) {
  cohort <- direct_cohort(frames[[sample]], sample)
  method <- str_extract(sample, "[a-z]+$")
  y <- sapply(boot_outcomes, \(antigen) cohort[[antigen_outcomes[[antigen]]]])
  if (method == "rds") {
    nodes <- rds_network(cohort, antigen_outcomes)
    design <- as.rds.data.frame(nodes, id = "id_participant_calc", recruiter.id = "id_referring", network.size = "degree", population.size = population_size[[sample]])
    w <- sapply(boot_outcomes, \(antigen) {
      field <- antigen_outcomes[[antigen]]
      seed <- derived_seed(sample_seed(sample), paste0("rds_outcome/", if (field == "zero_dose") "penta_1p_clean" else field))
      rds_child_weights(design, cohort, field, population_size[[sample]], seed)
    })
    roots <- chain_root(nodes)
    unit <- roots[match(as.character(cohort$id_participant_calc), as.character(nodes$id_participant_calc))]
    units <- unique(roots)
  } else {
    w <- matrix(if (method == "vbs") vbs_weight(cohort, sample) else 1, nrow(cohort), length(boot_outcomes))
    unit <- switch(method, hh = cohort$cluster_id, cbe = cohort$cluster, vbs = cohort$record_ix)
    units <- unique(unit)
  }
  usable <- !is.na(y) & !is.na(w)
  sums <- \(m) rowsum(ifelse(usable, m, 0), group = factor(unit, levels = units))
  list(num = sums(y * w), den = sums(w), estimate = colSums(sums(y * w)) / colSums(sums(w)))
})

# draw the replicates: resample units with replacement and recompute each outcome
boot_draw <- function(num, den, seed, reps = 10000) {
  withr::with_seed(seed, {
    out <- matrix(NA_real_, reps, ncol(num), dimnames = list(NULL, colnames(num)))
    for (b in seq_len(reps)) {
      counts <- tabulate(sample.int(nrow(num), nrow(num), replace = TRUE), nbins = nrow(num))
      out[b, ] <- ifelse(drop(counts %*% den) > 0, drop(counts %*% num) / drop(counts %*% den), NA_real_)
    }
    out
  })
}

# rescale each outcome's replicates to the primary standard error, kept within 0 and 1
boot_replicates <- imap(boot_units, \(sums, sample) {
  draws <- boot_draw(sums$num, sums$den, derived_seed(study_seed, paste0("rmse_difference_bootstrap/", str_replace(sample, "phase1", "phase_1"))))
  primary <- filter(estimates, .data$sample == .env$sample)
  se <- primary$standard_error[match(boot_outcomes, primary$antigen)]
  sapply(seq_along(boot_outcomes), \(j) {
    x <- draws[, j]
    pmin(pmax(sums$estimate[[j]] + (x - mean(x, na.rm = TRUE)) * se[[j]] / sd(x, na.rm = TRUE), 0), 1)
  }) |>
    `colnames<-`(boot_outcomes)
})

# RMSE difference and its interval for each pair, country set and outcome group
boot_pairs <- tibble(method_a = c("cbe", "cbe", "rds"), method_b = c("rds", "vbs", "vbs"))
rmse_bootstrap <-
  expand_grid(scope = c("primary_antigens", "zero_dose"), boot_pairs, countries = list(eth = "eth", som = "som", pooled = c("eth", "som"))) |>
  mutate(country = names(countries)) |>
  pmap(\(scope, method_a, method_b, countries, country) {
    antigens <- if (scope == "zero_dose") "Zero-dose" else c("Penta3", "OPV3", "MCV1")
    cells <- filter(performance, country %in% countries, antigen %in% antigens)
    replicate_rmse <- \(method) {
      d <- map(countries, \(k) map(antigens, \(a) 100 * (boot_replicates[[str_c(k, "_phase1_", method)]][, a] - boot_replicates[[str_c(k, "_phase1_hh")]][, a]))) |>
        flatten() |>
        reduce(cbind)
      sqrt(rowMeans(matrix(d, nrow = nrow(boot_replicates[[1]]))^2))
    }
    point <- \(method) sqrt(mean(filter(cells, method == str_to_upper(.env$method))$difference_pp^2))
    diffs <- replicate_rmse(method_a) - replicate_rmse(method_b)
    tibble(
      scope, method_a, method_b, country,
      difference = point(method_a) - point(method_b),
      lower = quantile(diffs[is.finite(diffs)], 0.025, names = FALSE),
      upper = quantile(diffs[is.finite(diffs)], 0.975, names = FALSE),
      boundary = any(cells$boundary[cells$method %in% str_to_upper(c(method_a, method_b))])
    )
  }) |>
  list_rbind()

# Tables ----

## Table 9 ----
# Per country and pooled: the summary over the three primary antigens, with the
# zero-dose difference and verdict alongside. A dagger marks a verdict resting
# on a count-based interval.
# primary-antigen summary and the zero-dose columns
primary <- summary |> filter(scope == "primary_antigens")
zero_dose_pooled <- summary |> filter(scope == "zero_dose", country == "pooled") |> select(method, zero_cells = n_cells, zero_not_different = n_not_different, zero_boundary = n_not_different_boundary)
zero_dose_country <- performance |> filter(antigen == "Zero-dose") |> select(country, method, difference_pp, not_different, boundary)

# format the cells as printed
table_09 <-
  primary |>
  mutate(country = as.character(country)) |>
  left_join(zero_dose_country, by = c("country", "method")) |>
  left_join(zero_dose_pooled, by = "method") |>
  mutate(
    dagger = if_else(n_not_different_boundary > 0, "†", ""),
    zero_difference = if_else(country == "pooled", "–", format_signed(difference_pp)),
    zero_verdict = if_else(
      country == "pooled",
      paste0(zero_not_different, " of ", zero_cells, if_else(zero_boundary > 0, "†", "")),
      paste0(if_else(not_different, "Yes", "No"), if_else(boundary, "†", ""))
    )
  ) |>
  transmute(
    Country = if_else(method == "CBE", c(country_labels, pooled = "Pooled")[country], ""),
    Method = as.character(method),
    `RMSE (pp)` = format_number(rmse_pp),
    `Noise floor (pp)` = format_number(noise_floor_pp),
    `Mean CI width (pp)` = format_number(mean_ci_width_pp),
    `Actionable (antigens with CI ≤ 10 pp, of 3)` = paste(n_ci_width_le_10, "of", n_cells),
    `Not distinguishable from HH (of 3)` = paste0(n_not_different, " of ", n_cells, dagger),
    `Zero-dose: difference from HH (pp)` = zero_difference,
    `Zero-dose: not distinguishable` = zero_verdict
  )

save_table(table_09, "output/tables/table_09.csv")

## Annex Table G11 ----
# Difference in RMSE against the household survey between two methods (pp,
# first minus second, 95% interval), by country and pooled; a dagger marks a
# comparison involving an estimate on a count-based interval.
# one row per outcome group and pair, the country sets as columns
table_G11 <-
  rmse_bootstrap |>
  mutate(
    Outcome = c(primary_antigens = "Primary antigens", zero_dose = "Zero-dose prevalence")[scope],
    `Pair (first minus second)` = str_c(method_labels[method_a], " − ", method_labels[method_b]),
    value = str_c(format_difference(difference, lower, upper), if_else(boundary, "†", "")),
    country = c(eth = "Ethiopia", som = "Somalia", pooled = "Pooled")[country]
  ) |>
  select(Outcome, `Pair (first minus second)`, country, value) |>
  pivot_wider(names_from = country, values_from = value)

save_table(table_G11, "output/tables/table_G11.csv")

# Figures ----

## Figure 10 ----
# Scorecard of accuracy, precision and calibration for every method and
# outcome; calibration is read only where the difference from the household
# survey is beyond the accuracy tolerance.
# signed tile label with an explicit plus sign
tile_signed <- function(x) {
  r <- round(x, 1)
  case_when(r == 0 ~ "0.0", r < 0 ~ paste0("−", sprintf("%.1f", abs(r))), .default = paste0("+", sprintf("%.1f", r)))
}

# verdict and label of each tile
tolerance <- if_else(performance$antigen == "Zero-dose", 5, 10)
beyond <- abs(performance$difference_pp) > tolerance + 1e-9
tiles <-
  bind_rows(
    performance |> mutate(
      measure = "Accuracy", label = tile_signed(difference_pp),
      verdict = if_else(beyond, "fails", "meets")
    ),
    performance |> mutate(
      measure = "Precision", label = sprintf("%.1f", ci_width_pp),
      verdict = case_when(antigen == "Zero-dose" ~ "no_threshold", ci_width_pp <= 10 + 1e-9 ~ "meets", .default = "fails")
    ),
    performance |> mutate(
      measure = "Calibration",
      label = case_when(boundary ~ "count-\nbased", not_different ~ "ns", .default = "diff"),
      verdict = case_when(!beyond ~ "not_read", boundary ~ "not_readable", not_different ~ "beyond_ns", .default = "fails")
    )
  ) |>
  mutate(
    country = factor(country_labels[country], levels = c("Ethiopia", "Somalia")),
    measure = factor(measure, levels = c("Accuracy", "Precision", "Calibration")),
    outcome = antigen,
    row = factor(method, levels = rev(c("CBE", "RDS", "VBS"))),
    verdict = factor(verdict, levels = verdict_levels)
  )

# key rows: the general fills, the calibration fills and the tile text
count_based_beyond <- any(tiles$verdict == "not_readable")
key_rows <- list(
  list(key_item("Standard met", "meets"), key_item("Standard not met", "fails"), key_item("Zero-dose precision: no standard", "no_threshold")),
  list(
    key_item("Beyond tolerance, distinguishable", "fails", "diff"),
    key_item("Beyond tolerance, not distinguishable", "beyond_ns", "ns"),
    key_item("Within tolerance, not read", "not_read", "diff / ns / count-based")
  )
)
if (count_based_beyond) key_rows <- c(key_rows, list(list(key_item("Beyond tolerance, count-based interval", "not_readable", "count-based"))))
key_rows <- c(key_rows, list(list(key_item("diff, distinguishable from the household survey; ns, not distinguishable; count-based, count-based interval"))))
key <- scorecard_key(key_rows)

# draw the scorecard
verdicts <- c("meets", "fails", "beyond_ns", "not_read", "no_threshold", if (count_based_beyond) "not_readable")
figure_10 <- scorecard_plot(tiles, verdicts, key)
save_figure(figure_10, "output/figures/figure_10.png", height_cm = 6.75 + key$height_cm)

## Figure 12 ----
# Each method's difference from the household survey in Ethiopia against that
# in Somalia, for the three primary antigens, within a +/-10 pp square.
# pair the Ethiopia and Somalia differences of each method and antigen
scatter <-
  performance |>
  filter(antigen != "Zero-dose") |>
  mutate(antigen = droplevels(antigen)) |>
  select(country, method, antigen, difference_pp) |>
  pivot_wider(names_from = country, values_from = difference_pp) |>
  arrange(method, antigen) |>
  mutate(label = paste(method, antigen), method = str_to_lower(method))

# axis limits reach past the +/-10 pp square; grid redrawn in white over the shaded quadrants
limit <- max(ceiling(max(abs(c(scatter$eth, scatter$som))) / 2.5) * 2.5, 12.5)
breaks <- seq(-floor(limit / 5) * 5, floor(limit / 5) * 5, by = 5)
positive <- breaks[breaks > 0]
negative <- breaks[breaks < 0]
white_grid <- bind_rows(
  tibble(x = 0, xend = limit, y = positive, yend = positive), tibble(x = positive, xend = positive, y = 0, yend = limit),
  tibble(x = -limit, xend = 0, y = negative, yend = negative), tibble(x = negative, xend = negative, y = -limit, yend = 0)
)
pad <- limit * 0.03
corners <- tibble(
  x = c(limit - pad, -limit + pad, -limit + pad, limit - pad),
  y = c(limit - pad, -limit + pad, limit - pad, -limit + pad),
  label = c("Over in both", "Under in both", "Under in Ethiopia,\nover in Somalia", "Over in Ethiopia,\nunder in Somalia"),
  hjust = c(1, 0, 0, 1), vjust = c(1, 0, 1, 0)
)

# label positions beside each point so none touches an axis or the square
offsets <- tibble(
  label = c("RDS OPV3", "RDS MCV1", "CBE Penta3", "RDS Penta3", "VBS Penta3", "VBS OPV3", "VBS MCV1", "CBE MCV1"),
  dx = c(-0.7, 0.4, 0.7, -1.0, 0.9, 0.9, 0, 0), dy = c(-0.6, 1.1, 0, 0, 0, 0, -1.5, -1.5), hj = c(1, 1, 0, 1, 0, 0, 0.5, 0.5)
)
scatter <-
  scatter |>
  left_join(offsets, by = "label") |>
  mutate(dx = replace_na(dx, 0), dy = replace_na(dy, 1.5), hj = replace_na(hj, 0.5))

# draw the scatter
figure_12 <-
  ggplot(scatter, aes(eth, som)) +
  annotate("rect", xmin = 0, xmax = limit, ymin = 0, ymax = limit, fill = "#F2F2F4") +
  annotate("rect", xmin = -limit, xmax = 0, ymin = -limit, ymax = 0, fill = "#F2F2F4") +
  geom_segment(data = white_grid, aes(x = x, xend = xend, y = y, yend = yend), colour = "white", linewidth = 0.3, inherit.aes = FALSE) +
  annotate("rect", xmin = -10, xmax = 10, ymin = -10, ymax = 10, fill = NA, colour = ink[["muted"]], linewidth = 0.3) +
  annotate(
    "text", x = -10 + pad / 2, y = 10 - pad / 2, label = "±10 pp in both countries", hjust = 0, vjust = 1,
    size = 2.3, colour = ink[["muted"]], family = report_font, fontface = "italic"
  ) +
  geom_text(
    data = corners, aes(x, y, label = label, hjust = hjust, vjust = vjust), inherit.aes = FALSE,
    size = 2.3, colour = ink[["muted"]], family = report_font, fontface = "italic", lineheight = 0.9
  ) +
  geom_hline(yintercept = 0, colour = ink[["secondary"]], linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = ink[["secondary"]], linewidth = 0.4) +
  geom_polygon(aes(group = method, fill = method), alpha = 0.10, colour = NA, show.legend = FALSE) +
  geom_point(aes(fill = method, shape = antigen), colour = "white", size = 2.9, stroke = 0.6) +
  geom_text_repel(
    aes(label = label, hjust = hj), colour = ink[["secondary"]], family = report_font, size = 2.4,
    bg.colour = "white", bg.r = 0.15, nudge_x = scatter$dx, nudge_y = scatter$dy, force = 0, force_pull = 0,
    min.segment.length = Inf, max.overlaps = Inf, seed = 1, show.legend = FALSE
  ) +
  method_scale(c("cbe", "rds", "vbs"), "fill") +
  scale_shape_manual(values = c(Penta3 = 21, OPV3 = 24, MCV1 = 22), name = "Antigen") +
  guides(
    fill = guide_legend(order = 1, override.aes = list(shape = 21, colour = "white", size = 2.9)),
    shape = guide_legend(order = 2, override.aes = list(fill = ink[["secondary"]], colour = "white", size = 2.9))
  ) +
  scale_x_continuous(limits = c(-limit, limit), breaks = breaks, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-limit, limit), breaks = breaks, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    x = "Ethiopia: difference from the household survey (percentage points)",
    y = "Somalia: difference from the household survey (percentage points)"
  ) +
  theme_report()

save_figure(figure_12, "output/figures/figure_12.png", height_cm = 12.5)
