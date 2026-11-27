# 12_field_effort.R
# Field effort: team hours, coordinator time and the length of the interview
# and its sections
# Reproduces: Figures 19 and 20, Annex Tables M1-M4, Annex Figure M1

# Setup ----

pacman::p_load(tidyverse)

source("code/functions.R")

# Data ----

frames <- read_rds("data/clean/frames.rds")
hours <- read_rds("data/clean/hours.rds")
tallies <- read_rds("data/clean/tallies.rds")
record_counts <- read_rds("data/clean/kobo_record_counts.rds")

# Analysis ----

## Team hours ----
# Total team hours including travel, split into enumerator time on site,
# enumerator travel and supervisor time. CBE recorded no supervisor time.
# sum the hours of each sample by component
effort <- team_hours(hours)

## Core interview duration ----
# Median and interquartile range of the minutes from the end of consent to the
# start of the network module, among directly reporting respondents with a
# valid duration. CBE was recorded on paper, without timestamps.
# samples with interview timestamps (CBE was on paper)
timestamped <- c(
  "eth_phase1_hh", "eth_phase1_rds", "eth_phase1_vbs", "som_phase1_hh", "som_phase1_rds", "som_phase1_vbs",
  "eth_phase2_rds", "eth_phase2_vbs", "som_phase2_rds", "som_phase2_vbs"
)

# median and quartiles per sample
duration <-
  timestamped |>
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
  list_rbind(names_to = "sample") |>
  left_join(sample_parts(timestamped), by = "sample")

## Field staff ----
# Enumerators and supervisors tracked, and enumerator person-days. CBE sheets
# record no supervisor.
# count staff and person-days per sample
staff <- staff_counts(hours)

## Coordinator hours ----
# Phase 1 coordinators' weekly hours by activity. Hours logged for both
# countries are split equally between them, and hours logged for no single
# method equally across the four; AAU's data-analysis hours (qualitative
# analysis) are left out.
# apportion each timesheet row to country and method, then sum by activity
coordinator_hours <-
  hours |>
  filter(source == "coordinator timesheet", role == "coordinator", !(organisation == "AAU" & activity == "data analysis")) |>
  mutate(share = if_else(country == "both", 0.5, 1) * if_else(method == "cross-cutting", 0.25, 1)) |>
  left_join(tibble(country = "both", to_country = c("eth", "som")), by = "country", relationship = "many-to-many") |>
  left_join(tibble(method = "cross-cutting", to_method = method_order), by = "method", relationship = "many-to-many") |>
  mutate(country = coalesce(to_country, country), method = coalesce(to_method, method)) |>
  summarise(hours = sum(hours * share), .by = c(country, method, activity))

## Team hours per respondent and per contact ----
# Total team hours including travel per directly reporting respondent, and in
# minutes per contact: households approached (HH and CBE), persons screened
# (RDS) and travellers screened (VBS).
# directly reporting respondents and contacts per sample
contacts <-
  sample_contacts(frames, tallies, record_counts) |>
  mutate(n_respondents = map_int(sample, \(s) nrow(analytic_respondents(frames[[s]], s)))) |>
  left_join(sample_parts(study_samples), by = "sample")

## Interview sections ----
# Minutes in the background questions (end of consent to the end of the
# socio-demographic module), the direct report (to the start of the network
# module) and the network module (to the end of the form), and their sums; a
# section longer than 30 minutes (a paused form) or negative is left out, with
# any sum it feeds.
# section minutes of a set of respondents
section_minutes <- function(respondents) {
  minutes <- function(from, to) {
    m <- as.numeric(difftime(to, from, units = "mins"))
    if_else(m < 0 | m > 30, NA_real_, m)
  }
  tibble(
    socdem = minutes(respondents$t_end_consent, respondents$t_end_socdem),
    direct_module = minutes(respondents$t_end_socdem, respondents$t_start_nsum),
    network_module = minutes(respondents$t_start_nsum, respondents$t_end_survey)
  ) |>
    mutate(total_core = socdem + direct_module, total_full = total_core + network_module)
}

# median and quartiles of each section among directly reporting respondents
sections <-
  timestamped |>
  set_names() |>
  map(\(sample) {
    section_minutes(analytic_respondents(frames[[sample]], sample)) |>
      pivot_longer(everything(), names_to = "section", values_to = "minutes") |>
      summarise(
        median = median(minutes, na.rm = TRUE),
        p25 = quantile(minutes, 0.25, na.rm = TRUE, names = FALSE),
        p75 = quantile(minutes, 0.75, na.rm = TRUE, names = FALSE),
        .by = section
      )
  }) |>
  list_rbind(names_to = "sample") |>
  left_join(sample_parts(timestamped), by = "sample")

## Network-module time ----
# The network module's share of interview time and of total team hours, over
# every interview it was asked in (VBS caregivers outside the 12-23 month
# window included). Respondent-hours are the valid durations plus the median
# for each interview without one.
# respondent-hours of one section
section_hours <- function(minutes) (sum(minutes, na.rm = TRUE) + sum(is.na(minutes)) * median(minutes, na.rm = TRUE)) / 60

# network-module time of each timestamped sample
network_time <-
  timestamped |>
  set_names() |>
  map(\(sample) {
    respondents <- distinct(frames[[sample]], record_ix, .keep_all = TRUE)
    enrolled <- filter(respondents, consent == "Enrolled")
    interviewed <- if (str_detect(sample, "vbs")) filter(enrolled, participant_status == "Enrolled with Eligible Child") else enrolled
    m <- section_minutes(interviewed)
    tibble(
      n_interviews = nrow(interviewed),
      n_enrolled = nrow(enrolled),
      n_timed = sum(!is.na(m$network_module)),
      median_network = median(m$network_module, na.rm = TRUE),
      median_full = median(m$total_full, na.rm = TRUE),
      network_hours = section_hours(m$network_module),
      core_hours = section_hours(m$total_core)
    )
  }) |>
  list_rbind(names_to = "sample") |>
  left_join(sample_parts(timestamped), by = "sample") |>
  left_join(select(effort, phase, country, method, total_hours), by = c("phase", "country", "method"))

# Tables ----

## Annex Table M1 ----
# Staff tracked, person-days and team hours by sample. A double dagger marks
# supervisor hours from the Somalia weekly trackers, which carry no travel.
# one row per sample in report order
table_M01 <-
  staff |>
  left_join(effort, by = c("phase", "country", "method")) |>
  mutate(method = factor(method, levels = method_order)) |>
  arrange(phase, country, method) |>
  transmute(
    `Phase, country` = phase_country(phase, country),
    Method = method_labels[as.character(method)],
    `Enumerators tracked` = format_count(n_enumerators),
    `Supervisors tracked` = if_else(method == "cbe", "not recorded", format_count(n_supervisors)),
    `Person-days` = format_count(n_person_days),
    `On-site hours` = format_number(enumerator_on_site),
    `Enumerator travel hours` = format_number(enumerator_travel),
    `Supervisor hours` = case_when(
      method == "cbe" ~ "not recorded",
      country == "som" & method %in% c("rds", "vbs") ~ str_c(format_number(supervisor), "‡"),
      .default = format_number(supervisor)
    ),
    `Total-team hours` = format_number(total_hours)
  )

save_table(table_M01, "output/tables/table_M01.csv")

## Annex Table M2 ----
# Phase 1 coordinator hours by activity, with preparation as a share of the
# total.
# one row per country and method
table_M02 <-
  coordinator_hours |>
  pivot_wider(names_from = activity, values_from = hours) |>
  mutate(total = preparation + `field work` + `data analysis` + dissemination, method = factor(method, levels = method_order)) |>
  arrange(country, method) |>
  mutate(Country = if_else(row_number() == 1, country_labels[country], ""), .by = country) |>
  transmute(
    Country,
    Method = method_labels[as.character(method)],
    Preparation = format_number(preparation),
    Fieldwork = format_number(`field work`),
    `Data analysis` = format_number(`data analysis`),
    Dissemination = format_number(dissemination),
    Total = format_number(total),
    `Preparation (% of total)` = format_number(100 * preparation / total)
  )

save_table(table_M02, "output/tables/table_M02.csv")

## Annex Table M3 ----
# Team hours per directly reporting respondent and team minutes per contact; a
# dagger marks fewer than 100 directly reporting respondents.
# one row per sample in report order
table_M03 <-
  contacts |>
  left_join(effort, by = c("phase", "country", "method")) |>
  mutate(method = factor(method, levels = method_order)) |>
  arrange(phase, country, method) |>
  transmute(
    `Phase, country` = phase_country(phase, country),
    Method = method_labels[as.character(method)],
    `Total team hours` = format_count(total_hours),
    `Directly reporting respondents` = format_count(n_respondents),
    `Team hours per respondent` = str_c(format_number(total_hours / n_respondents, 2), if_else(n_respondents < 100, "†", "")),
    Contacts = format_count(n_contacts),
    `Team minutes per contact` = format_number(60 * total_hours / n_contacts)
  )

save_table(table_M03, "output/tables/table_M03.csv")

## Annex Table M4 ----
# Interviews timed, median network-module and full-interview minutes, and the
# module's share of interview time and of total team hours; a double dagger
# marks samples whose timed interviews cover part of the enrolled respondents
# only.
# one row per timestamped sample in report order
table_M04 <-
  network_time |>
  mutate(method = factor(method, levels = method_order)) |>
  arrange(phase, country, method) |>
  transmute(
    `Phase, country` = phase_country(phase, country),
    Method = method_labels[as.character(method)],
    `Interviews timed` = format_count(n_timed),
    `Median network module (min)` = format_number(median_network),
    `Median full interview (min)` = format_number(median_full),
    `Share of interview time` = str_c(format_number(100 * network_hours / (core_hours + network_hours)), "%"),
    `Share of total team hours` = str_c(format_number(100 * network_hours / total_hours), "%", if_else(n_interviews < n_enrolled, "‡", ""))
  )

save_table(table_M04, "output/tables/table_M04.csv")

# Figures ----

## Figure 19 ----
# Stacked team hours by method, one panel per country and phase, with the
# total as a label.
# hours by component, the first component darkest
components <- c(enumerator_on_site = "Enumerator on site", enumerator_travel = "Enumerator travel", supervisor = "Supervisor")

# one row per sample and component, with the panel labels
effort_long <-
  effort |>
  pivot_longer(all_of(names(components)), names_to = "component", values_to = "component_hours") |>
  mutate(
    component = factor(components[component], levels = rev(components)),
    method = factor(method, levels = method_order),
    country_f = factor(country_labels[country], levels = country_labels),
    panel_phase = factor(phase_labels[as.character(phase)], levels = phase_labels)
  )

# bar totals
effort_totals <- distinct(effort_long, country_f, panel_phase, method, total_hours)

# components of a bar in the method colour, lighter for each later component,
# with a grey legend key
stack_scale <- function(labels, alphas) {
  list(
    scale_fill_manual(values = method_colours, guide = "none"),
    scale_alpha_manual(
      values = set_names(alphas, labels), breaks = labels, name = NULL,
      guide = guide_legend(override.aes = list(fill = ink[["secondary"]]))
    )
  )
}

# create stacked bar chart
figure_19 <-
  ggplot(effort_long, aes(method, component_hours, fill = method, alpha = component)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.5) +
  geom_text(
    data = effort_totals, aes(method, total_hours, label = sprintf("%.0f", total_hours)), inherit.aes = FALSE,
    vjust = -0.4, size = 2.5, family = report_font, colour = ink[["primary"]]
  ) +
  stack_scale(unname(components), c(1, 0.6, 0.3)) +
  facet_grid(. ~ country_f + panel_phase, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = method_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = scales::label_comma()) +
  labs(x = NULL, y = "Team hours, including travel") +
  theme_report() +
  theme(panel.grid.major.x = element_blank())

save_figure(figure_19, "output/figures/figure_19.png", height_cm = 8.5)

## Figure 20 ----
# Median (point) and interquartile range (bar) of the core interview, with the
# number of valid durations, one panel per phase.
# one row per country and method, in report order from the top
duration_rows <-
  duration |>
  mutate(
    method = factor(method, levels = method_order),
    panel_phase = factor(phase_labels[as.character(phase)], levels = phase_labels)
  ) |>
  arrange(phase, country, method) |>
  mutate(row = str_c(country_labels[country], " ", method_labels[as.character(method)])) |>
  mutate(row = factor(row, levels = rev(unique(row))))

# create interval plot
figure_20 <-
  ggplot(duration_rows, aes(y = row, colour = method)) +
  geom_linerange(aes(xmin = p25, xmax = p75), linewidth = 1) +
  geom_point(aes(median, fill = method), shape = 21, colour = "white", size = 2.6, stroke = 0.6) +
  geom_text(aes(x = p75, label = str_c("n = ", n_valid)), hjust = -0.15, size = 2.3, family = report_font, colour = ink[["secondary"]]) +
  facet_wrap(vars(panel_phase), ncol = 1, scales = "free_y") +
  scale_colour_manual(values = method_colours, guide = "none") +
  scale_fill_manual(values = method_colours, guide = "none") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Core interview minutes (median and interquartile range)", y = NULL) +
  theme_report() +
  theme(strip.text.x = element_text(hjust = 0, face = "bold", size = 9), panel.grid.major.y = element_blank())

save_figure(figure_20, "output/figures/figure_20.png", height_cm = 10)

## Annex Figure M1 ----
# Median (point) and interquartile range (bar) of each interview section, with
# the median of the whole interview without the network module (open circle)
# and the median of the whole interview printed, one row of panels per country
# and phase.
# column labels of the four sections
section_labels <- c(
  socdem = "Background", direct_module = "Direct report\n(own child)", network_module = "Network module\n(other children)",
  total_full = "Whole interview\n(○ without network module)"
)

# rows by country within phase, methods from the top in report order
section_rows <- \(d) mutate(
  d,
  group = factor(str_c(country_labels[country], "\nPhase ", phase), levels = str_c(rep(country_labels, 2), "\nPhase ", rep(1:2, each = 2))),
  method_f = factor(method_labels[method], levels = rev(method_labels))
)
section_bars <-
  sections |>
  filter(section %in% names(section_labels)) |>
  section_rows() |>
  mutate(section_f = factor(section_labels[section], levels = section_labels))
core_points <-
  sections |>
  filter(section == "total_core") |>
  section_rows() |>
  mutate(section_f = factor(section_labels[["total_full"]], levels = section_labels))

# create interval plot by section
figure_M1 <-
  ggplot(section_bars, aes(y = method_f)) +
  geom_segment(aes(x = p25, xend = p75, yend = method_f, colour = I(method_colours[method])), linewidth = 1.6, alpha = 0.35, lineend = "round") +
  geom_point(aes(x = median, fill = I(method_colours[method])), shape = 21, colour = "white", size = 2.1, stroke = 0.5) +
  geom_point(data = core_points, aes(x = median, colour = I(method_colours[method])), shape = 21, fill = "white", size = 1.7, stroke = 0.6) +
  geom_text(
    data = filter(section_bars, section == "total_full"), aes(x = p75 + 0.5, label = sprintf("%.1f", median)), hjust = 0,
    family = report_font, size = 2.2, colour = ink[["secondary"]]
  ) +
  facet_grid(group ~ section_f, switch = "y", scales = "free", space = "free") +
  scale_x_continuous(breaks = seq(0, 20, 5), expand = expansion(add = c(0.3, 1.6))) +
  expand_limits(x = 0) +
  labs(x = "Minutes (median and interquartile range)", y = NULL) +
  theme_report() +
  theme(
    strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", lineheight = 1.1),
    strip.text.x = element_text(lineheight = 1.05, vjust = 0), panel.grid.major.y = element_blank(),
    panel.spacing.x = unit(1, "lines")
  )

save_figure(figure_M1, "output/figures/figure_M1.png", height_cm = 10)
