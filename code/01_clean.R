# 01_clean.R
# Cleaning: one analysis frame per sample, and the supporting tables
# Reproduces: no report item; writes the clean files in data/clean/ that
# 02-12 read

# Setup ----

pacman::p_load(tidyverse)

source("code/functions.R")

dir.create("data/clean", showWarnings = FALSE, recursive = TRUE)

# Data ----

# read one file per form; the second Somalia Phase 2 RDS file keeps the records
# of the two enumerators excluded from the primary analysis
raw <- map(set_names(study_samples), \(sample) read_rds(file.path("data/raw", paste0(sample, ".rds"))))
raw$som_phase2_rds_all_enumerators <- read_rds("data/raw/som_phase2_rds_all_enumerators.rds")

# supporting tables: RDS tracking sheet, cluster weights and places, record
# counts, tally sheets, card photographs, staff hours, budget lines,
# research-team cost, exchange rates and price index
tracking <- read_rds("data/raw/som_phase2_rds_tracking.rds")
cluster_weights <- read_rds("data/raw/cluster_weights.rds")
cluster_places <- read_rds("data/raw/cluster_places.rds")
kobo_record_counts <- read_rds("data/raw/kobo_record_counts.rds")
tallies <- read_rds("data/raw/tallies.rds")
card_photo_status <- read_rds("data/raw/card_photo_status.rds")
time_motion <- read_rds("data/raw/time_motion.rds")
cost_line_items <- read_rds("data/raw/cost_line_items.rds")
cost_shared_items <- read_rds("data/raw/cost_shared_items.rds")
research_team_cost <- read_rds("data/raw/research_team_cost.rds")
fx_cpi <- read_rds("data/raw/fx_cpi.rds")

# Analysis ----

## Respondents ----
# Participation status from the consent answer and the number of eligible
# children, and the respondent's relation to the children.
# classify consent and caregiver relation
clean_respondents <- function(form) {
  main <- form$main
  main$consent <- classify_consent(main, progressed = main$record_ix %in% form$child$record_ix)
  if ("relation_child" %in% names(main)) main$relation_child <- classify_relation(main$relation_child)
  main
}

# clean every respondent table (all forms except CBE)
main <- map(raw[!str_detect(names(raw), "cbe")], clean_respondents)

# VBS Phase 1: interview venue (respondents were met in two towns) and earlier
# ESTIMATE interviews
for (sample in c("eth_phase1_vbs", "som_phase1_vbs")) {
  country <- str_sub(sample, 1, 3)
  main[[sample]]$int_location <- classify_venue(main[[sample]]$int_location, country)
  main[[sample]] <- classify_previous_survey(main[[sample]], country)
}

# VBS: whether the respondent had an eligible child
for (sample in str_subset(names(main), "vbs")) {
  main[[sample]]$participant_status <- classify_participant_status(main[[sample]]$hh_child_num)
}

## Children reported directly ----
# Dose counts and antigen indicators for each child. On the household form,
# children who never received a vaccine skipped the dose questions and are
# coded unvaccinated, and the source of each answer is completed.
# classify dose counts and derive antigen indicators
clean_children <- function(child) {
  child |>
    classify_doses() |>
    derive_antigens()
}

# clean each form's child table
child <- map(raw[names(main)], \(form) clean_children(form$child))

# complete never-vaccinated children and answer sources (HH)
for (sample in c("eth_phase1_hh", "som_phase1_hh")) {
  child[[sample]] <-
    child[[sample]] |>
    complete_never_vaccinated() |>
    code_answer_source(holder_field = "has_card", household = TRUE)
}

# join each respondent to their children
frames <- map(set_names(names(main)), \(sample) {
  left_join(main[[sample]], select(child[[sample]], record_ix, !any_of(names(main[[sample]]))), by = "record_ix")
})

## Community-based enumeration ----
# One record per child: relation (Phase 1 and Somalia Phase 2 forms), dose
# counts and antigen indicators. The Somalia Phase 2 form asked about card
# possession, so the source of a reported non-receipt is completed there. CBE
# forms were filled on paper and keyed in later; the date on each record is
# the field day recovered for it (missing where it could not be recovered).
# clean each CBE form
for (sample in str_subset(study_samples, "cbe")) {
  cbe <- raw[[sample]]
  if ("relation_child" %in% names(cbe)) cbe$relation_child <- classify_relation(cbe$relation_child)
  cbe <- clean_children(cbe)
  if (sample == "som_phase2_cbe") cbe <- code_answer_source(cbe, holder_field = "vaccination_card")
  frames[[sample]] <- cbe
}

## Children reported through the network module ----
# The same dose and antigen coding, applied to each child a respondent reported
# about in their network.
# classify proxy-reported doses and join each respondent to their network children
network <- map(set_names(setdiff(names(main), "som_phase2_rds_all_enumerators")), \(sample) {
  children <-
    raw[[sample]]$network |>
    classify_doses(type = "network") |>
    derive_antigens()
  left_join(main[[sample]], select(children, record_ix, !any_of(names(main[[sample]]))), by = "record_ix")
})

## Analytic exclusions ----
# Flag the children left out of the directly reported analysis; the full frame
# is kept for the fieldwork tables.
# start with every child in the analysis
frames <- map(frames, \(frame) mutate(frame, analytic = TRUE))

# Ethiopia Phase 2 VBS screened for caregivers of children under five, so its
# analysis keeps children aged 12-23 completed months at interview (age in
# months = days / 30.4375, given on the form as age_months)
# flag children outside 12-23 months
eth_p2_vbs <- frames$eth_phase2_vbs
age_months <- eth_p2_vbs$age_months
screened <- eth_p2_vbs$consent %in% "Enrolled" & eth_p2_vbs$participant_status %in% "Enrolled with Eligible Child"
eth_p2_vbs$analytic <- !screened | (!is.na(age_months) & age_months >= 12 & age_months < 24)
frames$eth_phase2_vbs <- eth_p2_vbs

# one Somalia Phase 1 VBS form lists nine children born over eleven months; its
# children are left out of the directly reported analysis (the caregiver stays
# in the network module)
# find the form listing nine children and flag its children
listed_nine <-
  raw$som_phase1_vbs$child |>
  count(record_ix) |>
  filter(n == 9) |>
  pull(record_ix)
frames$som_phase1_vbs$analytic[frames$som_phase1_vbs$record_ix %in% listed_nine] <- FALSE

## Field staff hours ----
# One row per person, day (or week) and activity: the time-motion segments
# recorded by enumerators and supervisors, the hours Somalia supervisors logged
# in their tracker, and the weekly timesheets of the coordinating staff.
# time-motion segments
segment_hours <-
  time_motion$segments |>
  transmute(
    source = "time-motion sheet", phase, country, method, role, person = person_id, date,
    activity = activity_category, segment, segment_order, hours, structural_absent, flag_negative, flag_long, flag_any
  )

# Somalia supervisor tracker
tracker_hours <-
  time_motion$supervisor_tracker |>
  transmute(source = "supervisor tracker", phase, country = "som", method, role = "supervisor", person, activity = "field work", hours)

# coordinator weekly timesheets
coordinator_hours <-
  time_motion$coordinator_weeks |>
  transmute(
    source = "coordinator timesheet", country = str_to_lower(country), method = str_to_lower(method),
    role = str_to_lower(role), person, organisation, date = week_start, activity = str_to_lower(activity), hours,
    hours_as_entered, activity_inferred, method_inferred
  )

# stack the three sources
hours <- bind_rows(segment_hours, tracker_hours, coordinator_hours)

## Cluster places ----
# The named place of each fielded HH and CBE cluster, with the place's first
# letter in upper case as printed.
# keep the fielded clusters
cluster_places <-
  cluster_places |>
  filter(fielded) |>
  mutate(place_name = str_replace(place_name, "^no named", "No named")) |>
  select(country, phase, method, cluster_id, place_name, place_type)

## Record counts ----
# Records in each KoboToolbox export before cleaning, one row per form.
# keep the form and its count
kobo_record_counts <- select(kobo_record_counts, country, phase, method, n_records)

## Tally sheets ----
# One row per supervisor tally sheet with the cluster number and the household
# counts: households found absent in the Ethiopia Phase 1 household survey, and
# households found absent, ineligible, refusing or completed in each CBE sample.
# keep the counts as numbers and order the sheets
tallies <-
  tallies |>
  mutate(across(c(absent, ineligible, refused, completed), as.numeric)) |>
  arrange(phase, country, method, cluster)

## Card photographs ----
# One row per household-survey child with at least one coded photograph: the
# number of photos (and of vaccination-card photos) and the status the photos
# show for Penta1, Penta3, OPV3 and MCV1 (yes, no or unknown).
# keep the status columns as plain text
card_photo_status <-
  card_photo_status |>
  as_tibble() |>
  mutate(across(c(penta1, penta3, opv3, mcv1), as.character))

## Budget line items ----
# One row per budget line, country, phase and method, in the currency the
# partner budgeted in: birr (AAU), US dollars (SORDI and the tablet setting)
# or Swiss francs by calendar-year column (UNIGE). Phase 1 lines are at 2025
# prices.
# lower-case the country and method codes used across the clean files
cost_line_items <-
  cost_line_items |>
  as_tibble() |>
  mutate(country = str_to_lower(country), method = str_to_lower(method), study_phase = as.integer(study_phase))

## Shared budget lines ----
# The budget lines AAU and SORDI booked once for a country and phase, which
# each method's line items repeat: country, phase, organisation, item and the
# amount in the original currency.
# lower-case the country codes used across the clean files
cost_shared_items <-
  cost_shared_items |>
  as_tibble() |>
  mutate(country = str_to_lower(country), study_phase = as.integer(study_phase))

## Research-team cost ----
# The UNIGE research team's Phase 1 timesheet hours at their base day rates
# (8-hour day, US dollars at 2025 prices, before overhead), summed over the
# team by the country (Ethiopia, Somalia or both) and method (or cross-cutting)
# the hours were logged against.
# keep the cost as a number
research_team_cost <- research_team_cost |> as_tibble() |> mutate(day_rate_cost_usd = as.numeric(day_rate_cost_usd))

## Exchange rates and price index ----
# Mean mid-market rates over each fieldwork window or calendar year, and the
# US CPI-U means that take 2025 prices to 2026 prices.
# keep the rates as numbers
fx_cpi <- fx_cpi |> as_tibble() |> mutate(value = as.numeric(value))

# Save ----

# order the samples and write the clean files
frames <- frames[c(study_samples, "som_phase2_rds_all_enumerators")]
network <- network[intersect(study_samples, names(network))]

write_rds(frames, "data/clean/frames.rds")
write_rds(network, "data/clean/network.rds")
write_rds(tracking, "data/clean/rds_tracking.rds")
write_rds(cluster_weights, "data/clean/cluster_weights.rds")
write_rds(cluster_places, "data/clean/cluster_places.rds")
write_rds(kobo_record_counts, "data/clean/kobo_record_counts.rds")
write_rds(tallies, "data/clean/tallies.rds")
write_rds(card_photo_status, "data/clean/card_photo_status.rds")
write_rds(hours, "data/clean/hours.rds")
write_rds(cost_line_items, "data/clean/cost_line_items.rds")
write_rds(cost_shared_items, "data/clean/cost_shared_items.rds")
write_rds(research_team_cost, "data/clean/research_team_cost.rds")
write_rds(fx_cpi, "data/clean/fx_cpi.rds")
