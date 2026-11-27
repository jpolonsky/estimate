# 03_samples.R
# Study areas, fieldwork periods, recruitment and sample composition
# Reproduces: Table 7, Figures 2, 6, 7a, 7b, 8 and 9, Annex Tables D1 and F2-F5

# Setup ----

pacman::p_load(sf, ggrepel, ggspatial, ggtext, patchwork, tidyverse)

source("code/functions.R")

sf_use_s2(FALSE) |> suppressMessages()

# Data ----

frames <- read_rds("data/clean/frames.rds")
cluster_weights <- read_rds("data/clean/cluster_weights.rds")
tallies <- read_rds("data/clean/tallies.rds")
kobo_record_counts <- read_rds("data/clean/kobo_record_counts.rds")

# map layers and the towns and field bases drawn on the map
map_layers <- map(
  set_names(c("countries", "admin_units", "settlements", "roads", "rivers")),
  \(layer) st_read("data/geo/map_layers.gpkg", layer = layer, quiet = TRUE)
)
places <- read_csv("data/geo/map_places.csv", show_col_types = FALSE)

# the fourteen samples, in report order
samples <- study_samples
frames <- frames[samples]

# Analysis ----

## Respondents and cohorts ----
# One row per respondent (every CBE record is one), and the directly reported
# children each estimate is based on.
# take one row per respondent
respondents <- imap(frames, \(frame, sample) if (str_detect(sample, "cbe")) frame else distinct(frame, record_ix, .keep_all = TRUE))

# keep each sample's directly reported children
cohorts <- imap(frames, direct_cohort)

## Fieldwork periods ----
# First and last interview date of each sample (for CBE, the field day
# recovered for each record).
# first and last date of each sample
fieldwork <-
  imap(respondents, \(d, sample) tibble(sample = sample, start = min(as.Date(d$date), na.rm = TRUE), end = max(as.Date(d$date), na.rm = TRUE))) |>
  list_rbind() |>
  left_join(sample_parts(samples), by = "sample")

## Recruitment ----
# Respondents by disposition: absent (from the tally sheets, where kept),
# screened ineligible, refused and enrolled; approached is their sum. CBE
# forms record enrolled households only, so their ineligible, absent and
# refused households come from the tally sheets. The Somalia Phase 1 CBE
# sheets of 8 of the 30 clusters were not received, and their counts are
# imputed from the mean of the other clusters.
# sum the tally sheets of each sample, imputing the clusters without a sheet
tally_totals <- function(sample) {
  parts <- sample_parts(sample)
  sheets <- filter(tallies, country == parts$country, phase == parts$phase, method == parts$method)
  if (nrow(sheets) == 0) return(tibble(absent = NA_real_, ineligible = NA_real_, refused = NA_real_))
  n_missing <- if (parts$method == "cbe") n_distinct(cohorts[[sample]]$cluster) - n_distinct(sheets$cluster) else 0
  sheets |>
    summarise(across(c(absent, ineligible, refused), \(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE))) |>
    mutate(across(everything(), \(x) x + round(x / n_distinct(sheets$cluster) * n_missing)))
}

# count each sample's respondents by disposition
dispositions <- imap(respondents, \(d, sample) {
  if (str_detect(sample, "cbe")) {
    tally <- tally_totals(sample)
    return(tibble(sample, absent = tally$absent, ineligible = tally$ineligible, refused = tally$refused, enrolled = nrow(d), reached_consent = nrow(d) + tally$refused))
  }
  tally <- if (sample == "eth_phase1_hh") tally_totals(sample) else tibble(absent = NA_real_)
  tibble(
    sample,
    absent = tally$absent,
    ineligible = sum(d$consent == "Ineligible"),
    refused = sum(d$consent == "Refused"),
    enrolled = sum(d$consent == "Enrolled"),
    # caregivers who consented but reported no eligible child also reached consent
    reached_consent = enrolled + refused + sum(d$consent == "Ineligible" & !is.na(d$hh_child_num))
  )
}) |>
  list_rbind()

# directly reporting caregivers, their children and the share of children
# with a recorded status for all three primary antigens
achieved <- imap(cohorts, \(cohort, sample) {
  tibble(
    sample,
    direct = if (str_detect(sample, "cbe")) nrow(cohort) else n_distinct(cohort$record_ix),
    children = nrow(cohort),
    complete = mean(!is.na(cohort$penta_3p_clean) & !is.na(cohort$opv_3p_clean) & !is.na(cohort$measles_1p_clean))
  )
}) |>
  list_rbind()

# sampling units: clusters visited (and selected), seeds enrolled, venues with
# a directly reporting caregiver
units <- imap(cohorts, \(cohort, sample) {
  parts <- sample_parts(sample)
  switch(
    parts$method,
    hh = tibble(sample, unit = "clusters", visited = n_distinct(cohort$cluster_id)),
    cbe = tibble(sample, unit = "clusters", visited = n_distinct(cohort$cluster)),
    rds = tibble(sample, unit = "seeds", visited = sum(distinct(cohort, record_ix, .keep_all = TRUE)$seed_participant == "Yes", na.rm = TRUE)),
    vbs = tibble(sample, unit = "venues", visited = n_distinct(cohort$int_location))
  )
}) |>
  list_rbind() |>
  left_join(
    cluster_weights |>
      count(country, phase, method, name = "selected") |>
      transmute(sample = str_c(country, "_phase", phase, "_", method), selected),
    by = "sample"
  )

# assemble one row per sample
recruitment <-
  sample_parts(samples) |>
  left_join(dispositions, by = "sample") |>
  left_join(achieved, by = "sample") |>
  left_join(units, by = "sample") |>
  mutate(
    approached = coalesce(absent, 0) + ineligible + refused + enrolled,
    refusal_pct = 100 * refused / reached_consent,
    imputed = sample == "som_phase1_cbe"
  )

## Participant flow ----
# The steps from approach to children directly reported. Eligible after
# screening is approached minus absent and ineligible. Enrolled VBS caregivers
# who are not directly reporting either had no eligible child (network module
# only) or had their child records excluded: children outside 12-23 months in
# Ethiopia Phase 2 VBS, and one Somalia Phase 1 VBS form listing nine children.
# count enrolled VBS caregivers without an eligible child, and those whose
# child records were excluded
vbs_losses <-
  imap(respondents[str_detect(samples, "vbs")], \(d, sample) {
    eligible <- filter(d, consent == "Enrolled", participant_status == "Enrolled with Eligible Child")
    tibble(
      sample,
      network_only = sum(d$consent == "Enrolled") - nrow(eligible),
      excluded = sum(!eligible$record_ix %in% cohorts[[sample]]$record_ix)
    )
  }) |>
  list_rbind()

# stage counts
flow_stages <-
  recruitment |>
  transmute(
    sample, country, phase, method,
    `Approached` = approached,
    `Eligible after screening` = approached - coalesce(absent, 0) - ineligible,
    `Enrolled caregivers` = enrolled,
    `Directly reporting caregivers` = direct,
    `Children directly reported` = children
  ) |>
  pivot_longer(Approached:`Children directly reported`, names_to = "stage_label", values_to = "n") |>
  mutate(stage = match(stage_label, unique(stage_label)))

# losses at each step, in the side box after that step
flow_losses <-
  recruitment |>
  left_join(vbs_losses, by = "sample") |>
  transmute(
    sample, country, phase, method,
    `1_Absent` = absent,
    `1_Ineligible` = ineligible,
    `2_Refused` = refused,
    `3_Network only` = if_else(network_only > 0, network_only, NA),
    `3_Outside 12–23 mo` = if_else(sample == "eth_phase2_vbs" & excluded > 0, excluded, NA),
    `3_Form excluded` = if_else(sample == "som_phase1_vbs" & excluded > 0, excluded, NA)
  ) |>
  pivot_longer(starts_with(c("1_", "2_", "3_")), names_to = "item", values_to = "n", values_drop_na = TRUE) |>
  separate_wider_delim(item, "_", names = c("step", "item")) |>
  mutate(step = as.integer(step))

## Sample composition ----
# Unweighted share of directly reporting respondents (of children, for child
# sex) in each category of seven characteristics. Wealth is "High asset" when
# the respondent owns more than half of the assets asked about in every
# sample of that country and phase.
# group religions
religion_group <- c("Muslim" = "Muslim", "Orthodox" = "Christian", "Protestant" = "Christian", "Catholic" = "Christian", "Other" = "Other")

# directly reporting respondents with the derived characteristics
composition_respondents <- imap(cohorts, \(cohort, sample) {
  d <- if (str_detect(sample, "cbe")) cohort else distinct(cohort, record_ix, .keep_all = TRUE)
  age <- if ("respondent_age" %in% names(d)) d$respondent_age else NA
  d$age_band <- case_when(between(age, 18, 24) ~ "18-24", between(age, 25, 34) ~ "25-34", age >= 35 ~ "35+")
  d$education_stage <- if ("education" %in% names(d)) unname(education_stages[as.character(d$education)]) else NA
  d$religion_grp <- if ("religion" %in% names(d)) unname(religion_group[as.character(d$religion)]) else NA
  d$wealth_grp <- wealth_group(d, common_assets(cohorts[samples], sample))
  d
})

# the seven characteristics and the field each is read from
characteristics <- c(
  respondent_sex = "respondent_sex", age_band = "age_band", residence = "residence", education = "education_stage",
  wealth_grp = "wealth_grp", religion = "religion_grp", child_sex = "child_sex"
)

# share in each category; a characteristic needs at least two categories
composition <-
  map(samples, \(sample) {
    imap(characteristics, \(field, variable) {
      d <- if (variable == "child_sex") cohorts[[sample]] else composition_respondents[[sample]]
      if (!field %in% names(d)) return(NULL)
      counts <- table(d[[field]][!is.na(d[[field]])])
      if (length(counts) < 2) return(NULL)
      tibble(sample, variable, category = names(counts), n = as.integer(counts), pct = 100 * n / sum(n))
    }) |>
      list_rbind()
  }) |>
  list_rbind() |>
  left_join(sample_parts(samples), by = "sample")

## Distance from the reference ----
# Total variation distance (half the sum of absolute differences in
# proportions) of each sample from the household survey in Phase 1, and from
# the equal-weight mean of the methods in Phase 2, averaged over the
# characteristics recorded in every sample compared.
# calculate total variation distance between two sets of proportions
tvd <- function(p, q) {
  categories <- union(names(p), names(q))
  0.5 * sum(abs(replace_na(p[categories], 0) - replace_na(q[categories], 0)))
}
proportions <- \(rows) set_names(rows$pct / 100, rows$category)

# Phase 1: distance of each method from the household survey
distance_phase1 <-
  composition |>
  filter(phase == 1) |>
  group_by(country, variable) |>
  group_modify(\(rows, key) {
    hh <- filter(rows, method == "hh")
    if (nrow(hh) == 0) return(tibble())
    rows |>
      filter(method != "hh") |>
      group_by(method) |>
      summarise(tvd = tvd(proportions(pick(category, pct)), proportions(hh)), .groups = "drop")
  }) |>
  ungroup() |>
  mutate(phase = 1)

# Phase 2: distance of each method from the mean of the methods
distance_phase2 <-
  composition |>
  filter(phase == 2) |>
  group_by(country, variable) |>
  filter(n_distinct(method) >= 2) |>
  group_modify(\(rows, key) {
    pooled <- rows |> group_by(category) |> summarise(pct = sum(pct) / n_distinct(rows$method))
    rows |>
      group_by(method) |>
      summarise(tvd = tvd(proportions(pick(category, pct)), proportions(pooled)), .groups = "drop")
  }) |>
  ungroup() |>
  mutate(phase = 2)

# average over the characteristics every compared sample recorded
mean_distance <-
  bind_rows(distance_phase1, distance_phase2) |>
  group_by(country, phase, variable) |>
  mutate(n_methods = n_distinct(method)) |>
  group_by(country, phase) |>
  filter(n_methods == max(n_methods)) |>
  group_by(country, phase, method) |>
  summarise(tvd = mean(tvd), n_characteristics = n(), .groups = "drop")

## Composition tests ----
# A chi-squared test of each characteristic between a method and the
# household survey (Phase 1) and between each pair of methods (Phase 2), on
# the unweighted counts; Phase 2 pairs also get their total variation
# distance.
# the recorded values of one characteristic in one sample
characteristic_values <- function(sample, variable) {
  d <- if (variable == "child_sex") cohorts[[sample]] else composition_respondents[[sample]]
  field <- characteristics[[variable]]
  if (!field %in% names(d)) return(NULL)
  values <- as.character(d[[field]])
  values <- values[!is.na(values)]
  if (n_distinct(values) < 2) NULL else values
}

# chi-squared p-value between two samples (with Yates's correction on a 2 x 2 table)
chisq_p <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NA_real_)
  counts <- table(rep(c("a", "b"), c(length(a), length(b))), c(a, b))
  suppressWarnings(chisq.test(counts)$p.value)
}

# Phase 1: each method against the household survey
composition_tests_phase1 <-
  distance_phase1 |>
  mutate(p_value = pmap_dbl(list(country, method, variable), \(country, method, variable) {
    chisq_p(characteristic_values(str_c(country, "_phase1_", method), variable), characteristic_values(str_c(country, "_phase1_hh"), variable))
  }))

# Phase 2: each pair of methods
composition_pairs <-
  tibble(country = rep(c("eth", "som"), each = 3), method_a = rep(c("cbe", "cbe", "rds"), 2), method_b = rep(c("rds", "vbs", "vbs"), 2)) |>
  expand_grid(variable = names(characteristics))
composition_tests_phase2 <-
  composition_pairs |>
  mutate(
    tvd = pmap_dbl(list(country, method_a, method_b, variable), \(country, method_a, method_b, variable) {
      rows <- filter(composition, phase == 2, country == .env$country, variable == .env$variable)
      a <- filter(rows, method == method_a)
      b <- filter(rows, method == method_b)
      if (nrow(a) == 0 || nrow(b) == 0) NA_real_ else tvd(proportions(a), proportions(b))
    }),
    p_value = pmap_dbl(list(country, method_a, method_b, variable), \(country, method_a, method_b, variable) {
      chisq_p(characteristic_values(str_c(country, "_phase2_", method_a), variable), characteristic_values(str_c(country, "_phase2_", method_b), variable))
    })
  )

## Records excluded in cleaning ----
# Records in each KoboToolbox export, those removed in cleaning (duplicates,
# test forms and the other exclusions of the data-cleaning rules) and those
# retained: every respondent screened, refused or enrolled, and for CBE every
# enrolled caregiver.
# compare the export with the retained records
cleaning <-
  kobo_record_counts |>
  mutate(sample = str_c(country, "_phase", phase, "_", method)) |>
  left_join(tibble(sample = samples, retained = map_int(respondents[samples], nrow)), by = "sample") |>
  mutate(excluded = n_records - retained) |>
  arrange(match(sample, samples))

## Cluster design ----
# The design parameters of each household-survey and CBE sample, as recorded
# in the cluster design files, with the clusters selected at design and
# visited in fieldwork.
# design parameters of each sample
cluster_design <- tribble(
  ~sample, ~surface, ~age_band, ~target, ~polygon, ~seed_rule, ~redraw,
  "eth_phase1_hh", "HRSL, 30 m", "Under five", "Approximately 120", "Growth cells", "Not recorded", "Not recorded",
  "eth_phase1_cbe", "HRSL, 30 m", "Under five", "Approximately 240", "Growth cells", "Not recorded", "Not recorded",
  "eth_phase2_cbe", "HRSL, 30 m", "Under five", "Approximately 240", "Growth cells", "Not recorded", "Not recorded",
  "som_phase1_hh", "WorldPop 2025 constrained, 25 m", "Under one", "Approximately 25", "Growth cells and the 100 nearest buildings", "Within settlements", "Not recorded",
  "som_phase1_cbe", "WorldPop 2025 constrained, 25 m", "Under one", "Approximately 50", "Growth cells and the 200 nearest buildings", "Within settlements", "Not recorded",
  "som_phase2_cbe", "WorldPop 2025 constrained, 25 m", "Under one", "Approximately 50", "Growth cells and the 200 nearest buildings",
  "At least 200 buildings within 500 m, inside the GRID3 settlement extents", "Overlap, or fewer than 180 buildings enclosed; at most 40 attempts"
)

# add the clusters selected and visited
cluster_design <- left_join(cluster_design, select(units, sample, selected, visited), by = "sample")

# Tables ----

## Table 7 ----
# Recruitment and achieved samples, a header row per phase; an asterisk marks
# the imputed Somalia Phase 1 CBE counts.
# format the recruitment counts
star <- \(x, imputed) if_else(imputed, str_c(x, "*"), x)
table_07_rows <-
  recruitment |>
  mutate(
    Country = if_else(method == if_else(phase == 1, "hh", "cbe"), country_labels[country], ""),
    Method = method_labels[method],
    Approached = star(format_count(approached), imputed),
    Absent = if_else(is.na(absent), "–", star(format_count(absent), imputed)),
    `Screened ineligible` = star(format_count(ineligible), imputed),
    Refused = star(format_count(refused), imputed),
    `Refusal (%)` = star(format_number(refusal_pct), imputed),
    `Enrolled caregivers` = format_count(enrolled),
    `Directly reporting caregivers` = format_count(direct),
    `Children directly reported` = format_count(children),
    `Complete on all three primary antigens (%)` = format_number(100 * complete),
    `Children per directly reporting caregiver` = format_number(children / direct, 2),
    `Sampling units (visited / selected)` = if_else(
      !is.na(selected) & selected != visited, str_c(visited, " of ", selected, " ", unit), str_c(visited, " ", unit)
    )
  )

# add a header row for each phase
table_07 <-
  map(1:2, \(p) {
    bind_rows(
      tibble(Country = str_c("Phase ", p)),
      table_07_rows |> filter(phase == p) |> select(Country:`Sampling units (visited / selected)`)
    )
  }) |>
  list_rbind() |>
  mutate(across(everything(), \(x) replace_na(as.character(x), "")))

save_table(table_07, "output/tables/table_07.csv")

## Annex Table D1 ----
# Cluster-design parameters of each household-survey and CBE sample.
# label the samples and keep the printed columns
table_D01 <-
  cluster_design |>
  left_join(sample_parts(samples), by = "sample") |>
  transmute(
    Sample = str_c(country_labels[country], " Phase ", phase, " ", method_labels[method]),
    `Population surface` = surface,
    `Age band of the target` = age_band,
    `Growth target (children)` = target,
    Polygon = polygon,
    `Seed eligibility` = seed_rule,
    `Redraw rule` = redraw,
    `Clusters selected` = selected,
    `Clusters visited` = visited
  )

save_table(table_D01, "output/tables/table_D01.csv")

## Annex Table F2 ----
# Records in the export, excluded in cleaning and retained, by phase, country
# and method.
# format the counts
table_F02 <-
  cleaning |>
  transmute(
    `Phase and country` = str_c("Phase ", phase, ", ", country_labels[country]),
    Method = method_labels[method],
    `Raw records` = format_count(n_records),
    `Excluded in cleaning` = format_count(excluded),
    Retained = format_count(retained)
  ) |>
  mutate(`Phase and country` = if_else(`Phase and country` == lag(`Phase and country`, default = ""), "", `Phase and country`))

save_table(table_F02, "output/tables/table_F02.csv")

## Annex Table F3 ----
# Phase 1 composition (% of respondents; child sex, % of children), one column
# per country and method under a country header; a dagger marks a value more
# than 10 pp from the household survey.
# the categories shown and their labels
f3_rows <- tibble(
  variable = c("respondent_sex", "age_band", "age_band", "age_band", "child_sex", "residence", "education", "education", "education", "education", "wealth_grp", "religion", "religion"),
  category = c("Female", "18-24", "25-34", "35+", "Female", "Urban", "No formal education", "Primary", "Secondary", "Higher", "High asset", "Muslim", "Christian"),
  Characteristic = c(
    "Respondent sex: female", "Age band: 18–24 years", "Age band: 25–34 years", "Age band: 35 years or more", "Child sex: female",
    "Residence: urban", "Education: no formal education", "Education: primary", "Education: secondary", "Education: higher",
    "Household wealth: high asset", "Religion: Muslim", "Religion: Christian"
  )
)

# format each share, marking those far from the household survey
f3_values <-
  composition |>
  filter(phase == 1) |>
  inner_join(f3_rows, by = c("variable", "category")) |>
  group_by(country, variable, category) |>
  mutate(far = method != "hh" & abs(pct - pct[method == "hh"]) > 10) |>
  ungroup() |>
  mutate(value = str_c(format_number(pct), if_else(far, " †", ""))) |>
  select(Characteristic, country, method, value)

# one column per country and method, with the country header row
f3_columns <- expand_grid(country = c("eth", "som"), method = method_order)
table_F03 <-
  map(seq_len(nrow(f3_columns)), \(k) {
    values <- filter(f3_values, country == f3_columns$country[[k]], method == f3_columns$method[[k]])
    replace_na(values$value[match(f3_rows$Characteristic, values$Characteristic)], "–")
  }) |>
  set_names(str_c(f3_columns$country, "_", f3_columns$method)) |>
  as_tibble() |>
  mutate(Characteristic = f3_rows$Characteristic, .before = 1)
f3_header <- c("", "Ethiopia", "", "", "", "Somalia", "", "", "")
table_F03 <- rbind(c("Characteristic", method_labels[f3_columns$method]), as.matrix(table_F03))
colnames(table_F03) <- f3_header

save_table(as_tibble(table_F03, .name_repair = "minimal"), "output/tables/table_F03.csv")

## Annex Tables F4 and F5 ----
# Total variation distance with the chi-squared p-value in parentheses, by
# characteristic: from the household survey in Phase 1 (with the mean over
# the characteristics compared) and between each pair of methods in Phase 2.
# format a distance and its p-value
format_tvd <- \(tvd, p) if_else(is.na(tvd), "–", str_c(sprintf("%.3f", tvd), " (", if_else(p < 0.001, "<0.001", sprintf("%.3f", p)), ")"))
characteristic_labels <- c(
  respondent_sex = "Respondent sex", age_band = "Age band", child_sex = "Child sex", residence = "Residence",
  education = "Education", wealth_grp = "Household wealth", religion = "Religion"
)

# Phase 1, against the household survey
f4_columns <- expand_grid(country = c("eth", "som"), method = c("cbe", "rds", "vbs"))
table_F04 <-
  map(seq_len(nrow(f4_columns)), \(k) {
    rows <- filter(composition_tests_phase1, country == f4_columns$country[[k]], method == f4_columns$method[[k]])
    i <- match(names(characteristic_labels), rows$variable)
    mean_row <- filter(mean_distance, phase == 1, country == f4_columns$country[[k]], method == f4_columns$method[[k]])
    c(format_tvd(rows$tvd[i], rows$p_value[i]), sprintf("%.3f", mean_row$tvd))
  }) |>
  set_names(str_c(country_labels[f4_columns$country], " ", method_labels[f4_columns$method])) |>
  as_tibble() |>
  mutate(Characteristic = c(characteristic_labels, "Mean over the characteristics compared"), .before = 1)

save_table(table_F04, "output/tables/table_F04.csv")

# Phase 2, between each pair of methods
table_F05 <-
  composition_tests_phase2 |>
  mutate(
    pair = str_c(country_labels[country], " ", method_labels[method_a], " and ", method_labels[method_b]),
    value = format_tvd(tvd, p_value),
    Characteristic = characteristic_labels[variable]
  ) |>
  select(Characteristic, pair, value) |>
  pivot_wider(names_from = pair, values_from = value) |>
  arrange(match(Characteristic, characteristic_labels))

save_table(table_F05, "output/tables/table_F05.csv")

# Figures ----

## Figure 2 ----
# Study areas: one map per study area (Ethiopia left, Somalia right, Phase 1
# above Phase 2) around a regional map and a key. The study area is the
# woreda in Ethiopia and, where study-area boundaries are not public, the
# enclosing district in Somalia (dashed).
# panel centre (degrees), half-width (km) and scale-bar corner
sites <- tibble(
  panel = c("eth1", "som1", "eth2", "som2"), num = 1:4, country = c("eth", "som", "eth", "som"), phase = c(1, 1, 2, 2),
  site = c("Meskan and Butajira", "Galkayo South", "Yilmana Densa", "Afmadow"),
  clon = c(38.335, 47.10, 37.47, 41.80), clat = c(8.03, 6.50, 11.355, 0.55), half_km = c(31, 92, 34, 122),
  scale_loc = c("br", "br", "tr", "tr")
)
study_admin <- list(eth1 = c("Meskan", "Butajira town"), eth2 = c("Yilmana Densa", "Adete town"), som1 = "Gaalkacyo", som2 = "Afmadow")
country_colours <- c(eth = "#B5443C", som = "#1F6F8B")
country_tints <- c(eth = "#F2E0DD", som = "#DAE9EF")
map_style <- list(
  land = "#F7F6F3", water = "#DDE8EE", border = "#CDC9C0", road = "#C9C3B8", main_road = "#AFA898", river = "#9DBCCB",
  settlement = c("Hamlet" = "#CDCDD2", "Small Settlement Area" = "#A4A4AC", "Built-up Area" = "#7C7C86"), settlement_som = "#9A9AA2"
)
laea <- \(lon, lat) sprintf("+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m", lat, lon)

# outline of each study area
study_areas <-
  map(sites$panel, \(k) {
    units <- filter(map_layers$admin_units, panel == k, name %in% study_admin[[k]])
    st_sf(panel = k, geometry = suppressMessages(st_union(units)))
  }) |>
  list_rbind() |>
  st_as_sf()

# map theme and scale bar
theme_map <- function() {
  theme_void(base_family = report_font) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = map_style$water, colour = NA),
      panel.border = element_rect(fill = NA, colour = "#B5B1A8", linewidth = 0.35),
      plot.title = element_markdown(size = 7.6, colour = ink[["primary"]], lineheight = 1.3, margin = margin(0, 0, 3, 0)),
      plot.title.position = "panel", plot.margin = margin(2, 3, 4, 3)
    )
}
scale_bar <- function(location) {
  annotation_scale(
    location = location, style = "ticks", width_hint = 0.28, text_family = report_font, text_cex = 0.5, line_width = 0.4,
    line_col = ink[["secondary"]], text_col = ink[["secondary"]], tick_height = 0.5, height = unit(1, "mm"),
    pad_x = unit(1.8, "mm"), pad_y = unit(1.6, "mm")
  )
}

# create one study-area map
site_map <- function(k) {
  s <- filter(sites, panel == k)
  colour <- country_colours[[s$country]]
  crs <- laea(s$clon, s$clat)
  h <- s$half_km * 1000
  # Lake Tana is land in the countries layer, so the woredas carry the land around Bahir Dar
  countries <- if (k == "eth2") filter(map_layers$countries, admin != "Ethiopia") else map_layers$countries
  settlements <- filter(map_layers$settlements, panel == k)
  settlements$fill <- if (s$country == "eth") unname(map_style$settlement[settlements$type]) else map_style$settlement_som
  roads <- filter(map_layers$roads, panel == k)
  major <- roads$highway %in% c("trunk", "primary")
  area <- filter(study_areas, panel == k)
  pl <- filter(places, panel == k)
  xy <- st_coordinates(st_transform(st_as_sf(pl, coords = c("lon", "lat"), crs = 4326), crs))
  pl <- mutate(
    pl,
    x = xy[, 1], y = xy[, 2], lx = x + label_dx * 0.0255 * h, ly = y + label_dy * 0.085 * h,
    hj = case_when(label_dx < 0 ~ 1, label_dx > 0 ~ 0, .default = 0.5), face = if_else(role == "town", "plain", "bold")
  )
  badge <- c(x = -0.89 * h, y = 0.89 * h)
  ggplot() +
    geom_sf(data = countries, fill = map_style$land, colour = map_style$border, linewidth = 0.3) +
    geom_sf(data = filter(map_layers$admin_units, panel == k), fill = map_style$land, colour = map_style$border, linewidth = 0.2) +
    (if (s$country == "eth") geom_sf(data = area, fill = country_tints[[s$country]], colour = NA)) +
    geom_sf(data = roads[!major, ], colour = map_style$road, linewidth = 0.18) +
    geom_sf(data = roads[major, ], colour = map_style$main_road, linewidth = 0.4) +
    geom_sf(data = filter(map_layers$rivers, panel == k), colour = map_style$river, linewidth = 0.55) +
    geom_sf(data = settlements, aes(fill = I(fill), colour = I(fill)), linewidth = if (s$country == "eth") 0 else 0.12) +
    (if (s$country == "eth") {
      geom_sf(data = area, fill = NA, colour = colour, linewidth = 0.55)
    } else {
      geom_sf(data = area, fill = NA, colour = ink[["primary"]], linewidth = 0.45, linetype = "32")
    }) +
    geom_point(data = filter(pl, role == "town"), aes(x, y), shape = 21, fill = "white", colour = ink[["secondary"]], size = 1.5, stroke = 0.5) +
    geom_point(data = filter(pl, role == "vbs_venue"), aes(x, y), shape = 24, fill = method_colours[["vbs"]], colour = "white", size = 2.4, stroke = 0.55) +
    geom_point(data = filter(pl, role == "field_base"), aes(x, y), shape = 22, fill = ink[["primary"]], colour = "white", size = 2.4, stroke = 0.55) +
    geom_point(data = filter(pl, role == "field_base_and_vbs"), aes(x, y), shape = 22, fill = method_colours[["vbs"]], colour = ink[["primary"]], size = 2.4, stroke = 0.65) +
    # place names with a white halo at fixed offsets
    geom_text_repel(
      data = pl, aes(x = lx, y = ly, label = name, hjust = hj, fontface = face), family = report_font, size = 2.35, colour = ink[["primary"]],
      bg.colour = "white", bg.r = 0.2, seed = 1, force = 0, force_pull = 0, max.time = 0, max.iter = 0, min.segment.length = Inf, point.size = NA
    ) +
    annotate("point", x = badge[["x"]], y = badge[["y"]], shape = 21, fill = colour, colour = "white", size = 4.6, stroke = 0.6) +
    annotate("text", x = badge[["x"]], y = badge[["y"]], label = s$num, colour = "white", fontface = "bold", size = 2.4, family = report_font) +
    coord_sf(crs = crs, default_crs = crs, xlim = c(-h, h), ylim = c(-h, h), expand = FALSE, datum = NA) +
    scale_bar(s$scale_loc) +
    labs(title = sprintf(
      "**%s**<br><span style='color:%s'>%s · Phase %d, %s</span>", s$site, ink[["muted"]],
      country_labels[[s$country]], s$phase, if (s$phase == 1) "accessible" else "inaccessible"
    )) +
    theme_map()
}

# create the regional map, with Somaliland drawn as part of Somalia
region_map <- function() {
  focus_of <- c(Ethiopia = "eth", Somalia = "som", Somaliland = "som")
  focus <-
    map_layers$countries |>
    filter(admin %in% names(focus_of)) |>
    mutate(country = focus_of[admin]) |>
    group_by(country) |>
    summarise(.groups = "drop") |>
    suppressMessages()
  centres <- suppressWarnings(st_point_on_surface(study_areas)) |> left_join(select(sites, panel, num, country), by = "panel")
  as_points <- \(d) st_as_sf(d, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  capitals <- tibble(name = c("Addis Ababa", "Mogadishu"), lon = c(38.74, 45.334), lat = c(9.03, 2.041))
  capital_labels <- tibble(name = capitals$name, lon = capitals$lon + c(0.45, -0.45), lat = capitals$lat + 0.05, hj = c(0, 1))
  country_names <- tibble(name = c("ETHIOPIA", "SOMALIA"), lon = c(40.4, 45.2), lat = c(6.9, 3.9), country = c("eth", "som"))
  neighbours <- tibble(
    name = c("Eritrea", "Djibouti", "Kenya", "South\nSudan", "Sudan", "Yemen"),
    lon = c(38.3, 42.75, 38.2, 33.4, 34.4, 46.8), lat = c(15.0, 11.7, 1.4, 6.0, 12.6, 15.0)
  )
  waters <- tibble(name = c("Red Sea", "Gulf of Aden", "Indian\nOcean", "Lake\nTurkana"), lon = c(42.1, 47.6, 49.6, 36.0), lat = c(14.4, 12.2, 3.0, 3.1))
  ggplot() +
    geom_sf(data = map_layers$countries, fill = map_style$land, colour = map_style$border, linewidth = 0.3) +
    geom_sf(data = focus, aes(fill = I(country_tints[country])), colour = NA) +
    geom_sf(data = focus, aes(colour = I(country_colours[country])), fill = NA, linewidth = 0.45) +
    geom_sf_text(data = as_points(country_names), aes(label = name, colour = I(country_colours[country])), family = report_font, fontface = "bold", size = 3) +
    geom_sf_text(data = as_points(neighbours), aes(label = name), family = report_font, size = 2.1, colour = ink[["muted"]], lineheight = 0.9) +
    geom_sf_text(data = as_points(waters), aes(label = name), family = report_font, fontface = "italic", size = 2.1, colour = "#6A8796", lineheight = 0.9) +
    geom_sf(data = as_points(capitals), shape = 21, fill = "white", colour = ink[["primary"]], size = 1.4, stroke = 0.55) +
    geom_sf_text(data = as_points(capital_labels), aes(label = name, hjust = hj), family = report_font, size = 2.1, colour = ink[["secondary"]]) +
    geom_sf(data = centres, aes(fill = I(country_colours[country])), shape = 21, colour = "white", size = 4.6, stroke = 0.6) +
    geom_sf_text(data = centres, aes(label = num), family = report_font, fontface = "bold", size = 2.4, colour = "white") +
    coord_sf(crs = laea(42.5, 7), default_crs = 4326, xlim = c(32.4, 51.6), ylim = c(-2.0, 15.4), expand = FALSE, datum = NA) +
    scale_bar("bl") +
    theme_map()
}

# create the key from the map symbols
map_key <- function() {
  keys <- tibble(
    col = c(1, 1, 1, 1, 2, 2, 2), row = c(1, 2, 3, 4, 1, 2, 3),
    label = c(
      "Field base", "VBS venue town", "Field base and VBS venue", "Settlement footprint",
      "Study woreda (Ethiopia)", "Study district (Somalia)", "Main road / river"
    )
  )
  lx <- c(0.2, 5.2)
  keys$x <- lx[keys$col] + 0.45
  small_settlement <- map_style$settlement[["Small Settlement Area"]]
  ggplot() +
    annotate("point", x = lx[1], y = -1, shape = 22, fill = ink[["primary"]], colour = "white", size = 2.6, stroke = 0.55) +
    annotate("point", x = lx[1], y = -2, shape = 24, fill = method_colours[["vbs"]], colour = "white", size = 2.6, stroke = 0.55) +
    annotate("point", x = lx[1], y = -3, shape = 22, fill = method_colours[["vbs"]], colour = ink[["primary"]], size = 2.6, stroke = 0.65) +
    annotate("rect", xmin = lx[1] - 0.08, xmax = lx[1] + 0.12, ymin = -4.1, ymax = -3.9, fill = small_settlement) +
    annotate("rect", xmin = lx[1] - 0.22, xmax = lx[1] - 0.12, ymin = -3.98, ymax = -3.88, fill = small_settlement) +
    annotate("rect", xmin = lx[2] - 0.22, xmax = lx[2] + 0.22, ymin = -1.2, ymax = -0.8, fill = country_tints[["eth"]], colour = country_colours[["eth"]], linewidth = 0.5) +
    annotate("segment", x = lx[2] - 0.25, xend = lx[2] + 0.25, y = -2, yend = -2, colour = ink[["primary"]], linewidth = 0.45, linetype = "32") +
    annotate("segment", x = lx[2] - 0.25, xend = lx[2] + 0.02, y = -3, yend = -3, colour = map_style$main_road, linewidth = 0.45) +
    annotate("segment", x = lx[2] + 0.06, xend = lx[2] + 0.28, y = -3, yend = -3, colour = map_style$river, linewidth = 0.6) +
    geom_text(data = keys, aes(x = x, y = -row, label = label), hjust = 0, family = report_font, size = 2.4, colour = ink[["primary"]]) +
    scale_x_continuous(limits = c(0, 11), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-4.5, -0.5), expand = c(0, 0)) +
    theme_void()
}

# arrange the four study areas around the regional map and key
site_maps <- map(set_names(sites$panel), site_map)
figure_02 <- wrap_plots(
  site_maps$eth1, site_maps$eth2, site_maps$som1, site_maps$som2, region_map(), map_key(),
  design = c(area(1, 1, 6, 1), area(7, 1, 12, 1), area(1, 3, 6, 3), area(7, 3, 12, 3), area(1, 2, 9, 2), area(10, 2, 12, 2)),
  widths = c(4.3, 7.4, 4.3)
)
save_figure(figure_02, "output/figures/figure_02.png", height_cm = 11)

## Figure 6 ----
# Fieldwork periods, one calendar per phase, countries as row strips and
# methods as rows; each bar runs to the end of its last day.
# prepare the periods for plotting
timeline <-
  fieldwork |>
  mutate(
    country_f = factor(country_labels[country], levels = country_labels),
    method_f = factor(method_labels[method], levels = rev(method_labels)),
    panel = factor(phase_labels[as.character(phase)], levels = phase_labels),
    start_label = str_remove(format(start, "%d %b"), "^0"),
    end_label = str_remove(format(end, "%d %b"), "^0")
  )

# create the timeline of one phase
timeline_panel <- function(phase_label) {
  d <- filter(timeline, panel == phase_label)
  span <- as.numeric(max(d$end) - min(d$start))
  ggplot(d, aes(y = method_f)) +
    geom_segment(aes(x = start, xend = end + 1, yend = method_f, colour = I(method_colours[method])), linewidth = 3.6, lineend = "butt") +
    geom_text(aes(x = start - 0.012 * span - 0.8, label = start_label), hjust = 1, family = report_font, size = 2.2, colour = ink[["secondary"]]) +
    geom_text(aes(x = end + 0.012 * span + 1.8, label = end_label), hjust = 0, family = report_font, size = 2.2, colour = ink[["secondary"]]) +
    facet_grid(country_f ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", expand = expansion(mult = c(0.09, 0.09))) +
    labs(x = NULL, y = NULL, title = phase_label) +
    theme_report() +
    theme(
      strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      plot.title = element_text(face = "bold", size = 9, hjust = 0, colour = ink[["primary"]], margin = margin(0, 0, 4, 0)),
      plot.title.position = "plot", panel.grid.major.y = element_blank(), panel.spacing.y = unit(0.5, "lines")
    )
}

# stack the phases, panel heights following the number of rows
figure_06 <- wrap_plots(map(phase_labels, timeline_panel), ncol = 1, heights = as.numeric(table(timeline$panel)))
save_figure(figure_06, "output/figures/figure_06.png", height_cm = 11.5)

## Figures 7a and 7b ----
# Participant flow diagrams, one per phase: countries stacked down the page,
# methods as columns, the five stages down each column and the recorded
# losses in a side box to the right of the step they leave. One plot unit is
# one centimetre.
# lay out the boxes of one phase
flow_layout <- function(phase) {
  s <- filter(flow_stages, phase == .env$phase)
  l <- filter(flow_losses, phase == .env$phase)
  countries <- intersect(names(country_labels), s$country)
  methods <- intersect(method_order, s$method)
  left <- 2.6
  col_w <- (15.9 - left) / length(methods)
  main_w <- 0.33 * col_w
  side_w <- 0.575 * col_w
  pitch <- 1.2
  box_h <- 0.6
  head_h <- 0.95
  top_pad <- 0.15
  block_h <- head_h + pitch * 5 + 0.1
  height <- length(countries) * block_h + top_pad + 0.1
  tops <- set_names(height - top_pad - (seq_along(countries) - 1) * block_h, countries)
  columns <- expand_grid(country = countries, method = methods) |>
    mutate(x0 = left + (match(method, methods) - 1) * col_w + 0.12, xc = x0 + main_w / 2, y_top = tops[country])
  boxes <-
    s |>
    inner_join(columns, by = c("country", "method")) |>
    mutate(
      yc = y_top - head_h - box_h / 2 - (stage - 1) * pitch, label = format_count(n),
      xmin = xc - main_w / 2, xmax = xc + main_w / 2, ymin = yc - box_h / 2, ymax = yc + box_h / 2
    )
  sides <-
    l |>
    group_by(country, method, step) |>
    summarise(text = str_c(item, " ", format_count(n), collapse = "\n"), n_items = n(), .groups = "drop") |>
    inner_join(columns, by = c("country", "method")) |>
    mutate(
      ymid = y_top - head_h - box_h / 2 - (step - 0.5) * pitch, h = 0.3 * n_items + 0.16,
      xmin = xc + main_w / 2 + 0.2, xmax = xmin + side_w, ymin = ymid - h / 2, ymax = ymid + h / 2
    )
  list(
    boxes = boxes, sides = sides, heads = mutate(columns, x = x0, y = y_top - 0.42), height = height, left = left,
    blocks = tibble(country = countries, ymax = tops),
    stage_labels = expand_grid(country = countries, stage = 1:5) |>
      mutate(y = tops[country] - head_h - box_h / 2 - (stage - 1) * pitch, label = scales::label_wrap(18)(unique(flow_stages$stage_label))[stage])
  )
}

# create the flow diagram of one phase
flow_diagram <- function(phase) {
  g <- flow_layout(phase)
  boxes <- mutate(g$boxes, colour = method_colours[method])
  arrows <-
    boxes |>
    arrange(country, method, stage) |>
    group_by(country, method) |>
    reframe(x = head(xc, -1), y = head(ymin, -1), yend = tail(ymax, -1))
  heads <- mutate(g$heads, colour = method_colours[method], label = method_labels[method])
  country_heads <- g$blocks |> mutate(label = country_labels[country], y = heads$y[match(country, heads$country)])
  rules <- country_heads[-1, ]
  arrow_head <- grid::arrow(length = unit(1.2, "mm"), type = "closed", angle = 22)
  ggplot() +
    geom_segment(data = rules, aes(x = 0.1, xend = 15.9, y = ymax + 0.02, yend = ymax + 0.02), colour = ink[["grid"]], linewidth = 0.5) +
    geom_segment(data = arrows, aes(x = x, xend = x, y = y, yend = yend), colour = ink[["muted"]], linewidth = 0.35, arrow = arrow_head) +
    geom_segment(data = g$sides, aes(x = xc, xend = xmin, y = ymid, yend = ymid), colour = ink[["muted"]], linewidth = 0.35, arrow = arrow_head) +
    geom_rect(data = g$sides, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#F4F4F6", colour = NA) +
    geom_text(
      data = g$sides, aes(x = xmin + 0.1, y = ymid, label = text), hjust = 0, vjust = 0.5, size = 2.35, lineheight = 1.05,
      colour = ink[["secondary"]], family = report_font
    ) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = I(colour)), fill = "white", linewidth = 0.45) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(colour)), alpha = 0.10, colour = NA) +
    geom_text(data = boxes, aes(x = xc, y = yc, label = label), size = 2.6, fontface = "bold", colour = ink[["primary"]], family = report_font) +
    geom_point(data = heads, aes(x = x + 0.12, y = y, colour = I(colour)), size = 2.2) +
    geom_text(data = heads, aes(x = x + 0.38, y = y, label = label), hjust = 0, size = 3, fontface = "bold", colour = ink[["primary"]], family = report_font) +
    geom_text(
      data = g$stage_labels, aes(x = g$left - 0.12, y = y, label = label), hjust = 1, size = 2.4, lineheight = 0.95,
      colour = ink[["secondary"]], family = report_font
    ) +
    geom_text(data = country_heads, aes(x = 0.1, y = y, label = label), hjust = 0, size = 3.3, fontface = "bold", colour = ink[["primary"]], family = report_font) +
    coord_cartesian(xlim = c(0, 16), ylim = c(0, g$height), expand = FALSE) +
    theme_void(base_family = report_font) +
    theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "white", colour = NA))
}

# save one diagram per phase
walk(1:2, \(p) save_figure(flow_diagram(p), str_glue("output/figures/figure_07{letters[p]}.png"), height_cm = flow_layout(p)$height))

## Figures 8 and 9 ----
# Sample composition, one dot chart per country: a row per category, grouped
# by characteristic. Phase 1 marks the household survey with a diamond and
# Phase 2 the mean of the methods with a tick; a ring marks a category more
# than 10 pp from that reference. The legend gives each method's mean
# distance.
# categories shown and their labels
categories <- tibble(
  variable = c("respondent_sex", "age_band", "age_band", "age_band", "child_sex", "residence", "education", "education", "education", "education", "wealth_grp", "religion", "religion"),
  category = c("Female", "18-24", "25-34", "35+", "Female", "Urban", "No formal education", "Primary", "Secondary", "Higher", "High asset", "Muslim", "Christian"),
  group = c("Respondent sex", "Age band", "Age band", "Age band", "Child sex", "Residence", "Education", "Education", "Education", "Education", "Wealth", "Religion", "Religion"),
  short = c("Female", "18–24", "25–34", "35+", "Female", "Urban", "None", "Primary", "Secondary", "Higher", "High asset", "Muslim", "Christian")
) |>
  mutate(label = str_c(variable, ": ", category))

# compare each category with its reference
composition_points <- function(phase) {
  rows <-
    composition |>
    filter(phase == .env$phase) |>
    inner_join(categories, by = c("variable", "category"))
  reference <- if (phase == 1) {
    rows |> filter(method == "hh") |> select(country, label, reference_pct = pct)
  } else {
    rows |> group_by(country, label) |> summarise(reference_pct = mean(pct), .groups = "drop")
  }
  rows |>
    left_join(reference, by = c("country", "label")) |>
    mutate(
      flagged = method != "hh" & !is.na(reference_pct) & abs(pct - reference_pct) > 10,
      label = factor(label, levels = rev(categories$label)),
      group = factor(group, levels = unique(categories$group)),
      method = factor(method, levels = method_order)
    )
}

# create the dot chart of one country
composition_panel <- function(rows, country, phase, show_y) {
  r <- filter(rows, country == .env$country) |> mutate(country_label = country_labels[[.env$country]])
  distance <- filter(mean_distance, country == .env$country, phase == .env$phase)
  methods <- intersect(method_order, r$method)
  legend_labels <- map_chr(methods, \(m) {
    value <- distance$tvd[distance$method == m]
    if (m == "hh" || length(value) == 0) method_labels[[m]] else sprintf("%s (%.3f)", method_labels[[m]], value)
  })
  others <- filter(r, method != "hh")
  all_labels <- mutate(categories, label = factor(label, levels = levels(r$label)), group = factor(group, levels = levels(r$group)), country_label = r$country_label[[1]])
  p <- ggplot(r, aes(x = pct, y = label)) +
    # keep a row for every category, recorded or not
    geom_blank(data = all_labels, aes(x = 50, y = label), inherit.aes = FALSE) +
    geom_point(data = others, aes(fill = method), shape = 21, colour = "white", size = 2.8, stroke = 0.6) +
    geom_point(data = filter(others, flagged), shape = 1, size = 4.3, stroke = 0.6, colour = ink[["primary"]])
  if (phase == 1) {
    p <- p + geom_point(data = filter(r, method == "hh"), aes(fill = method), shape = 23, colour = "white", size = 3.6, stroke = 0.6)
  } else {
    reference <- distinct(r, label, group, country_label, reference_pct) |> mutate(key = "ref")
    p <- p + geom_point(data = reference, aes(x = reference_pct, y = label, fill = key), shape = 124, size = 4.2, colour = ink[["secondary"]], inherit.aes = FALSE)
  }
  keys <- c(methods, if (phase == 2) "ref")
  key_shapes <- c(if (phase == 1) c(23, rep(21, length(methods) - 1)) else rep(21, length(methods)), if (phase == 2) 124)
  key_colours <- c(rep("white", length(methods)), if (phase == 2) ink[["secondary"]])
  p <- p +
    scale_fill_manual(
      values = c(method_colours[methods], ref = ink[["secondary"]]), breaks = keys, labels = c(legend_labels, if (phase == 2) "Mean of methods"),
      name = "Method (mean TVD)",
      guide = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", override.aes = list(shape = key_shapes, colour = key_colours, size = 2.6))
    ) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0.02, 0)) +
    scale_y_discrete(labels = set_names(all_labels$short, as.character(all_labels$label))) +
    facet_grid(group ~ country_label, scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Respondents in category (%)", y = NULL) +
    theme_report(base_size = 8.5) +
    theme(
      legend.position = "bottom", legend.title = element_text(face = "bold", size = 8),
      legend.key.spacing.y = unit(0, "pt"), legend.key.height = unit(0.35, "cm"),
      legend.text = element_text(size = 8), panel.grid.major.x = element_line(colour = ink[["grid"]], linewidth = 0.3),
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.text.x = element_text(face = "bold", size = 9.5)
    )
  # mark categories a country did not collect
  missing <- filter(all_labels, !label %in% r$label)
  if (nrow(missing) > 0) {
    p <- p + geom_text(data = missing, aes(x = 50, y = label, label = "not collected"), inherit.aes = FALSE, size = 2.4, colour = ink[["muted"]], family = report_font, fontface = "italic")
  }
  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(), strip.text.y.left = element_blank(), strip.background.y = element_blank())
  }
  p + theme(panel.widths = unit(5.5, "cm"))
}

# place the two countries side by side
composition_figure <- function(phase) {
  rows <- composition_points(phase)
  wrap_plots(composition_panel(rows, "eth", phase, TRUE), composition_panel(rows, "som", phase, FALSE), nrow = 1)
}
save_figure(composition_figure(1), "output/figures/figure_08.png", height_cm = 12.5)
save_figure(composition_figure(2), "output/figures/figure_09.png", height_cm = 12.5)
