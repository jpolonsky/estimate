# estimate

Analysis code for ESTIMATE, a study of whether alternative survey methods can validly estimate childhood vaccination coverage (children aged 12 to 23 months) in conflict-affected and inaccessible populations. Community-based enumeration (CBE), respondent-driven sampling (RDS) and venue-based sampling (VBS) were compared with a household survey (HH) in accessible areas of Ethiopia and Somalia (Phase 1), then conducted without a household survey in inaccessible areas (Phase 2), in 2025 and 2026 (University of Geneva, Addis Ababa University and SORDI; funded by the Gates Foundation).

The scripts reproduce the computed tables and figures of the final report and its annexes, in the order they appear:

| Script | Reproduces |
|---|---|
| `code/01_clean.R` | Clean analysis files in `data/clean/` |
| `code/02_estimates.R` | Coverage estimates for every sample (`data/clean/estimates.rds`) |
| `code/03_samples.R` | Figure 2 (study areas), Figure 6 (fieldwork periods), Table 7 (recruitment), Figures 7a and 7b (participant flow), Figures 8 and 9 (sample composition); Annex Table D1 (cluster design), Annex Tables F2 to F5 (exclusions and sample composition) |
| `code/04_coverage.R` | Table 8, Figure 11 (Phase 1 coverage); Annex Tables G1 to G10 (missingness, clusters, subgroups, places, dropout), Annex Figures G1 to G5 (all antigens, timeliness, age at dose, cluster coverage, dropout) |
| `code/05_validity.R` | Table 9, Figure 10 (Phase 1 scorecard), Figure 12 (difference from HH in both countries); Annex Table G11 (paired bootstrap of RMSE differences) |
| `code/06_phase2.R` | Table 10, Figure 13 (Phase 2 scorecard), Figure 14 (Phase 2 coverage), Figure 15 (pairwise differences), Table 12 (predicted against observed differences); Annex Tables J1 and J2 (pairwise differences, RDS record sets) |
| `code/07_rds_vbs.R` | Annex Tables E1 to E5 (VBS travel and weights, design effects, RDS recruitment and homophily), Annex Figures E1 to E4 (VBS venues, RDS convergence, recruitment chains, coverage within chains) |
| `code/08_sensitivity.R` | Table 11 (unknown dose counts, Somalia Phase 2 RDS), Table 13 (sensitivity analyses); Annex Tables H1 and H3 to H6 (card documentation, photo-corrected coverage, unknown doses, RDS calibration, mobility weights), Annex Tables I1 to I3 (adjusted estimates) |
| `code/09_network.R` | Figure 16 (network module against direct report); Annex Tables K1 to K8 (network module) |
| `code/10_stopping.R` | Stopping points and estimate paths (`data/clean/stopping.rds`), reported by `07_rds_vbs.R` and `11_cost.R` |
| `code/11_cost.R` | Table 14 (resource requirements), Figure 17 (cost by category), Figure 18 (ratios to HH), Figure 21 (stopping-point trajectory); Annex Tables L1 to L9 (cost detail), Annex Figure L1 (cost by margin of error) |
| `code/12_field_effort.R` | Figure 19 (team hours), Figure 20 (interview duration); Annex Tables M1 to M4 (staff, coordinator hours, hours per respondent, network module time), Annex Figure M1 (duration by section) |

`code/functions.R` holds the helpers shared by the scripts. Tables are written to `output/tables/` as UTF-8 CSV files laid out and rounded as printed in the report (`table_08.csv`, `table_G04.csv`), and figures to `output/figures/` as PNG files at the report's size (`figure_11.png`, `figure_E2.png`). Main-report scripts are numbered in report order; `07_rds_vbs.R` reproduces annex items only.

Items not reproduced:

| Item | Reason |
|---|---|
| Tables 1 to 6 | Table 1 is the hand-assembled verdict grid; Tables 2 to 6 describe the study design and methods (text, not computed) |
| Table 15 | Realised qualitative sample; the qualitative study's records are not shared |
| Table 16 | Estimates from external published surveys |
| Table 17 | Recommended approach by circumstance, hand-assembled from the findings |
| Figures 1, 5 and 22 | Conceptual diagrams (study design, cost schema, decision tree) |
| Figure 3 | External image |
| Figure 4 | Needs household GPS coordinates, which are not in the shared data |
| Annex Tables A1, B1 to B3, C1 | Text tables: analysis list, form modules, tally and time-motion sheets, indicator coding |
| Annex Table C2 | R package list; `renv.lock` records the packages |
| Annex Table F1 | Hand-assembled from gazetteer place names |
| Annex Table H2 | Comes from a separate card-photograph audit, whose data are not shared |
| Annex Table N1 | Qualitative participants |
| Annex Table O1, Annex Figures O1 to O4 | Ethics approvals and scanned approval letters |

## Running

Run the scripts in order from the repository root with R 4.5 or later, after placing the input files in `data/raw/` (see below), each in its own R session. `01_clean.R` and `02_estimates.R` come first; after that, each of `03` to `12` runs on its own from `data/clean/`, except that `07_rds_vbs.R` and `11_cost.R` also need `data/clean/stopping.rds` from `10_stopping.R`.

```
Rscript code/01_clean.R
Rscript code/02_estimates.R
Rscript code/10_stopping.R
Rscript code/03_samples.R
Rscript code/04_coverage.R
Rscript code/05_validity.R
Rscript code/06_phase2.R
Rscript code/07_rds_vbs.R
Rscript code/08_sensitivity.R
Rscript code/09_network.R
Rscript code/11_cost.R
Rscript code/12_field_effort.R
```

Most scripts take seconds to a minute. Two use every core but one through `parallel` (PSOCK sessions): `08_sensitivity.R` takes about 12 minutes on 15 cores and about an hour on one (the raking bootstraps, multiple imputation and a 2,000-replicate recruitment-tree bootstrap for each Phase 1 RDS sample), and `10_stopping.R` about 30 minutes on 15 cores and several hours on one (it refits each method's estimator for Penta3, OPV3 and MCV1 along 100 to 200 random enrolment orders per Phase 1 sample). Every random step is seeded, so the outputs are the same on every run.

## Software

The analysis was run with R 4.5.2 on macOS. Packages are loaded with `pacman::p_load()` at the top of each script; `renv.lock` records the version of every package used and its dependencies (`renv::restore()` installs them). The main packages are `tidyverse`, `survey`, `RDS`, `mice`, `sf`, `patchwork`, `ggtext`, `ggh4x`, `ggrepel`, `ggspatial`, `ragg` and `systemfonts`.

The figures use Aptos, the report font, which ships with Microsoft Office: `functions.R` registers it from the Office applications' font folder on macOS when it is there, and uses Arial otherwise. With Arial the figures are otherwise the same, but text widths and positions differ slightly from the report.

## Inputs

`data/raw/` (not included; see Data access) holds the de-identified raw files:

| File | Content |
|---|---|
| `<country>_phase<n>_<method>.rds` | One per KoboToolbox form, 14 in all (`eth` and `som`; Phase 1 `hh`, `cbe`, `rds`, `vbs`; Phase 2 `cbe`, `rds`, `vbs`). HH, RDS and VBS forms are a list of three tables: `main` (one row per respondent), `child` (one row per child reported directly) and `network` (one row per child reported through the network module). CBE forms are one table, one row per child; CBE was collected on paper and keyed in later, and its `date` is the recovered field day. Dates of birth and dose dates are not included: each child has `age_months` (age in completed months at interview) and, for each dose slot (`opv_d1` to `opv_d4`, `penta_d1` to `penta_d3`, `measles_d1`, `measles_d2`), `<slot>_age_days` (age at that dose in days, to the day in CBE and Phase 2 and to the month in the Phase 1 HH, RDS and VBS forms) and `<slot>_dated` (a date was recorded). The interview date is kept, since the fieldwork timeline and the enrolment order use it. |
| `som_phase2_rds_all_enumerators.rds` | The Somalia Phase 2 RDS form including the two enumerators excluded from the primary analysis |
| `som_phase2_rds_tracking.rds` | Somalia Phase 2 RDS supervisor tracking sheet: participant code, referral code and date |
| `cluster_weights.rds` | HH and CBE cluster sampling frames, with the clusters fielded |
| `tallies.rds` | Supervisor tally sheets: absent, ineligible, refusing and completed households per cluster |
| `card_photo_status.rds` | Vaccination status read from the card photographs of Phase 1 HH children |
| `time_motion.rds` | Time-motion segments, supervisor tracker hours and coordinator timesheets (staff as role and number; the UNIGE research team's hours summed over its members) |
| `cost_line_items.rds` | Budget line items in the currency each partner budgeted in. Lines paying named roles (research team, supervisors, field guides and facilitators, enumerators) are summed into one personnel line per partner, role group, country, phase, method, cost category, variability and scope; these group lines carry no allocation basis or activity phase |
| `cost_shared_items.rds` | The budget lines AAU and SORDI booked once for a country and phase, which each method's line items repeat |
| `research_team_cost.rds` | The UNIGE research team's Phase 1 timesheet hours costed at the team's day rates (US dollars at 2025 prices, before overhead), summed over the team by country and method |
| `fx_cpi.rds` | Exchange rates and US CPI-U values used to convert costs to constant 2026 US dollars, with their sources |
| `cluster_places.rds` | The named place (public gazetteer name) of each HH and CBE design cluster |
| `kobo_record_counts.rds` | Records in each KoboToolbox export before cleaning |

Two Ethiopia Phase 2 VBS venues whose names reveal an institutional setting are released as "Institutional venue 1" and "Institutional venue 2"; no output names them. Corrections to the raw records (data-entry errors, deletions, referral repairs), the assignment of household-survey clusters from GPS points and the removal of identifiers were applied before release, so the public code contains no record-level corrections. Record, respondent, enumerator and staff identifiers are integer codes.

`data/geo/` (included) holds the map layers of Figure 2; their sources and licences are in `data/geo/SOURCES.md`. The report drew the Ethiopian woreda boundaries from FEWS NET; this release uses OCHA's openly licensed boundaries instead, so woreda lines in Figure 2 may differ slightly from the printed report.

## Clean files

`01_clean.R`, `02_estimates.R` and `10_stopping.R` write these files to `data/clean/`:

| File | Content and key variables |
|---|---|
| `frames.rds` | List of 15 data frames, one per sample (plus `som_phase2_rds_all_enumerators`), one row per child reported directly (or per respondent without an eligible child). `record_ix` (respondent; `record_id` for CBE), `cluster_id` (HH) or `cluster` (CBE), `date`, `consent` (Enrolled, Refused, Ineligible), `participant_status` (VBS), `analytic` (child in the analysis), `id_participant_calc`, `id_referring`, `seed_participant`, `nsum_child_total` (RDS), `int_location` and the `*_travel` answers (VBS), respondent characteristics, `age_months` and the age at each dose (`*_d<k>_age_days`, `*_d<k>_dated`), and the antigen indicators `penta_1p_clean`, `penta_3p_clean`, `opv_1p_clean`, `opv_3p_clean`, `measles_1p_clean`, `measles_2p_clean` (1, 0 or missing) with their source (`*_source`: card or recall) |
| `network.rds` | List of 10 data frames (HH, RDS and VBS samples), one row per child reported through the network module, with the same antigen indicators and `nsum_caregiver_knows` |
| `estimates.rds` | One row per sample and outcome (Penta1, Penta3, OPV1, OPV3, MCV1, MCV2, Zero-dose): `estimate`, `standard_error`, `ci_lower`, `ci_upper` (proportions), `n_children`, `n_observed`, `n_respondents`, `n_network_nodes` |
| `stopping.rds` | List: `stopping` (Penta3 stopping point per Phase 1 sample: `stopping_basis`, `stopping_n` and its 10th and 90th percentiles over orders), `realised` (estimate and 95% CI along the realised enrolment order: Penta3 in every Phase 1 sample, and OPV3, MCV1 and zero-dose in RDS) and `order_stopping` (stability, precision and stopping point of each random order for Penta3, OPV3 and MCV1) |
| `hours.rds` | One row per person, day (or week) and activity: `source`, `phase`, `country`, `method`, `role`, `person`, `activity`, `hours` |
| `cost_line_items.rds`, `cost_shared_items.rds`, `fx_cpi.rds` | As in `data/raw/`, with lower-case country and method codes |
| `research_team_cost.rds` | `country` (eth, som or both), `method` (or cross-cutting), `day_rate_cost_usd` |
| `cluster_places.rds` | Fielded clusters only: `country`, `phase`, `method`, `cluster_id`, `place_name`, `place_type` |
| `kobo_record_counts.rds` | `country`, `phase`, `method`, `n_records` |
| `cluster_weights.rds`, `tallies.rds`, `card_photo_status.rds`, `rds_tracking.rds` | As in `data/raw/` |

The clean questionnaires and data dictionaries of the HH, RDS and VBS forms (English, one per country and phase, with every field, its skip logic and answer options) are deposited with the data and are available with them on request.

## Data access

The data contain information on children's vaccination and on respondents' households, villages and social networks in conflict-affected areas, and cannot be deposited in a public repository. De-identified data are available on reasonable request from the corresponding author at the University of Geneva, Dr Jonathan Polonsky (<Jonathan.Polonsky@unige.ch>), subject to the approvals of the study's ethics committees.

## Citation

If you use this code, please cite the report (citation details will be added on publication) and the code itself:

Jonathan Polonsky. (2026). jpolonsky/estimate [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23673550

## Licence

The code is released under the MIT licence (see `LICENSE`). The map layers in `data/geo/` keep their own licences, listed in `data/geo/SOURCES.md`.
