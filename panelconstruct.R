# =============================================================================
# panelconstruct.R
# FEMA disaster exclusions for the Tennessee grocery tax holiday panel.
# Applies Donor Exclusion 2 and the Week Exclusion from the panel-construction
# document.
#
# Inputs:  FEMA Disaster Declarations Summaries (CSV)
#          Census county population estimates, 2020-2023 (CSV)
# Outputs: dropped_states.csv, dropped_weeks.csv
# =============================================================================


# ---- Paths (edit these) -----------------------------------------------------

disaster_path   <- "path/to/disaster_declarations.csv"
population_path <- "path/to/county_population_estimates.csv"
output_dir      <- "path/to/output"


# ---- Packages ---------------------------------------------------------------

library(httr)
library(jsonlite)
library(dplyr)
library(lubridate)
library(tidyverse)


# ---- Parameters -------------------------------------------------------------

states     <- c("TN", "AR", "FL", "ID", "MS", "MO", "MT", "NE", "ND")
sample_beg <- as.Date("2022-02-28")  # Monday of week ending 2022-03-06
sample_end <- as.Date("2024-03-03")  # last Sunday of Post
w_beg      <- as.Date("2023-02-06")  # Monday of week ending 2023-02-12 (Announced)
during_beg <- as.Date("2023-07-31")  # Monday of the week ending Aug 6
during_end <- as.Date("2023-10-29")
threshold  <- 0.15


# ---- Load and clean FEMA declarations ---------------------------------------

disaster_data <- read.csv(disaster_path)

# Major disaster declarations in sample states that overlap the sample window
disaster_data <- disaster_data |>
  filter(declarationType == "DR") |>
  filter(incidentType != "Biological") |>
  mutate(incidentBeginDate = as.Date(incidentBeginDate)) |>
  mutate(incidentEndDate   = as.Date(incidentEndDate)) |>
  filter(state %in% states) |>
  filter(incidentBeginDate <= sample_end & incidentEndDate >= sample_beg)

# Build 5-digit FIPS; keep county designations, plus statewide designations
# that carry an IH, IA, or PA program
disaster_data <- disaster_data |>
  mutate(
    fipsStateCode  = str_pad(fipsStateCode,  width = 2, pad = "0"),
    fipsCountyCode = str_pad(fipsCountyCode, width = 3, pad = "0"),
    fips           = paste0(fipsStateCode, fipsCountyCode)
  ) |>
  filter(
    fipsCountyCode != "000" |
      (designatedArea == "Statewide" &
         (ihProgramDeclared == 1 | iaProgramDeclared == 1 | paProgramDeclared == 1))
  ) |>
  distinct(disasterNumber, state, fips, fipsCountyCode,
           incidentBeginDate, incidentEndDate)


# ---- Link population data ---------------------------------------------------

population_data <- read.csv(population_path)

population_data <- population_data |>
  mutate(
    STATE  = str_pad(STATE,  width = 2, pad = "0"),
    COUNTY = str_pad(COUNTY, width = 3, pad = "0"),
    fips   = paste0(STATE, COUNTY)
  )

county_pop <- population_data |>
  filter(SUMLEV == 50) |>
  select(fips, county_pop = POPESTIMATE2022)

state_pop <- population_data |>
  filter(SUMLEV == 40) |>
  select(st = STATE, state_pop = POPESTIMATE2022)

# Population share of each state covered by each declaration
shares_data <- disaster_data |>
  left_join(county_pop, by = "fips") |>
  group_by(disasterNumber, state, incidentBeginDate, incidentEndDate) |>
  summarize(
    statewide = any(fipsCountyCode == "000"),
    decl_pop  = sum(county_pop, na.rm = TRUE),
    st        = first(substr(fips, 1, 2)),
    .groups   = "drop"
  ) |>
  left_join(state_pop, by = "st") |>
  mutate(share = if_else(statewide, 1, decl_pop / state_pop))


# ---- Qualifying declarations ------------------------------------------------

# Sunday that ends the Monday-start week containing d
week_ending <- function(d) floor_date(d, "week", week_start = 1) + days(6)

qualifying <- shares_data |>
  filter(share >= 0.15) |>
  mutate(
    in_W = incidentEndDate >= w_beg,
    in_H = incidentBeginDate <= during_end & incidentEndDate >= during_beg
  )


# ---- Exclusions -------------------------------------------------------------

# Donor Exclusion 2: drop donor states with a qualifying disaster during the holiday
dropped_states <- qualifying |>
  filter(state != "TN", in_H) |>
  select(state, disasterNumber, incidentBeginDate, incidentEndDate, share, in_H)

# Week Exclusion: drop every sample week touched by a qualifying disaster
# in the remaining states
dropped_weeks <- qualifying |>
  filter(!state %in% dropped_states$state) |>
  rowwise() |>
  mutate(week_end = list(seq(week_ending(incidentBeginDate),
                             week_ending(incidentEndDate), by = "week"))) |>
  ungroup() |>
  unnest(week_end) |>
  filter(week_end >= week_ending(sample_beg), week_end <= sample_end) |>
  distinct(week_end, state, disasterNumber)


# ---- Save -------------------------------------------------------------------

write.csv(dropped_states, file.path(output_dir, "dropped_states.csv"), row.names = FALSE)
write.csv(dropped_weeks,  file.path(output_dir, "dropped_weeks.csv"),  row.names = FALSE)
