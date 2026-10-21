# =============================================================================
# sensitivity.R
# Sensitivity of the During-period estimate and elasticity to the exempt
# share (s) and pass-through (p) assumptions used in the tax adjustment.
#
# NOTE: run after normaldid.R in the same R session. The "Period-level
# objects carried over" section uses objects created there (trunc_*, Y1, Y0,
# mean_estimates, placebo_est, placebo_mean, announced_mean, ...).
#
# Inputs:  Opportunity Insights spending data (.dta)
#          dropped_weeks.csv (from panelconstruct.R)
# Outputs: sensitivity table (LaTeX, printed to console)
# =============================================================================


# ---- Paths (edit these) -----------------------------------------------------

raw_data_path      <- "path/to/spending_data.dta"
dropped_weeks_path <- "path/to/dropped_weeks.csv"


# ---- One-time package setup (uncomment to install) --------------------------

# install.packages("haven")
# install.packages("tidyverse")
# install.packages("modelsummary")
# install.packages("fixest")
# install.packages("devtools")
# install.packages("https://github.com/synth-inference/synthdid/archive/refs/heads/master.tar.gz",
#                  repos = NULL, type = "source")


# ---- Packages ---------------------------------------------------------------

library(haven)
library(tidyverse)
library(modelsummary)
library(fixest)
library(devtools)
library(scinference)


# ---- Data -------------------------------------------------------------------

raw_data <- read_dta(raw_data_path)

dropped_weeks <- read.csv(dropped_weeks_path) |>
  mutate(week_end = as.Date(week_end)) |>
  distinct(week_end)


# ---- Data transformation ----------------------------------------------------

# Weekly (Sunday) observations for Tennessee (47) and donor states,
# Mar 2022 - Mar 2024, net of disaster-excluded weeks
us_data <- raw_data |>
  mutate(across(starts_with("spend_"),
                ~ as.numeric(na_if(as.character(.x), ".")))) |>
  mutate(date = make_date(year, month, day)) |>
  filter(date >= as.Date("2022-03-06"), date <= as.Date("2024-03-03")) |>
  filter(statefips %in% c(47, 5, 16, 28, 29, 30, 31, 38)) |>
  mutate(weekly = if_else(freq == "w", 1, 0)) |>
  mutate(weekly = if_else(wday(date) == 1, 1, 0)) |>
  filter(weekly == 1) |>
  anti_join(dropped_weeks, by = c("date" = "week_end"))

# Period indicators (same windows as normaldid.R)
siphoned_us_data <- us_data |>
  mutate(pre       = if_else((date >= as.Date("2022-03-06") &
                                as.Date("2023-02-12") > date), 1, 0)) |>
  mutate(announced = if_else((as.Date("2023-02-12") <= date &
                                as.Date("2023-04-23") > date), 1, 0)) |>
  mutate(passed    = if_else((as.Date("2023-04-23") <= date &
                                as.Date("2023-08-06") > date), 1, 0)) |>
  mutate(during    = if_else((as.Date("2023-08-06") <= date &
                                as.Date("2023-11-05") > date), 1, 0)) |>
  mutate(post      = if_else(as.Date("2023-11-05") <= date &
                               date <= as.Date("2024-03-03"), 1, 0)) |>
  mutate(tennessee = if_else(statefips == 47, 1, 0))

siphoned_us_data <- siphoned_us_data |>
  mutate(
    # share of each week's days inside the Aug 1 - Oct 31, 2023 holiday
    exposure = case_when(
      date == as.Date("2023-08-06") ~ 6/7,
      date == as.Date("2023-11-05") ~ 2/7,
      date >= as.Date("2023-08-13") & date <= as.Date("2023-10-29") ~ 1,
      TRUE ~ 0
    ),
    # only Tennessee's tax changed
    exposure = if_else(statefips == 47, exposure, 0)
  )

siphoned_us_data <- siphoned_us_data |>
  mutate(period = case_when(pre       == 1 ~ "Pre",
                            announced == 1 ~ "Announced",
                            passed    == 1 ~ "Passed",
                            during    == 1 ~ "During",
                            post      == 1 ~ "Post"),
         period = factor(period, levels = c("Pre", "Announced", "Passed", "During", "Post"))) |>
  filter(!is.na(period))

tau <- 0.065  # combined state + local rate on food removed in 2023


# ---- Inference functions ----------------------------------------------------

# Moving-block permutation p-value for H0: effect = theta0 over the T1 block
# (estimation_method and lsei_type are accepted but unused)
movingblock_signed <- function(Y1, Y0, T1, T0, theta0, estimation_method, lsei_type) {
  T01 <- T0 + T1

  Y1_0 <- Y1
  Y1_0[(T0 + 1):T01] <- Y1[(T0 + 1):T01] - theta0

  u.hat <- scinference:::did(Y1_0, Y0)$u.hat

  sub.size <- T1
  u.hat.c  <- c(u.hat, u.hat)
  S.vec    <- matrix(NA, T01, 1)
  for (s in 1:(T01)) {
    S.vec[s, 1] <- abs(sum(u.hat.c[s:(s + sub.size - 1)]))
  }
  p <- mean(S.vec >= S.vec[T0 + 1])
  return(p)
}

# Subset to Pre weeks plus the weeks of one period
# (takes wide, Y1, Y0 as arguments so it can be called inside run_scenario)
truncation <- function(block_name, wide, Y1, Y0) {
  pre_n        <- which(wide$period == "Pre")
  block_n      <- which(wide$period == block_name)
  complete_set <- c(pre_n, block_n)
  T0 <- length(pre_n)
  T1 <- length(block_n)
  list(Y1 = Y1[complete_set], Y0 = Y0[complete_set, ], T0 = T0, T1 = T1)
}

# Test inversion over ci_grid (estimation_method and lsei_type unused)
ci_binned <- function(Y1, Y0, T1, T0, estimation_method, alpha, ci_grid, lsei_type) {

  ps_temp <- rep(NA, length(ci_grid))
  for (ind in 1:length(ci_grid)) {
    ps_temp[ind] <- movingblock_signed(Y1, Y0, T1, T0, theta0 = ci_grid[ind])
  }
  ci_temp <- ci_grid[ps_temp > alpha]
  lb      <- min(ci_temp, na.rm = TRUE)
  ub      <- max(ci_temp, na.rm = TRUE)

  return(list(lb = lb, ub = ub))
}


# ---- Period-level objects carried over (requires normaldid.R) ---------------

announced1 <- movingblock_signed(trunc_announced$Y1, trunc_announced$Y0,
                                 trunc_announced$T1, trunc_announced$T0, theta0 = 0)
passed1    <- movingblock_signed(trunc_passed$Y1, trunc_passed$Y0,
                                 trunc_passed$T1, trunc_passed$T0, theta0 = 0)
post1      <- movingblock_signed(trunc_post$Y1, trunc_post$Y0,
                                 trunc_post$T1, trunc_post$T0, theta0 = 0)
preperiods <- movingblock_signed(Y1[1:35], Y0[1:35, ], T1 = 13, T0 = 22, theta0 = 0)

ci_binned_announced <- ci_binned(trunc_announced$Y1, trunc_announced$Y0,
                                 trunc_announced$T1, trunc_announced$T0,
                                 alpha = 0.05, ci_grid = seq(-0.3, 0.3, by = 0.001))
ci_binned_passed    <- ci_binned(trunc_passed$Y1, trunc_passed$Y0,
                                 trunc_passed$T1, trunc_passed$T0,
                                 alpha = 0.05, ci_grid = seq(-0.3, 0.3, by = 0.001))
ci_binned_post      <- ci_binned(trunc_post$Y1, trunc_post$Y0,
                                 trunc_post$T1, trunc_post$T0,
                                 alpha = 0.05, ci_grid = seq(-0.3, 0.3, by = 0.001))
ci_binned_pre       <- ci_binned(Y1[1:35], Y0[1:35, ], T1 = 13, T0 = 22,
                                 alpha = 0.05, ci_grid = seq(-0.3, 0.3, by = 0.001))

announced_est <- mean_estimates$mean[2]
passed_est    <- mean_estimates$mean[3]
post_est      <- mean_estimates$mean[5]

percent_placebo   <- 100 * placebo_est   / (1 + placebo_mean   - placebo_est)
percent_announced <- 100 * announced_est / (1 + announced_mean - announced_est)
percent_passed    <- 100 * passed_est    / (1 + passed_mean    - passed_est)
percent_during    <- 100 * during_est    / (1 + during_mean    - during_est)
percent_post      <- 100 * post_est      / (1 + post_mean      - post_est)


# ---- Scenario function ------------------------------------------------------

# Re-runs the During-period analysis for exempt share s and pass-through p
run_scenario <- function(s, p) {

  # Tax adjustment under (s, p)
  siphoned_us_data <- siphoned_us_data |>
    mutate(spend_grf_adj = (1 + spend_grf) * (1 + tau) /
             (1 + tau * (1 - (1 - p) * s * exposure)) - 1)

  # Event study
  event_study <- feols(spend_grf_adj ~ i(date, tennessee, ref = as.Date("2023-02-05")) |
                         statefips + date,
                       data = siphoned_us_data)

  during_mean <- siphoned_us_data |>
    filter(statefips == 47) |>
    filter(period == "During") |>
    summarize(mean = mean(spend_grf_adj)) |>
    pull(mean)

  # Level shift for re-centering
  tn_mean <- siphoned_us_data |>
    filter(statefips == 47) |>
    filter(date < as.Date("2023-02-12")) |>
    summarize(mean = mean(spend_grf_adj))

  control_mean <- siphoned_us_data |>
    filter(statefips != 47) |>
    filter(date < as.Date("2023-02-12")) |>
    summarize(mean = mean(spend_grf_adj))

  tn_last <- siphoned_us_data |>
    filter(statefips == 47) |>
    filter(date == as.Date("2023-02-05")) |>
    summarize(mean = mean(spend_grf_adj))

  control_last <- siphoned_us_data |>
    filter(statefips != 47) |>
    filter(date == as.Date("2023-02-05")) |>
    summarize(mean = mean(spend_grf_adj))

  dif <- tn_last - tn_mean - control_last + control_mean

  # Wide panel
  wide <- siphoned_us_data |>
    select(date, period, spend_grf_adj, statefips) |>
    pivot_wider(names_from = statefips, values_from = spend_grf_adj) |>
    arrange(date)

  Y1 <- wide[["47"]]

  Y0 <- wide |>
    select(-period, -date, -`47`) |>
    as.matrix()

  # Re-centered weekly estimates with period labels
  estimates_tbl <- tibble(term = names(coef(event_study)),
                          est  = coef(event_study)) |>
    mutate(date = as.Date(str_extract(term, "\\d{4}-\\d{2}-\\d{2}")))

  estimates_tbl <- estimates_tbl |>
    mutate(est = est + dif$mean) |>
    add_row(date = as.Date("2023-02-05"), est = dif$mean)

  estimates_tbl <- estimates_tbl |>
    left_join(wide |> select(date, period), by = "date")

  # During-period inference
  trunc_during <- truncation("During", wide, Y1, Y0)

  during1 <- movingblock_signed(trunc_during$Y1, trunc_during$Y0,
                                trunc_during$T1, trunc_during$T0, theta0 = 0)
  ci_binned_during <- ci_binned(trunc_during$Y1, trunc_during$Y0,
                                trunc_during$T1, trunc_during$T0,
                                alpha = 0.05, ci_grid = seq(-0.3, 0.3, by = 0.001))

  post_mean <- siphoned_us_data |>
    filter(statefips == 47) |>
    filter(period == "Post") |>
    summarize(mean = mean(spend_grf_adj)) |>
    pull(mean)

  # Point estimate, percent effect, elasticity
  during_est     <- mean(estimates_tbl$est[estimates_tbl$period == "During"])
  percent_during <- 100 * during_est / (1 + during_mean - during_est)
  elasticity     <- -(percent_during / 100) / ((1 - p) * s * tau / (1 + tau))

  tibble(s = s, p = p, est = during_est, pct_cf = percent_during,
         elasticity = elasticity, p_value = during1,
         lb = ci_binned_during$lb, ub = ci_binned_during$ub)
}


# ---- Sensitivity scenarios --------------------------------------------------

scenarios <- tribble(
  ~panel,             ~s,    ~p,
  "A. No adjustment", 0,     0,
  "B. Exempt share",  0.40,  0,
  "B. Exempt share",  0.50,  0,
  "B. Exempt share",  0.60,  0,
  "B. Exempt share",  0.70,  0,
  "B. Exempt share",  0.80,  0,
  "C. Pass-through",  0.60, -0.29,
  "C. Pass-through",  0.60,  0,
  "C. Pass-through",  0.60,  0.20
)

sens <- scenarios |>
  mutate(res = pmap(list(s, p), run_scenario)) |>
  select(-s, -p) |>
  unnest(res)

sens

datasummary_df(sens, output = "latex", fmt = 3)
