# =============================================================================
# normaldid.R
# Main analysis: Tennessee grocery tax holiday (Aug 1 - Oct 31, 2023).
# Builds the weekly state panel, adjusts spending for the tax change,
# estimates the event study, runs conformal (scinference) inference week by
# week and by period, and produces the summary, event-study, and results
# tables.
#
# Inputs:  Opportunity Insights spending data (.dta)
#          dropped_weeks.csv (from panelconstruct.R)
# Outputs: event_study.pdf, LaTeX tables printed to console
# =============================================================================


# ---- Paths (edit these) -----------------------------------------------------

raw_data_path      <- "path/to/spending_data.dta"
dropped_weeks_path <- "path/to/dropped_weeks.csv"
output_dir         <- "path/to/output"


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

# Period indicators
#   Pre:       2022-03-06 to before 2023-02-12
#   Announced: 2023-02-12 to before 2023-04-23
#   Passed:    2023-04-23 to before 2023-08-06
#   During:    2023-08-06 to before 2023-11-05
#   Post:      2023-11-05 to 2024-03-03
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


# ---- Group and period labels ------------------------------------------------

siphoned_us_data <- siphoned_us_data |>
  mutate(
    group  = factor(if_else(statefips == 47, "Tennessee", "Other states"),
                    levels = c("Tennessee", "Other states")),
    period = case_when(pre       == 1 ~ "Pre",
                       announced == 1 ~ "Announced",
                       passed    == 1 ~ "Passed",
                       during    == 1 ~ "During",
                       post      == 1 ~ "Post"),
    period = factor(period, levels = c("Pre", "Announced", "Passed", "During", "Post"))
  ) |>
  filter(!is.na(period))


# ---- Tax adjustment ---------------------------------------------------------

tau <- 0.065  # combined state + local rate on food removed in 2023
s   <- 0.6    # exempt share of grocery-store spending

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
    exposure = if_else(statefips == 47, exposure, 0),
    # convert tax-inclusive index to a pre-tax-equivalent index
    spend_grf_adj = (1 + spend_grf) * (1 + tau) /
      (1 + tau * (1 - s * exposure)) - 1
  )


# ---- Summary statistics -----------------------------------------------------

datasummary(period ~ (spend_grf + spend_grf_adj) * group * (Mean + SD + N),
            data = siphoned_us_data, output = "latex")

# Diagnostic plot: Tennessee grocery spending, Apr - Aug 2022
raw_data |>
  mutate(spend_grf = as.numeric(na_if(as.character(spend_grf), ".")),
         date      = make_date(year, month, day)) |>
  filter(statefips == 47, wday(date) == 1,
         date >= as.Date("2022-04-01"), date <= as.Date("2022-08-31")) |>
  ggplot(aes(date, spend_grf)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = as.Date("2022-06-12"), linetype = "dashed")


# ---- DiD / event-study estimation -------------------------------------------

siphoned_us_data <- siphoned_us_data |>
  mutate(treats = if_else(date < as.Date("2023-02-12"), "pre", as.character(date)))

event_study <- feols(spend_grf_adj ~ i(date, tennessee, ref = as.Date("2023-02-05")) |
                       statefips + date,
                     data = siphoned_us_data)

# Level shift: (TN last pre-week - TN pre mean) - (control last pre-week - control pre mean)
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


# ---- Weekly conformal inference ---------------------------------------------

# Wide panel: one row per week, one column per state
wide <- siphoned_us_data |>
  select(date, period, spend_grf_adj, statefips) |>
  pivot_wider(names_from = statefips, values_from = spend_grf_adj) |>
  arrange(date)

Y1 <- wide[["47"]]  # treated unit (Tennessee)

Y0 <- wide |>       # donor units
  select(-period, -date, -`47`) |>
  as.matrix()

# Joint run over all 46 post-announcement weeks (gives weekly CIs)
res <- scinference(
  Y1,
  Y0,
  T1 = 46,
  T0 = 35,
  alpha = 0.05,
  estimation_method = "did",
  ci = TRUE,
  ci_grid = seq(-0.3, 0.3, by = 0.001)
)

# One-week-at-a-time p-values: 35 pre-weeks plus post-week k
p_weekly = numeric(46)

for (k in 1:46) {
  rows <- c(1:35, 35 + k)
  out <- scinference(
    Y1[rows], Y0[rows, ], T0 = 35, T1 = 1, alpha = 0.05,
    estimation_method = "did",
  )
  p_weekly[k] <- out$p_val
}

weekly_inference <- tibble(
  date = wide$date[36:81],
  pval = p_weekly,
  lb   = res$lb,
  ub   = res$ub
)


# ---- Event-study figure -----------------------------------------------------

estimates_tbl <- tibble(term = names(coef(event_study)),
                        est  = coef(event_study)) |>
  mutate(date = as.Date(str_extract(term, "\\d{4}-\\d{2}-\\d{2}")))

# Re-center on the pre-period mean and add back the reference week
estimates_tbl <- estimates_tbl |>
  mutate(est = est + dif$mean) |>
  add_row(date = as.Date("2023-02-05"), est = dif$mean)

weekly_estimates <- estimates_tbl |>
  left_join(weekly_inference, by = "date")

fig_event <- ggplot(weekly_estimates, mapping = aes(x = date, y = est)) +
  # shading: During (darkest), Passed, Announced (lightest)
  annotate("rect", xmin = as.Date("2023-08-06"), xmax = as.Date("2023-10-29"),
           ymin = -Inf, ymax = Inf, alpha = 0.45) +
  annotate("rect", xmin = as.Date("2023-04-23"), xmax = as.Date("2023-08-06"),
           ymin = -Inf, ymax = Inf, alpha = 0.3) +
  annotate("rect", xmin = as.Date("2023-02-12"), xmax = as.Date("2023-04-23"),
           ymin = -Inf, ymax = Inf, alpha = 0.15) +
  scale_y_continuous(labels = scales::label_number(scale = 100)) +
  scale_x_date(date_breaks = "6 months") +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 6, linewidth = 0.4, na.rm = TRUE) +
  geom_point() +
  labs(x = NULL, y = "Effect on adjusted grocery spending (pp of 2020 level)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "event_study.pdf"), plot = fig_event,
       width = 6.5, height = 4.5, units = "in")


# ---- Binned (by-period) inference -------------------------------------------

# Subset to Pre weeks plus the weeks of one period
truncation <- function(block_name) {
  pre_n        <- which(wide$period == "Pre")
  block_n      <- which(wide$period == block_name)
  complete_set <- c(pre_n, block_n)
  T0 <- length(pre_n)
  T1 <- length(block_n)
  list(Y1 = Y1[complete_set], Y0 = Y0[complete_set, ], T0 = T0, T1 = T1)
}

trunc_during    <- truncation("During")
trunc_announced <- truncation("Announced")
trunc_passed    <- truncation("Passed")
trunc_post      <- truncation("Post")
trunc_pre       <- truncation("Pre")

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

announced1 <- movingblock_signed(trunc_announced$Y1, trunc_announced$Y0,
                                 trunc_announced$T1, trunc_announced$T0, theta0 = 0)
during1    <- movingblock_signed(trunc_during$Y1, trunc_during$Y0,
                                 trunc_during$T1, trunc_during$T0, theta0 = 0)
passed1    <- movingblock_signed(trunc_passed$Y1, trunc_passed$Y0,
                                 trunc_passed$T1, trunc_passed$T0, theta0 = 0)
post1      <- movingblock_signed(trunc_post$Y1, trunc_post$Y0,
                                 trunc_post$T1, trunc_post$T0, theta0 = 0)

# Placebo: last 13 pre-weeks vs first 22 pre-weeks
preperiods <- movingblock_signed(Y1[1:35], Y0[1:35, ], T1 = 13, T0 = 22, theta0 = 0)


# ---- Binned confidence intervals --------------------------------------------

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

ci_binned_during    <- ci_binned(trunc_during$Y1, trunc_during$Y0,
                                 trunc_during$T1, trunc_during$T0,
                                 alpha = 0.05, ci_grid = seq(-0.3, 0.3, by = 0.001))
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


# ---- Period-average estimates -----------------------------------------------

mini_wide <- wide |>
  select(date, period)

newweekly_estimates <- weekly_estimates |>
  left_join(mini_wide, by = "date")

mean_estimates <- newweekly_estimates |>
  group_by(period) |>
  summarize(mean = mean(est))

small_weekly <- weekly_estimates |>
  filter(date < as.Date("2023-02-12")) |>
  arrange(date)

placebo_est <- mean(tail(small_weekly$est, 13)) - mean(head(small_weekly$est, 22))


# ---- Percent effects and elasticity -----------------------------------------

# mean_estimates rows: 1 Pre, 2 Announced, 3 Passed, 4 During, 5 Post
announced_est <- mean_estimates$mean[2]
passed_est    <- mean_estimates$mean[3]
during_est    <- mean_estimates$mean[4]
post_est      <- mean_estimates$mean[5]

# Tennessee period means of adjusted spending
announced_mean <- siphoned_us_data |>
  filter(statefips == 47) |>
  filter(period == "Announced") |>
  summarize(mean = mean(spend_grf_adj)) |>
  pull(mean)

passed_mean <- siphoned_us_data |>
  filter(statefips == 47) |>
  filter(period == "Passed") |>
  summarize(mean = mean(spend_grf_adj)) |>
  pull(mean)

during_mean <- siphoned_us_data |>
  filter(statefips == 47) |>
  filter(period == "During") |>
  summarize(mean = mean(spend_grf_adj)) |>
  pull(mean)

post_mean <- siphoned_us_data |>
  filter(statefips == 47) |>
  filter(period == "Post") |>
  summarize(mean = mean(spend_grf_adj)) |>
  pull(mean)

placebo_mean <- siphoned_us_data |>
  filter(statefips == 47, period == "Pre") |>
  arrange(date) |>
  slice_tail(n = 13) |>
  summarize(mean = mean(spend_grf_adj)) |>
  pull(mean)

# Effect as a percent of the counterfactual level
percent_placebo   <- 100 * placebo_est   / (1 + placebo_mean   - placebo_est)
percent_announced <- 100 * announced_est / (1 + announced_mean - announced_est)
percent_passed    <- 100 * passed_est    / (1 + passed_mean    - passed_est)
percent_during    <- 100 * during_est    / (1 + during_mean    - during_est)
percent_post      <- 100 * post_est      / (1 + post_mean      - post_est)

percent_estimates <- c(percent_placebo, percent_announced, percent_passed,
                       percent_during, percent_post)

elasticity <- -(percent_during / 100) / (s * tau / (1 + tau))


# ---- Results table ----------------------------------------------------------

results <- tibble(
  block   = c("Placebo", "Announced", "Passed", "During", "Post"),
  est     = c(placebo_est, mean_estimates$mean[2:5]),
  percent = percent_estimates,
  weeks   = c(as.character(13),
              as.character(trunc_announced$T1),
              as.character(trunc_passed$T1),
              as.character(trunc_during$T1),
              as.character(trunc_post$T1)),
  p_value = c(preperiods, announced1, passed1, during1, post1),
  lb      = c(ci_binned_pre$lb, ci_binned_announced$lb, ci_binned_passed$lb,
              ci_binned_during$lb, ci_binned_post$lb),
  ub      = c(ci_binned_pre$ub, ci_binned_announced$ub, ci_binned_passed$ub,
              ci_binned_during$ub, ci_binned_post$ub)
)

results <- results |>
  mutate(Elasticity = if_else(block == "During", elasticity, NA_real_))

results <- results |>
  mutate(CI = sprintf("[%.3f, %.3f]", lb, ub)) |>
  select(-lb, -ub)

results <- results |>
  rename(Period    = block,
         Estimate  = est,
         Weeks     = weeks,
         "p-value" = p_value,
         "95% CI"  = CI)

datasummary_df(results, output = "latex",
               title = "Effects of the Tennessee grocery tax holiday by period",
               fmt = 3)
