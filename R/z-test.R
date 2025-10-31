# Scratchpad for interactive testing before integration in a script
#
#
source("R/shared_variables.R", local = TRUE)

library(dplyr)
library(EpiModelHIV)

context <- "hpc"
source("R/netsim_settings.R", local = TRUE)

d_calib <- readRDS(fs::path(calib_dir, "merged_tibbles", "df__empty_scenario.rds"))

d_outs <- EpiModelHIV::mutate_calibration_targets(d_calib) |>
  mutate(sim = as.integer(as.factor(paste0(batch_number, "_", sim)))) |>
  as.epi.data.frame()

races <- c("B", "H", "W")
calib_plot_infos <- list(
  cc.dx = list(
    names = paste0("cc.dx.", races),
    ylab = "Proportion",
    text_offset = 0.01,
    fmt_target = scales::percent_format(0.1)
  ),
  cc.linked1m = list(
    names = paste0("cc.linked1m.", races),
    ylab = "Proportion",
    text_offset = 0.005,
    fmt_target = scales::percent_format(0.1)
  ),
  cc.vsupp = list(
    names = paste0("cc.vsupp.", races),
    ylab = "Proportion",
    text_offset = 0.005,
    fmt_target = scales::percent_format(0.1)
  ),
  i.prev.dx = list(
    names = paste0("i.prev.dx.", races),
    ylab = "Proportion",
    text_offset = 0.01,
    fmt_target = scales::percent_format(0.1)
  ),
  ir100.sti = list(
    names = c("ir100.gc", "ir100.ct"),
    ylab = "Infection Rate per 100 PYAR",
    text_offset = 0.3,
    fmt_target = scales::number_format(0.1)
  ),
  cc.prep = list(
    names = paste0("cc.prep.", races),
    ylab = "Proportion",
    text_offset = 0.005,
    fmt_target = scales::percent_format(0.1)
  ),
  disease.mr100 = list(
    names = "disease.mr100",
    ylab = "Proportion",
    text_offset = 0.01,
    fmt_target = scales::percent_format(0.1)
  ),
  num = list(
    names = "num",
    ylab = "Population",
    text_offset = 500,
    fmt_target = scales::number_format(1)
  )
)

make_calib_plot <- function(d, plot_info) {
  targets <- EpiModelHIV::get_calibration_targets()
  targets["num"] <- 1e5
  colors <-  c("steelblue", "firebrick", "seagreen")
  text_pos <- max(d$time) - 500
  par(mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
  offset <- plot_info$text_offset
  cur_targs <- plot_info$names

  plot(
    d,
    xaxt = "n",
    y = cur_targs,
    legend = TRUE,
    ylab = plot_info$ylab,
    xlab = "Calibration Years"
  )
  axis(1, seq(0, max(d$time), 10 * year_steps),
       labels = seq(0, max(d$time), 10 * year_steps) / year_steps)

  x <- round(colMeans(tail(d_outs[, cur_targs], 52)), 3)
  abline(h = targets[cur_targs], col = colors, lty = 2)
  for (i in seq_along(plot_info$names)) {
    v <- plot_info$fmt_target(x[i])
    text(text_pos, targets[cur_targs[i]] + offset, v, col = colors[i])
  }
}

targets <- EpiModelHIV::get_calibration_targets()
races_names <- c("B", "H", "W")
races <- 1:3

med_iqr <- function(x, fmtr) {
  vs <- quantile(x, c(0.5, 0.25, 0.75)) |> fmtr()
  paste0(vs[1], " [", vs[2], "-",  vs[3], "]")
}

col_names <- paste0(c("gono", "chla"), " IR100")
row_names <- c(
  "Target Statistic: STI IR100",
  "Simulations: STI IR100 (med [IQR])",
  "Calibrated Transmission Urethral Probability (per act)"
)

fmtr <- scales::label_number(0.1)
tar_names <- c("ir100.gono", "ir100.chla")
tar <- targets[tar_names] |> fmtr()
p_names <- c("gono.uret.prob", "chla.uret.prob")
prms <- unlist(param[p_names]) |> scales::label_number(0.01)()

sim_vals <- d_outs |>
  filter(time > max(time) - year_steps) |>
  select(sim, all_of(tar_names)) |>
  group_by(sim) |>
  summarise(across(everything(), mean)) |>
  select(-sim) |>
  summarise(across(everything(), \(x) med_iqr(x, fmtr))) |>
  as.character()

tbl <- rbind(tar, sim_vals, prms)
rownames(tbl) <- row_names
colnames(tbl) <- col_names

knitr::kable(tbl)