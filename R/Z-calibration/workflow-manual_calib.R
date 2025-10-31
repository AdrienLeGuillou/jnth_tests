## HPC Workflow: Manual Calibration 1
##
## Define a workflow to proposal parameters for calibration. This runs the first
## part of the model (before the restart point). The values are assessed with
## the script 1-manual_calib_assess.R

# Restart R before running this script (Ctrl_Shift_F10 / Cmd_Shift_0)

# Settings ---------------------------------------------------------------------
library(slurmworkflow)
library(EpiModelHPC)
library(EpiModelHIV)
library(dplyr)

hpc_context <- TRUE
source("R/shared_variables.R", local = TRUE)
source("R/Z-calibration/z-context.R", local = TRUE)
source("R/hpc_configs.R", local = TRUE)

max_cores <- 8

# Process ----------------------------------------------------------------------
source("R/netsim_settings.R", local = TRUE)

# Control settings
control <- control_msm(
  nsteps = calibration_end,
  .tracker.list = EpiModelHIV::make_calibration_trackers()
)

# Workflow creation ------------------------------------------------------------
wf <- make_em_workflow("calibration_1", override = TRUE)

# Using scenarios --------------------------------------------------------------

# Define calibration scenarios
# insert test values here
n_scenarios <- 10
rrs <- c(0.5, 1, 2)
rrs1 <- rep(rrs, each = 3)
rrs2 <- rep(rrs, 3)
scenarios_df <- tibble(
  # .scenario.id = as.character(seq_len(n_scenarios)),
  .scenario.id = paste0("hiv_", rrs1, "__prep_", rrs2),
  .at = 1,
  substance.effects = 0,
  hiv.test.rate_1 = param$hiv.test.rate[1] * rrs1,
  hiv.test.rate_2 = param$hiv.test.rate[2] * rrs1,
  hiv.test.rate_3 = param$hiv.test.rate[3] * rrs1,
  prep.start.rate_1 = param$prep.start.rate[1] * rrs2,
  prep.start.rate_2 = param$prep.start.rate[2] * rrs2,
  prep.start.rate_3 = param$prep.start.rate[3] * rrs2,
)
scenarios_list <- EpiModel::create_scenario_list(scenarios_df)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_est, param, init, control,
    scenarios_list = scenarios_list,
    output_dir = calib_dir,
    n_rep = 128,
    n_cores = max_cores,
    max_array_size = 500,
    setup_lines = hpc_node_setup
  ),
  sbatch_opts = list(
    "mail-type" = "FAIL,TIME_LIMIT",
    "cpus-per-task" = max_cores,
    "time" = "04:00:00",
    "mem-per-cpu" = "5G"
  )
)

# Process calibrations
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_merge_netsim_scenarios_tibble(
      sim_dir = calib_dir,
      output_dir = fs::path(calib_dir, "merged_tibbles"),
      steps_to_keep = Inf, # keep everything
      cols = dplyr::everything(),
      n_cores = max_cores,
      setup_lines = hpc_node_setup
    ),
    sbatch_opts = list(
      "mail-type" = "END",
      "cpus-per-task" = max_cores,
      "time" = "02:00:00",
      "mem-per-cpu" = "5G"
    )
)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call_script(
    r_script = "R/Z-calibration/process_calibs.R",
    args = list(hpc_context = TRUE),
    setup_lines = hpc_node_setup
  ),
  sbatch_opts = list(
    "cpus-per-task" = max_cores,
    "time" = "02:00:00",
    "mem-per-cpu" = "5G"
  )
)