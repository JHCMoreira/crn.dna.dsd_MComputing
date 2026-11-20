if (!isNamespace(topenv())) {
    library(pracma)   # for signal processing utilities
    library(xtable)   # for generating LaTeX tables
}


# ' Estimate time constants for capacitor voltage response
#'#' @param time                  numeric vector of time points
#' @param voltage_input         numeric vector of input voltage (excitation)
#' @param capacitor_voltage_output numeric vector of capacitor voltage output 
#' @return list with:
#'  - tau_63: time at which capacitor voltage reaches 63.2%
#'  - tau_95: time at which capacitor voltage reaches 95%
#'  - tau_d36.8: time at which capacitor voltage discharges to 36.8% of its peak value
#'  - t_threshold: time when input voltage first exceeds 0.01 V

estimate_tau_capacitor_voltage <- function(time, voltage_input, capacitor_voltage_output) {
  # Ensure inputs are the same length
  if (!all(length(time) == length(voltage_input),
           length(time) == length(capacitor_voltage_output))) {
    stop("time, voltage_input and capacitor_voltage_output must be the same length")
  }

  # Find first point where excitation exceeds 0.01 V
  i0 <- which(voltage_input > 0.01)[1]
  if (is.na(i0)) stop("voltage_input never exceeds 0.01 V")

  # Reference time at threshold crossing
  t0 <- time[i0]

  # Truncate vectors from i0 onward
  t_full <- time[i0:length(time)]
  Vin     <- voltage_input[i0:length(time)]
  Vc      <- capacitor_voltage_output[i0:length(time)]

  # Compute step levels
  V0 <- min(Vin)
  V1 <- max(Vin)
  ΔV <- V1 - V0

  # Helper: linear interpolation of crossing time
  get_t_at <- function(t, y, target) {
    idx <- which(diff(sign(y - target)) != 0)
    if (length(idx) == 0) return(NA_real_)
    i <- idx[1]
    t1 <- t[i]; t2 <- t[i+1]
    y1 <- y[i]; y2 <- y[i+1]
    return(t1 + (target - y1) * (t2 - t1) / (y2 - y1))
  }

  # Charge targets (63.2% and 95% of ΔV)
  target_63  <- V0 + 0.632 * ΔV
  target_95  <- V0 + 0.95  * ΔV

  # Find relative times of crossings within truncated data
  t_rel_63 <- get_t_at(t_full - t0, Vc, target_63)
  t_rel_95 <- get_t_at(t_full - t0, Vc, target_95)

  # Convert to absolute times
  t_abs_63 <- if (!is.na(t_rel_63)) t0 + t_rel_63 else NA
  t_abs_95 <- if (!is.na(t_rel_95)) t0 + t_rel_95 else NA

  # DISCHARGE: identify peak index in Vin
  i_peak <- which.max(Vin)
  t_peak <- t_full[i_peak]
  Vp     <- Vc[i_peak]
  Vfinal <- Vc[length(Vc)]
  ΔVd    <- Vp - Vfinal
  target_36.8 <- Vfinal + 0.368 * ΔVd

  # Post-peak segment
  t_post_rel <- (t_full - t_peak)[i_peak:length(t_full)]
  Vc_post    <- Vc[i_peak:length(Vc)]
  t_rel_36.8 <- get_t_at(t_post_rel, Vc_post, target_36.8)

  # Absolute discharge time
  t_abs_36.8 <- if (!is.na(t_rel_36.8)) t_peak + t_rel_36.8 else NA

  list(
    tau_63     = t_rel_63,
    tau_95     = t_rel_95,
    tau_d36.8  = t_rel_36.8,
    t_threshold = t0
  )
}

#' Estimate the fundamental AC period of two waveforms
#'
#' @param timing        numeric vector of time points
#' @param source_input  not used here (placeholder for your input signal)
#' @param model_vals    numeric vector of model signal values
#' @param sim_vals      numeric vector of simulation signal values
#' @param start_time    time after which to begin analysis (defaults to 10)
#' @return              list(period_model, period_sim)
estimate_ac_cycle <- function(timing, source_input, model_vals, sim_vals, start_time = 10) {
  # slice to times >= start_time
  idx <- which(timing >= start_time)
  t  <- timing[idx]
  m  <- model_vals[idx]
  s  <- sim_vals[idx]
  
  # helper: find upward mean‐crossings
  find_period <- function(t, y) {
    ym <- mean(y)
    # indices where we cross from below to above the mean
    cross_up <- which(diff(y < ym) == 1)
    times   <- t[cross_up]
    # periods between successive crossings
    per     <- diff(times)
    mean(per)
  }
  
  period_model <- find_period(t, m)
  period_sim   <- find_period(t, s)
  
  list(period_model = period_model,
       period_sim   = period_sim)
}


#' Compute phase shift between current and voltage for model & sim
#'
#' @param timing    numeric vector of time points
#' @param curr_vals numeric vector of current values
#' @param volt_vals numeric vector of voltage values
#' @param start_time time after which to begin analysis (defaults to 10)
#' @return          list(phase_model, phase_sim) in degrees
compute_phase_shift <- function(timing, 
                                model_current, model_voltage, 
                                sim_current,   sim_voltage, 
                                start_time = 10) {
  idx <- which(timing >= start_time)
  t   <- timing[idx]
  
  # helper: compute phase shift of y2 relative to y1
  # by comparing their first upward zero-crossing times
  phase_deg <- function(t, y1, y2) {
    th1 <- mean(y1) 
    th2 <- mean(y2)
    # find first up-crossing after start
    i1 <- which(diff(y1 < th1) == 1)[1]
    i2 <- which(diff(y2 < th2) == 1)[1]
    dt  <- t[i2] - t[i1]
    # estimate period from y1
    per <- 1e-9
    if (length(which(diff(y1 < th1) == 1)) >= 2) {
      times <- t[which(diff(y1 < th1) == 1)]
      per   <- mean(diff(times))
    }
    (dt / per) * 360
  }
  
  phase_model <- phase_deg(t, model_current[idx], model_voltage[idx])
  phase_sim   <- phase_deg(t, sim_current[idx],   sim_voltage[idx])
  
  list(phase_model = phase_model,
       phase_sim   = phase_sim)
}


#' Compare sim vs model for phase and amplitude
#'
#' @param timing        numeric vector of time points
#' @param model_current numeric vector of model current
#' @param model_voltage numeric vector of model voltage
#' @param sim_current   numeric vector of sim current
#' @param sim_voltage   numeric vector of sim voltage
#' @param start_time    time after which to begin analysis (defaults to 10)
#' @return              list with:
#'   - phase_diff_current:   phase(model_current vs sim_current)
#'   - phase_diff_voltage:   phase(model_voltage vs sim_voltage)
#'   - amp_diff_current:     max(sim_current)-max(model_current), and same for minima
#'   - amp_diff_voltage:     likewise for voltage
compare_model_sim <- function(timing, 
                              model_current, model_voltage, 
                              sim_current,   sim_voltage, 
                              start_time = 10) {
  idx <- which(timing >= start_time)
  t   <- timing[idx]
  
  # reuse phase shift helper (y2 relative to y1)
  phase_deg <- function(t, y1, y2) {
    th1 <- mean(y1); th2 <- mean(y2)
    i1  <- which(diff(y1 < th1) == 1)[1]
    i2  <- which(diff(y2 < th2) == 1)[1]
    dt  <- t[i2] - t[i1]
    per <- 1e-9
    ups <- which(diff(y1 < th1) == 1)
    if (length(ups) >= 2) per <- mean(diff(t[ups]))
    (dt / per) * 360
  }
  
  # phase diffs
  phase_diff_current <- phase_deg(t, model_current[idx], sim_current[idx])
  phase_diff_voltage <- phase_deg(t, model_voltage[idx], sim_voltage[idx])
  
  # amplitude diffs
  amp_diff_current <- c(
    max_diff = max(sim_current[idx]) - max(model_current[idx]),
    min_diff = min(sim_current[idx]) - min(model_current[idx])
  )
  amp_diff_voltage <- c(
    max_diff = max(sim_voltage[idx]) - max(model_voltage[idx]),
    min_diff = min(sim_voltage[idx]) - min(model_voltage[idx])
  )
  
  list(
    phase_diff_current = phase_diff_current,
    phase_diff_voltage = phase_diff_voltage,
    amp_diff_current   = amp_diff_current,
    amp_diff_voltage   = amp_diff_voltage
  )
}

analyze_transient_metrics <- function(timing,
                                      v_in,
                                      vc_model,
                                      vc_sim,
                                      t0,
                                      t1,
                                      resistance = NA_real_,
                                      inductance = NA_real_,
                                      capacitance = NA_real_) {

  idx <- which(timing >= t0 & timing <= t1)
  t_slice <- timing[idx]
  v_in_slice <- v_in[idx]
  v_model <- vc_model[idx]
  v_sim   <- vc_sim[idx]

  safe_tail_median <- function(x, fraction = 0.1) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    n_tail <- max(3, ceiling(length(x) * fraction))
    stats::median(tail(x, n_tail), na.rm = TRUE)
  }

  estimate_step <- function(t, vin) {
    valid_idx <- is.finite(t) & is.finite(vin)
    t <- t[valid_idx]
    vin <- vin[valid_idx]

    if (length(t) == 0 || length(vin) == 0) {
      return(list(step_start_index = 1L, vin_initial = NA_real_, vin_final = NA_real_))
    }

    vin_initial <- safe_tail_median(head(vin, max(3, ceiling(length(vin) * 0.1))))
    vin_final <- safe_tail_median(vin)
    step_mag <- vin_final - vin_initial

    if (!is.finite(step_mag) || abs(step_mag) < 1e-9) {
      return(list(step_start_index = 1L, vin_initial = vin_initial, vin_final = vin_final))
    }

    threshold <- vin_initial + 0.05 * step_mag
    if (step_mag > 0) {
      start_hits <- which(vin >= threshold)
    } else {
      start_hits <- which(vin <= threshold)
    }

    step_start_index <- if (length(start_hits) > 0) start_hits[1] else 1L
    list(
      step_start_index = step_start_index,
      vin_initial = vin_initial,
      vin_final = vin_final
    )
  }

  crossing_time <- function(t, v, target) {
    if (length(t) < 2 || length(v) < 2) {
      return(NA_real_)
    }

    lower <- v[-length(v)]
    upper <- v[-1]
    hits <- which((lower <= target & upper >= target) | (lower >= target & upper <= target))
    if (length(hits) == 0) {
      return(NA_real_)
    }

    i <- hits[1]
    if (v[i + 1] == v[i]) {
      return(t[i])
    }

    t[i] + (target - v[i]) * (t[i + 1] - t[i]) / (v[i + 1] - v[i])
  }

  safe_relative_error <- function(model_value, sim_value) {
    # Handled potential div by zero and removed the 0.1 arbitrary cutoff
    if (!is.finite(model_value) || !is.finite(sim_value) || sim_value == 0) {
      return(NA_real_)
    }
    abs(model_value - sim_value) / abs(sim_value) * 100
  }

  step_info <- estimate_step(t_slice, v_in_slice)
  step_start_index <- step_info$step_start_index
  v_in_final <- step_info$vin_final
  v_in_initial <- step_info$vin_initial

  target_final <- v_in_final

  # Retained just for reference, but no longer forced into output
  theoretical_zeta <- if (is.finite(resistance) && is.finite(inductance) && is.finite(capacitance) && inductance > 0 && capacitance > 0) {
    (resistance / 2) * sqrt(capacitance / inductance)
  } else {
    NA_real_
  }

  compute_metrics <- function(t, v) {
    if (length(t) == 0 || length(v) == 0) {
      return(list(
        peak_time = NA_real_,
        overshoot = NA_real_,
        rise_time = NA_real_,
        settling_time = NA_real_,
        damping_ratio = NA_real_,
        final_value = NA_real_
      ))
    }

    valid_idx <- is.finite(t) & is.finite(v)
    t <- t[valid_idx]
    v <- v[valid_idx]

    if (step_start_index > length(v)) {
      step_start_index <- 1L
    }

    v_start <- if (step_start_index > 1L) v[step_start_index - 1L] else v[1]
    v_final <- target_final
    v_peak <- max(v)
    peak_index <- which.max(v)
    t_peak <- t[peak_index]
    
    # Calculate overshoot and cap it at 0 minimum
    overshoot <- if (is.finite(v_final) && v_final != 0) {
      ov <- (v_peak - v_final) / abs(v_final) * 100
      max(0, ov) # Prevents negative overshoot
    } else {
      NA_real_
    }

    t_10 <- crossing_time(t, v, 0.1 * v_final)
    t_90 <- crossing_time(t, v, 0.9 * v_final)
    rise_time <- if (is.na(t_10) || is.na(t_90)) NA_real_ else t_90 - t_10

    within_bounds <- abs(v - v_final) <= 0.02 * abs(v_final)
    settling_time <- NA_real_
    for (i in seq_along(t)) {
      if (t[i] > t_peak && all(within_bounds[i:length(v)])) {
        settling_time <- t[i]
        break
      }
    }

    # Estimate empirical damping ratio from overshoot (for underdamped only)
    zeta_est <- NA_real_
    if (is.finite(overshoot) && overshoot > 0) {
      ov_frac <- overshoot / 100
      zeta_est <- sqrt((log(ov_frac)^2) / (pi^2 + log(ov_frac)^2))
    }

    return(list(
      peak_time = t_peak,
      overshoot = overshoot,
      rise_time = rise_time,
      settling_time = settling_time,
      damping_ratio = zeta_est, # Returns calculated zeta instead of theoretical
      final_value = safe_tail_median(v) # Use actual signal final value
    ))
  }

  m_model <- compute_metrics(t_slice, v_model)
  m_sim   <- compute_metrics(t_slice, v_sim)

  # Calculate Steady-State Error (%) compared to the target step input
  calc_ss_error <- function(final_val, target) {
    if (is.finite(final_val) && is.finite(target) && target != 0) {
      abs(final_val - target) / abs(target) * 100
    } else {
      NA_real_
    }
  }
  
  ss_err_model <- calc_ss_error(m_model$final_value, target_final)
  ss_err_sim   <- calc_ss_error(m_sim$final_value, target_final)

  result <- data.frame(
    Metric = c("Peak Time", "Max Overshoot (%)", "Rise Time", "Settling Time", "Damping Ratio", "Steady-State Error (%)"),
    Model  = c(m_model$peak_time, m_model$overshoot, m_model$rise_time, m_model$settling_time, m_model$damping_ratio, ss_err_model),
    Sim    = c(m_sim$peak_time,   m_sim$overshoot,   m_sim$rise_time,   m_sim$settling_time,   m_sim$damping_ratio, ss_err_sim),
    AbsError = c(
      abs(m_model$peak_time - m_sim$peak_time),
      abs(m_model$overshoot - m_sim$overshoot),
      abs(m_model$rise_time - m_sim$rise_time),
      abs(m_model$settling_time - m_sim$settling_time),
      abs(m_model$damping_ratio - m_sim$damping_ratio),
      abs(ss_err_model - ss_err_sim) # Replaced hardcoded NA
    ),
    RelError = c(
      safe_relative_error(m_model$peak_time, m_sim$peak_time),
      safe_relative_error(m_model$overshoot, m_sim$overshoot),
      safe_relative_error(m_model$rise_time, m_sim$rise_time),
      safe_relative_error(m_model$settling_time, m_sim$settling_time),
      safe_relative_error(m_model$damping_ratio, m_sim$damping_ratio),
      safe_relative_error(ss_err_model, ss_err_sim) # Calculate relative difference in SS error
    )
  )

  return(result)
}
