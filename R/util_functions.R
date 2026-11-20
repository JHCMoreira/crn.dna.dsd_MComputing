if (!isNamespace(topenv())) {
    library(ggplot2) # plot()
    library(dplyr) # mutate()
}

React_circuit <- function(circuit, forced_concentrations = NULL, engine = 'desolve', verbose = FALSE, ...) {
  return(react4(
    species   = circuit$species,
    ci        = circuit$ci,
    reactions = circuit$reactions,
    ki        = circuit$ki,
    t         = circuit$t,
    engine = engine,
    verbose = verbose,
    forced_concentrations = forced_concentrations,
    ...
  ))
}

React_4domain <- function(circuit, forced_concentrations = NULL, engine = 'desolve',
                          qmax = 1e6, cmax = 1e7, alpha = 1, beta = 1,
                          auto_buffer = TRUE, dna_kinetics = NULL,
                          verbose = FALSE, ...) {
  return(react_4domain(
    species   = circuit$species,
    ci        = circuit$ci,
    reactions = circuit$reactions,
    ki        = circuit$ki,
    qmax      = qmax,
    cmax      = cmax,
    alpha     = alpha,
    beta      = beta,
    t         = circuit$t,
    auto_buffer = auto_buffer,
    dna_kinetics = dna_kinetics,
    forced_concentrations = forced_concentrations,
    engine = engine,
    verbose = verbose,
    ...
  ))
}

# mu_target_func1 <- function(t) {100}
# mu_target_func2 <- function(t) {
#     pulse_input(t, time = 20, width = 40, amplitude = 50)
# }
# fuzzy_pulse_data <- simulate_ou_process(
  # t = mu_target_func2, 
  # x0 = 50, 
  # theta = 0.5, 
  # mu = 100, 
  # sigma = 15,
  # seed = 42ss
# )
# fuzzy_input_func <- approxfun(x = fuzzy_pulse_data$time, y = fuzzy_pulse_data$value, rule = 2)

Translate_4domain <- function(circuit) {
  return(translate_4domain_crn(
    species   = circuit$species,
    ci        = circuit$ci,
    reactions = circuit$reactions,
    ki        = circuit$ki,
    qmax      = 1e6, # maximum strand displacement rate constant
    cmax      = 1e7, # 1e-4  # starting concentration of auxiliary complexes Gi and Ti
    alpha     = 1,   # DSD timescale versus CRN
    beta      = 1,   # DSD concentration scale versus CRN
    t         = circuit$t
  ))
}


React_stochastic <- function(circuit, volume = 10, seed = NULL,
                             forced_concentrations = NULL,
                             verbose = FALSE, ...) {
  return(react_stochastic_frates(
    species = circuit$species,
    ci = circuit$ci,
    reactions = circuit$reactions,
    ki = circuit$ki,
    t = circuit$t,
    forced_concentrations = forced_concentrations,
    verbose = verbose,
    volume = volume,
    seed = seed,
    ...
  ))
}

React_4domain_stochastic <- function(circuit, volume = 10, seed = NULL,
                                     forced_concentrations = NULL,
                                     qmax = 1e6, cmax = 1e7,
                                     alpha = 1, beta = 1,
                                     auto_buffer = TRUE,
                                     dna_kinetics = NULL,
                                     verbose = FALSE, ...) {
  translated <- translate_4domain_crn(
    species = circuit$species,
    ci = circuit$ci,
    reactions = circuit$reactions,
    ki = circuit$ki,
    qmax = qmax,
    cmax = cmax,
    alpha = alpha,
    beta = beta,
    t = circuit$t,
    auto_buffer = auto_buffer,
    dna_kinetics = dna_kinetics,
    apply_dna_kinetics = TRUE
  )

  behavior <- react_stochastic_frates(
    species = translated$species,
    ci = translated$ci,
    reactions = translated$reactions,
    ki = translated$ki,
    t = circuit$t,
    forced_concentrations = forced_concentrations,
    volume = volume,
    seed = seed,
    verbose = verbose,
    ...
  )

  list(
    behavior = behavior,
    species = translated$species,
    ci = translated$ci,
    reactions = translated$reactions,
    ki = translated$ki
  )
}


Plot_behavior <- function(
    result, circuit = NULL, gate_numbers = NULL, 
    min = NULL, max = NULL, 
    species = NULL, species_dotted = NULL, 
    title = "Behavior", intercept = FALSE,
    normalize = FALSE
) {
  
  if (normalize) {
    result_scaled <- result
    # Apply min-max scaling to all columns except the first one (Time)
    result_scaled[, -1] <- lapply(result[, -1], function(x) {
      x_min <- min(x, na.rm = TRUE)
      x_max <- max(x, na.rm = TRUE)
      
      if (x_max == x_min) {
        return(x) # Prevent division by zero for constant columns
      }
      return((x - x_min) / (x_max - x_min))
    })
    
    result <- result_scaled
    y_label_text <- 'Normalized Concentration (M)'
  } else {
    y_label_text <- 'Concentration (M)'
  }

  species_to_plot <- c()
  
  if (!is.null(species) && length(species) > 0) {
    species_to_plot <- species 
  } else if (!is.null(gate_numbers) && length(gate_numbers) > 0 && !is.null(circuit)) {
    for (i in gate_numbers) {
      species_to_plot <- c(species_to_plot, circuit$gates[[i]]$species)
    }
    species_to_plot <- unique(species_to_plot)
  } else if (!is.null(circuit)) {
    species_to_plot <- circuit$species
  }
  
  dotted_species_to_plot <- NULL
  if (!is.null(species_dotted) && length(species_dotted) > 0) {
    dotted_species_to_plot <- species_dotted
  }

  g <- plot_behavior(result, 
                     species = species_to_plot,
                     species_dotted = dotted_species_to_plot,
                     title = title, 
                     x_label     = 'Time (s)',
                     y_label     = y_label_text,
                     legend_name = 'Species',
                     geom_list   = c('line', 'point'),
                     variable_line_type = FALSE,
                     variable_point_type = TRUE
  )

  if (intercept || !is.null(min)) {
    min_val <- ifelse(is.null(min), 0, min) 
    g <- g + geom_hline(yintercept = min_val, linetype = "dashed", color = "green", linewidth = 1)
  }
  
  if (intercept || !is.null(max)) {
    # If normalized, the absolute max should probably be 1. Otherwise, default to 10.
    max_val <- ifelse(is.null(max), ifelse(normalize, 1, 10), max) 
    g <- g + geom_hline(yintercept = max_val, linetype = "dashed", color = "red", linewidth = 1)
  }

  # print(g)
  return(g)
}

# reactions <- append_reaction(reactions, "a -_. b")
append_reaction <- function (reactions, new_reaction) {
  reactions[[length(reactions)+1]] <- new_reaction
  return (reactions)
}
