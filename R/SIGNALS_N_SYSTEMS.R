#' Digital-to-Analog Converter for Molecular Computing
#' 
#' Converts binary digital signals to analog concentrations
#' Uses weighted summation approach similar to R-2R ladder
#' 
#' @param name Module name
#' @param bit_inputs Vector of binary input species names (LSB to MSB)
#' @param output_analog Name of analog output species
#' @param bit_concentrations Vector of concentrations for each bit (high=1, low=0)
#' @param reference_voltage Maximum analog output concentration
#' @param rate Reaction rate constant
#' @return DAC module structure
Make_DAC <- function(name, bit_inputs, output_analog, 
                     bit_concentrations, reference_voltage, rate) {
  
  n_bits <- length(bit_inputs)
  
  # Create weighted intermediate species for each bit
  intermediates <- lapply(1:n_bits, function(i) {
    jn(name, '_weight', i)
  })
  
  species <- list(
    inputs = bit_inputs,
    intermediates = intermediates,
    output = output_analog,
    fuel = jn(name, '_fuel')
  )
  
  # Initial concentrations
  # Fuel species represents the weighting factors (2^0, 2^1, 2^2, ...)
  weights <- sapply(0:(n_bits-1), function(i) reference_voltage / (2^(n_bits)) * (2^i))
  
  ci <- c(bit_concentrations,  # inputs
          weights,              # intermediates (weights)
          0,                    # output starts at 0
          reference_voltage)    # fuel for conversion
  
  # Reactions: Each bit contributes weighted amount to output
  # Bit[i] + Weight[i] -> Output + Weight[i]  (catalytic production)
  reactions <- c()
  for (i in 1:n_bits) {
    reactions <- c(reactions,
      jn(bit_inputs[i], ' + ', intermediates[[i]], ' -> ',
         output_analog, ' + ', intermediates[[i]])
    )
  }
  
  # Decay reaction to allow output to settle
  reactions <- c(reactions,
    jn(output_analog, ' -> 0')
  )
  
  ki <- c(rep(rate, n_bits), rate * 0.1)  # slower decay
  
  dac_module <- list(
    name = name,
    species = species,
    reactions = reactions,
    ci = ci,
    ki = ki,
    n_bits = n_bits
  )
  
  return(dac_module)
}

#' Analog-to-Digital Converter for Molecular Computing
#' 
#' Converts analog concentration to binary digital signals
#' Uses flash ADC approach with threshold comparators
#' 
#' @param name Module name
#' @param input_analog Name of analog input species
#' @param bit_outputs Vector of binary output species names (LSB to MSB)
#' @param reference_voltage Maximum input concentration (full scale)
#' @param rate Reaction rate constant
#' @return ADC module structure
Make_ADC <- function(name, input_analog, bit_outputs, 
                     reference_voltage, rate) {
  
  n_bits <- length(bit_outputs)
  n_comparators <- 2^n_bits - 1
  
  # Create threshold references
  thresholds <- lapply(1:n_comparators, function(i) {
    jn(name, '_threshold', i)
  })
  
  # Create comparator outputs
  comparator_outs <- lapply(1:n_comparators, function(i) {
    jn(name, '_comp', i)
  })
  
  species <- list(
    input = input_analog,
    thresholds = thresholds,
    comparator_outs = comparator_outs,
    bit_outputs = bit_outputs,
    encoder_intermediates = lapply(1:n_bits, function(i) jn(name, '_enc', i))
  )
  
  # Threshold levels (evenly spaced)
  threshold_levels <- sapply(1:n_comparators, function(i) {
    reference_voltage * i / (2^n_bits)
  })
  
  ci <- c(
    0,                                    # input (set externally)
    threshold_levels,                      # thresholds
    rep(0, n_comparators),                # comparator outputs
    rep(0, n_bits),                       # bit outputs
    rep(reference_voltage, n_bits)        # encoder intermediates
  )
  
  reactions <- c()
  
  # Comparator reactions: if input > threshold, produce comparator output
  # Input + Threshold[i] -> ComparatorOut[i] + Threshold[i]
  for (i in 1:n_comparators) {
    reactions <- c(reactions,
      jn(input_analog, ' + ', thresholds[[i]], ' -> ',
         comparator_outs[[i]], ' + ', thresholds[[i]])
    )
  }
  
  # Priority encoder logic (thermometer to binary)
  # This is simplified - you'd need proper digital logic gates
  # For n=3 bits, you need logic to convert 7 comparator outputs to 3 bits
  
  # Example for 3-bit ADC (8 levels):
  # Bit[2] (MSB) = Comp[4] | Comp[5] | Comp[6] | Comp[7]
  # Bit[1]       = Comp[2] | Comp[3] | !Comp[4] & (Comp[6] | Comp[7])
  # Bit[0] (LSB) = Comp[1] | !Comp[2] & Comp[3] | ...
  
  # Simplified: Use thresholds directly mapped to bits
  for (bit in 1:n_bits) {
    # If input exceeds the threshold for this bit position, activate bit
    threshold_for_bit <- reference_voltage * (2^(bit-1)) / (2^n_bits)
    reactions <- c(reactions,
      jn(input_analog, ' + ', species$encoder_intermediates[[bit]], ' -> ',
         bit_outputs[bit], ' + ', species$encoder_intermediates[[bit]])
    )
    # Decay
    reactions <- c(reactions,
      jn(bit_outputs[bit], ' -> 0')
    )
  }
  
  ki <- c(rep(rate, length(reactions)))
  
  adc_module <- list(
    name = name,
    species = species,
    reactions = reactions,
    ci = ci,
    ki = ki,
    n_bits = n_bits
  )
  
  return(adc_module)
}

#' Register/Flip-Flop for Storing Molecular Concentrations
#' 
#' Stores value on clock edge, implements delay for pipeline stages
#' Essential for non-blocking assignments in Verilog-like behavior
#' 
#' @param name Register name
#' @param input_data Input species name
#' @param clock_enable Clock signal species name
#' @param output_data Output species name (stored value)
#' @param initial_value Initial stored concentration
#' @param rate Reaction rate
#' @return Register module structure
Make_Register <- function(name, input_data, clock_enable, output_data,
                          initial_value, rate) {
  
  species <- list(
    input = input_data,
    clock = clock_enable,
    storage = jn(name, '_storage'),  # internal storage
    output = output_data,
    transfer_gate = jn(name, '_gate')  # controls transfer
  )
  
  ci <- c(
    0,              # input (set externally)
    0,              # clock (set externally)
    initial_value,  # storage holds initial value
    initial_value,  # output reflects storage
    1               # transfer gate enabled
  )
  
  reactions <- c()
  
  # On clock edge: transfer input to storage
  # Clock + Input + Gate -> Storage + Gate
  reactions <- c(reactions,
    jn(clock_enable, ' + ', input_data, ' + ', species$transfer_gate, ' -> ',
       species$storage, ' + ', species$transfer_gate)
  )
  
  # Storage continuously drives output
  # Storage -> Output
  reactions <- c(reactions,
    jn(species$storage, ' -> ', output_data)
  )
  
  # Output decay (to allow updates)
  reactions <- c(reactions,
    jn(output_data, ' -> 0')
  )
  
  # Storage decay (slow, to maintain value between clocks)
  reactions <- c(reactions,
    jn(species$storage, ' -> 0')
  )
  
  ki <- c(
    rate * 10,      # fast transfer on clock
    rate * 5,       # fast output drive
    rate * 0.5,     # slow output decay
    rate * 0.01     # very slow storage decay (memory retention)
  )
  
  register_module <- list(
    name = name,
    species = species,
    reactions = reactions,
    ci = ci,
    ki = ki
  )
  
  return(register_module)
}

#' Phased Reaction Simulation
#' 
#' Extends react() to support phase-based execution
#' Allows different reaction sets to be active during different time windows
#' 
#' @param species Species vector
#' @param ci Initial concentrations
#' @param reaction_phases List of phase definitions, each containing:
#'   - phase_name: identifier
#'   - time_start: when phase begins (or species trigger)
#'   - time_end: when phase ends
#'   - reactions: reactions active during this phase
#'   - ki: rate constants for these reactions
#' @param t Time vector
#' @return Behavior matrix
react_phased <- function(species, ci, reaction_phases, t, verbose = FALSE) {
  
  # This is a conceptual implementation
  # You would need to integrate with deSolve's event system
  
  # Create an event function that switches reaction sets
  eventfunc <- function(t, y, parms) {
    # Determine which phase is active
    active_phase <- NULL
    for (phase in reaction_phases) {
      if (t >= phase$time_start && t < phase$time_end) {
        active_phase <- phase
        break
      }
    }
    
    if (!is.null(active_phase)) {
      # Switch to this phase's reactions
      # This would require dynamic reaction set switching
      # (not trivially supported by deSolve, needs custom implementation)
    }
    
    return(y)
  }
  
  # For now, you can approximate phases by having conditional reactions
  # that are modulated by phase-indicator species (like clock signals)
  
  # Example: Reaction is only active when clock_phase species > threshold
  # Rate = base_rate * clock_phase
  # This makes reactions phase-dependent
  
  # Call standard react() with phase-modulated reactions
  all_reactions <- c()
  all_ki <- c()
  
  for (phase in reaction_phases) {
    # Add phase-gated reactions
    for (i in seq_along(phase$reactions)) {
      # Modulate reaction with phase species
      gated_reaction <- paste0(phase$reactions[i], 
                               ' [WHEN ', phase$phase_name, ' > 0]')
      all_reactions <- c(all_reactions, gated_reaction)
      all_ki <- c(all_ki, phase$ki[i])
    }
  }
  
  # Note: The above is pseudocode. In practice, you'd implement gating
  # by adding the phase species as a catalyst in each reaction
  
  return(react(species, ci, all_reactions, all_ki, t, verbose))
}

#' Clock Generator with Multiple Phases
#' 
#' Creates clock signals with positive edge, settling time, and negative edge
#' Similar to Dalchau oscillator but with explicit phase outputs
#' 
#' @param name Clock name
#' @param period Oscillation period
#' @param phases Vector of phase names: c('pos_edge', 'hold', 'neg_edge')
#' @param phase_concentrations Peak concentrations for each phase
#' @param rate Base reaction rate
#' @return Clock module structure
Make_Clock_System <- function(name, period, phases, phase_concentrations, rate) {
  
  # Use three-species oscillator (similar to Dalchau)
  # But add logic to create distinct phase signals
  
  species <- list(
    oscillator = list(
      A = jn(name, '_osc_A'),
      B = jn(name, '_osc_B'),
      C = jn(name, '_osc_C')
    ),
    phases = phases,  # pos_edge, hold, neg_edge
    comparators = list(
      pos_trigger = jn(name, '_pos_trigger'),
      neg_trigger = jn(name, '_neg_trigger')
    )
  )
  
  ci <- c(
    phase_concentrations[1],  # A
    phase_concentrations[2],  # B
    phase_concentrations[3],  # C
    rep(0, length(phases)),   # phase outputs start at 0
    phase_concentrations[1]/2,  # pos_trigger threshold
    phase_concentrations[1]/2   # neg_trigger threshold
  )
  
  # Core oscillator (Dalchau-style)
  reactions <- c(
    # B + A -> B + B
    jn(species$oscillator$B, ' + ', species$oscillator$A, ' -> ',
       species$oscillator$B, ' + ', species$oscillator$B),
    # C + B -> C + C
    jn(species$oscillator$C, ' + ', species$oscillator$B, ' -> ',
       species$oscillator$C, ' + ', species$oscillator$C),
    # A + C -> A + A
    jn(species$oscillator$A, ' + ', species$oscillator$C, ' -> ',
       species$oscillator$A, ' + ', species$oscillator$A),
    # Degradation
    jn(species$oscillator$A, ' -> 0'),
    jn(species$oscillator$B, ' -> 0'),
    jn(species$oscillator$C, ' -> 0')
  )
  
  # Phase detection: rising edge of A triggers pos_edge
  reactions <- c(reactions,
    jn(species$oscillator$A, ' + ', species$comparators$pos_trigger, ' -> ',
       phases[1], ' + ', species$comparators$pos_trigger),
    # Decay pos_edge
    jn(phases[1], ' -> 0')
  )
  
  # Hold phase - active when pos_edge is decaying
  reactions <- c(reactions,
    jn(phases[1], ' -> ', phases[2]),
    jn(phases[2], ' -> 0')
  )
  
  # Neg edge - when A falls
  reactions <- c(reactions,
    jn(species$comparators$neg_trigger, ' + ', species$oscillator$B, ' -> ',
       phases[3], ' + ', species$comparators$neg_trigger),
    jn(phases[3], ' -> 0')
  )
  
  # Adjust rate based on desired period
  base_rate <- rate * (2 * 3.14159) / period
  ki <- c(rep(base_rate, 6),  # oscillator reactions
          rep(rate, length(reactions) - 6))
  
  clock_module <- list(
    name = name,
    species = species,
    reactions = reactions,
    ci = ci,
    ki = ki
  )
  
  return(clock_module)
}

# =============================================================================
# Digital Signal Processing Filter Components
# Based on: Jiang, Riedel, Parhi (2013)
# "Digital Signal Processing with Molecular Reactions"
# =============================================================================

#' Phase Generator for Filter Coordination
#' 
#' Generates two alternating phases to coordinate sequential filter operations
#' Uses the proven Dalchau oscillator design without modification
#' Phase 1: Sample/compute phase (derived from oscillator species A)
#' Phase 2: Hold/transfer phase (derived from oscillator species B)
#' 
#' @param name Module name
#' @param period Oscillation period for phase switching (not used directly, kept for API compatibility)
#' @param crange Concentration range for phase signals
#' @param rate Base reaction rate for oscillator
#' @return Phase generator module structure with phase1 and phase2 species
#' @examples Make_Phase_Generator('phase_gen', period = 100, crange = 10, rate = 1e-3)
Make_Phase_Generator <- function(name, period, crange, rate) {
  
  # Helper function for string concatenation
  jn <- function(...) { paste(..., sep = '') }
  
  # Use proven Dalchau oscillator design (no degradation, just autocatalysis)
  # This is the working implementation from ANALOG_GATE_LIB.R
  
  # Create oscillator species
  species <- list(
    # Oscillator core species (A, B, C form rock-paper-scissors cycle)
    osc_A = jn(name, '_osc_A'),
    osc_B = jn(name, '_osc_B'),
    osc_C = jn(name, '_osc_C'),
    # Phase outputs (directly use oscillator species as phases)
    phase1 = jn(name, '_osc_A'),  # Phase 1 is oscillator A
    phase2 = jn(name, '_osc_B')   # Phase 2 is oscillator B
  )
  
  # Initial concentrations - staggered to ensure oscillation
  ci <- c(
    crange * 0.9,  # osc_A starts high
    crange * 0.5,  # osc_B starts medium
    crange * 0.3   # osc_C starts lower
  )
  
  # Three-species oscillator reactions (Dalchau-style, proven to work)
  # Rock-paper-scissors: A beats C, B beats A, C beats B
  reactions <- c(
    # B + A -> B + B (B suppresses A)
    jn(species$osc_B, ' + ', species$osc_A, ' -> ', species$osc_B, ' + ', species$osc_B),
    # C + B -> C + C (C suppresses B)
    jn(species$osc_C, ' + ', species$osc_B, ' -> ', species$osc_C, ' + ', species$osc_C),
    # A + C -> A + A (A suppresses C)
    jn(species$osc_A, ' + ', species$osc_C, ' -> ', species$osc_A, ' + ', species$osc_A)
  )
  
  # All three reactions have same rate (from working Dalchau implementation)
  ki <- c(rate, rate, rate)
  
  phase_gen_module <- list(
    name = name,
    species = species,
    reactions = reactions,
    ci = ci,
    ki = ki,
    period = period
  )
  
  return(phase_gen_module)
}

#' Delay Unit (Z^-1 Operation)
#' 
#' Implements one sample delay for digital filters
#' Stores input value and outputs it in the next phase
#' Essential for FIR and IIR filter implementations
#' 
#' @param name Module name
#' @param nameInput Input species name
#' @param nameOutput Output species name (delayed signal)
#' @param phase1 Phase 1 species name (sampling phase)
#' @param phase2 Phase 2 species name (output phase)
#' @param cinput Initial input concentration
#' @param crange Concentration range
#' @param rate Reaction rate
#' @return Delay unit module structure
#' @examples Make_Delay_Unit('delay1', 'X_in', 'X_delayed', 'phase1', 'phase2', 0, 10, 1e-3)
Make_Delay_Unit <- function(name, nameInput, nameOutput, phase1, phase2, 
                            cinput, crange, rate) {
  
  # Helper function for string concatenation
  jn <- function(...) { paste(..., sep = '') }
  
  species <- list(
    input = nameInput,
    storage = jn(name, '_storage'),  # Internal storage species
    output = nameOutput,
    phase1 = phase1,
    phase2 = phase2
  )
  
  ci <- c(
    cinput,    # input
    0,         # storage (initially empty)
    0,         # output
    0,         # phase1 (set externally)
    0          # phase2 (set externally)
  )
  
  reactions <- c()
  
  # Phase 1: Sample input into storage
  # Input + Phase1 -> Storage + Phase1 (catalytic)
  reactions <- c(reactions,
    jn(nameInput, ' + ', phase1, ' -> ', species$storage, ' + ', phase1)
  )
  
  # Phase 2: Transfer storage to output
  # Storage + Phase2 -> Output + Phase2 (catalytic transfer)
  reactions <- c(reactions,
    jn(species$storage, ' + ', phase2, ' -> ', nameOutput, ' + ', phase2)
  )
  
  # Storage decay (slow to maintain value between phases)
  reactions <- c(reactions,
    jn(species$storage, ' -> 0')
  )
  
  # Output decay
  reactions <- c(reactions,
    jn(nameOutput, ' -> 0')
  )
  
  ki <- c(
    rate * 5,      # fast sampling
    rate * 5,      # fast output
    rate * 0.05,   # slow storage decay
    rate * 0.5     # moderate output decay
  )
  
  delay_module <- list(
    name = name,
    species = species,
    reactions = reactions,
    ci = ci,
    ki = ki
  )
  
  return(delay_module)
}

#' Coefficient Multiplier
#' 
#' Multiplies input signal by a fixed coefficient (tap weight)
#' Used in FIR and IIR filters for scaling taps
#' 
#' @param name Module name
#' @param nameInput Input species name
#' @param nameOutput Output species name
#' @param coefficient Multiplication coefficient (0 to 1)
#' @param cinput Initial input concentration
#' @param crange Concentration range
#' @param rate Reaction rate
#' @return Coefficient multiplier module structure
#' @examples Make_Coefficient_Multiplier('mult1', 'X', 'X_scaled', 0.5, 0, 10, 1e-3)
Make_Coefficient_Multiplier <- function(name, nameInput, nameOutput, coefficient,
                                        cinput, crange, rate) {
  
  # Helper function for string concatenation
  jn <- function(...) { paste(..., sep = '') }
  
  species <- list(
    input = nameInput,
    catalyst = jn(name, '_catalyst'),  # Represents the coefficient
    output = nameOutput
  )
  
  # Catalyst concentration represents the coefficient
  catalyst_conc <- coefficient * crange
  
  ci <- c(
    cinput,         # input
    catalyst_conc,  # catalyst (fixed coefficient)
    0               # output
  )
  
  reactions <- c()
  
  # Multiplication: Input + Catalyst -> Output + Catalyst
  # The rate is adjusted by catalyst concentration to achieve scaling
  reactions <- c(reactions,
    jn(nameInput, ' + ', species$catalyst, ' -> ', nameOutput, ' + ', species$catalyst)
  )
  
  # Output decay
  reactions <- c(reactions,
    jn(nameOutput, ' -> 0')
  )
  
  ki <- c(
    rate * 2,    # multiplication
    rate * 0.5   # output decay
  )
  
  mult_module <- list(
    name = name,
    species = species,
    reactions = reactions,
    ci = ci,
    ki = ki,
    coefficient = coefficient
  )
  
  return(mult_module)
}

#' Two-Tap FIR Filter
#' 
#' Implements a two-tap Finite Impulse Response filter
#' Transfer function: y[n] = a0*x[n] + a1*x[n-1]
#' 
#' This is a first-order FIR filter with two coefficients.
#' Can implement simple filtering operations like moving average,
#' differencing, or basic low-pass/high-pass filtering.
#' 
#' @param name Filter name
#' @param nameInput Input signal species name
#' @param nameOutput Output signal species name
#' @param coeff0 Coefficient for current sample (a0)
#' @param coeff1 Coefficient for delayed sample (a1)
#' @param crange Concentration range
#' @param rate Base reaction rate
#' @return Two-tap filter module with all components and nested species structure
#' @examples 
#' # Moving average filter
#' Make_TwoTap_Filter('fir1', 'X_in', 'Y_out', 0.5, 0.5, 10, 1e-3)
#' # Differencing filter  
#' Make_TwoTap_Filter('fir2', 'X_in', 'Y_out', 1.0, -1.0, 10, 1e-3)
Make_TwoTap_Filter <- function(name, nameInput, nameOutput, 
                               coeff0 = 0.5, coeff1 = 0.5,
                               crange = 10, rate = 1e-3) {
  
  # Helper function for string concatenation
  jn <- function(...) { paste(..., sep = '') }
  
  # Source the required library for analog gates (adders)
  # Note: This assumes ANALOG_GATE_LIB.R defines Make_Adder2In_Song
  
  # Create phase generator
  phase_gen <- Make_Phase_Generator(
    name = jn(name, '_phase_gen'),
    period = 100,
    crange = crange,
    rate = rate
  )
  
  # Create delay unit (z^-1)
  delay <- Make_Delay_Unit(
    name = jn(name, '_delay'),
    nameInput = nameInput,
    nameOutput = jn(name, '_x_delayed'),
    phase1 = phase_gen$species$phase1,
    phase2 = phase_gen$species$phase2,
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  # Create coefficient multipliers for both taps
  mult0 <- Make_Coefficient_Multiplier(
    name = jn(name, '_mult0'),
    nameInput = nameInput,
    nameOutput = jn(name, '_tap0'),
    coefficient = coeff0,
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  mult1 <- Make_Coefficient_Multiplier(
    name = jn(name, '_mult1'),
    nameInput = delay$species$output,
    nameOutput = jn(name, '_tap1'),
    coefficient = coeff1,
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  # Create adder to sum the two taps
  # Note: Using simplified addition reactions instead of Make_Adder2In_Song
  # to avoid complex dependencies
  adder_species <- list(
    tap0 = mult0$species$output,
    tap1 = mult1$species$output,
    output = nameOutput
  )
  
  adder_reactions <- c(
    # Tap0 -> Output (transfer)
    jn(mult0$species$output, ' -> ', nameOutput),
    # Tap1 -> Output (transfer)
    jn(mult1$species$output, ' -> ', nameOutput),
    # Output decay
    jn(nameOutput, ' -> 0')
  )
  
  adder_ki <- c(rate, rate, rate * 0.3)
  adder_ci <- c(0, 0, 0)
  
  # Combine all components
  all_species <- c(
    unlist(phase_gen$species),
    unlist(delay$species),
    unlist(mult0$species),
    unlist(mult1$species),
    adder_species$output
  )
  
  # Remove duplicates (shared species like input, phase1, phase2)
  all_species_unique <- unique(all_species)
  
  all_reactions <- c(
    phase_gen$reactions,
    delay$reactions,
    mult0$reactions,
    mult1$reactions,
    adder_reactions
  )
  
  all_ki <- c(
    phase_gen$ki,
    delay$ki,
    mult0$ki,
    mult1$ki,
    adder_ki
  )
  
  all_ci <- c(
    phase_gen$ci,
    delay$ci,
    mult0$ci,
    mult1$ci,
    adder_ci
  )
  
  # Create a mapping to handle duplicate species
  species_map <- list()
  ci_final <- c()
  
  for (sp in all_species_unique) {
    # Find all occurrences of this species
    indices <- which(all_species == sp)
    # Use the first non-zero concentration if available
    concs <- all_ci[indices]
    final_conc <- ifelse(any(concs != 0), max(concs), 0)
    ci_final <- c(ci_final, final_conc)
  }
  
  filter_module <- list(
    name = name,
    species = all_species_unique,
    species_nested = list(
      input = nameInput,
      output = nameOutput,
      phase_gen = phase_gen$species,
      delay = delay$species,
      mult0 = mult0$species,
      mult1 = mult1$species
    ),
    reactions = all_reactions,
    ci = ci_final,
    ki = all_ki,
    components = list(
      phase_gen = phase_gen,
      delay = delay,
      mult0 = mult0,
      mult1 = mult1
    ),
    coeff0 = coeff0,
    coeff1 = coeff1
  )
  
  return(filter_module)
}

#' Biquad (Second-Order IIR) Filter
#' 
#' Implements a second-order Infinite Impulse Response (IIR) filter
#' Transfer function: H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
#' 
#' This is the standard biquad filter structure used in digital signal processing.
#' Can implement various filter types by selecting appropriate coefficients:
#' - Low-pass filter: attenuates high frequencies
#' - High-pass filter: attenuates low frequencies  
#' - Band-pass filter: passes specific frequency range
#' - Notch filter: rejects specific frequency
#' 
#' The implementation uses the Direct Form II structure which is efficient
#' and minimizes the number of delay elements needed.
#' 
#' Reference: Jiang, Riedel, Parhi (2013) "Digital Signal Processing with Molecular Reactions"
#' Section on IIR filters and biquad structures
#' 
#' @param name Filter name
#' @param nameInput Input signal species name
#' @param nameOutput Output signal species name
#' @param b0 Numerator coefficient 0 (feedforward)
#' @param b1 Numerator coefficient 1 (feedforward, z^-1 term)
#' @param b2 Numerator coefficient 2 (feedforward, z^-2 term)
#' @param a1 Denominator coefficient 1 (feedback, z^-1 term)
#' @param a2 Denominator coefficient 2 (feedback, z^-2 term)
#' @param crange Concentration range
#' @param rate Base reaction rate
#' @return Biquad filter module with all components
#' @examples
#' # Low-pass Butterworth filter (fc = 0.1*fs, where fs is sampling frequency)
#' Make_Biquad_Filter('lpf', 'X', 'Y', 
#'                    b0=0.0201, b1=0.0402, b2=0.0201, 
#'                    a1=-1.5610, a2=0.6414, 
#'                    crange=10, rate=1e-3)
#' 
#' # High-pass Butterworth filter
#' Make_Biquad_Filter('hpf', 'X', 'Y',
#'                    b0=0.8008, b1=-1.6016, b2=0.8008,
#'                    a1=-1.5610, a2=0.6414,
#'                    crange=10, rate=1e-3)
#'
#' # Band-pass filter  
#' Make_Biquad_Filter('bpf', 'X', 'Y',
#'                    b0=0.0318, b1=0, b2=-0.0318,
#'                    a1=-1.5610, a2=0.6364,
#'                    crange=10, rate=1e-3)
Make_Biquad_Filter <- function(name, nameInput, nameOutput,
                               b0 = 0.0201, b1 = 0.0402, b2 = 0.0201,
                               a1 = -1.5610, a2 = 0.6414,
                               crange = 10, rate = 1e-3) {
  
  # Helper function for string concatenation
  jn <- function(...) { paste(..., sep = '') }
  
  # Biquad Direct Form II structure:
  # w[n] = x[n] - a1*w[n-1] - a2*w[n-2]   (internal state)
  # y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2] (output)
  
  # Create phase generator for coordination
  phase_gen <- Make_Phase_Generator(
    name = jn(name, '_phase_gen'),
    period = 100,
    crange = crange,
    rate = rate
  )
  
  # Internal state species (w[n])
  w_current = jn(name, '_w_current')
  w_z1 = jn(name, '_w_z1')  # w[n-1]
  w_z2 = jn(name, '_w_z2')  # w[n-2]
  
  # Create delay units for state variables
  delay_w1 <- Make_Delay_Unit(
    name = jn(name, '_delay_w1'),
    nameInput = w_current,
    nameOutput = w_z1,
    phase1 = phase_gen$species$phase1,
    phase2 = phase_gen$species$phase2,
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  delay_w2 <- Make_Delay_Unit(
    name = jn(name, '_delay_w2'),
    nameInput = w_z1,
    nameOutput = w_z2,
    phase1 = phase_gen$species$phase1,
    phase2 = phase_gen$species$phase2,
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  # Feedback path: -a1*w[n-1] and -a2*w[n-2]
  # Note: Using abs() for coefficients and handling sign in reaction logic
  # This allows the same multiplier module to work with signed coefficients
  # Alternative: implement a signed coefficient multiplier (future enhancement)
  mult_a1 <- Make_Coefficient_Multiplier(
    name = jn(name, '_mult_a1'),
    nameInput = w_z1,
    nameOutput = jn(name, '_a1_w_z1'),
    coefficient = abs(a1),
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  mult_a2 <- Make_Coefficient_Multiplier(
    name = jn(name, '_mult_a2'),
    nameInput = w_z2,
    nameOutput = jn(name, '_a2_w_z2'),
    coefficient = abs(a2),
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  # Feedforward path: b0*w[n], b1*w[n-1], b2*w[n-2]
  mult_b0 <- Make_Coefficient_Multiplier(
    name = jn(name, '_mult_b0'),
    nameInput = w_current,
    nameOutput = jn(name, '_b0_w'),
    coefficient = abs(b0),
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  mult_b1 <- Make_Coefficient_Multiplier(
    name = jn(name, '_mult_b1'),
    nameInput = w_z1,
    nameOutput = jn(name, '_b1_w_z1'),
    coefficient = abs(b1),
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  mult_b2 <- Make_Coefficient_Multiplier(
    name = jn(name, '_mult_b2'),
    nameInput = w_z2,
    nameOutput = jn(name, '_b2_w_z2'),
    coefficient = abs(b2),
    cinput = 0,
    crange = crange,
    rate = rate
  )
  
  # State update and output computation reactions
  feedback_species <- c(nameInput, mult_a1$species$output, mult_a2$species$output, w_current)
  feedforward_species <- c(mult_b0$species$output, mult_b1$species$output, mult_b2$species$output, nameOutput)
  
  # Reactions for state computation: w[n] = x[n] - a1*w[n-1] - a2*w[n-2]
  state_reactions <- c()
  
  # Add input to state (x[n] -> w[n])
  state_reactions <- c(state_reactions,
    jn(nameInput, ' -> ', w_current)
  )
  
  # Subtract feedback terms based on the equation: w[n] = x[n] - a1*w[n-1] - a2*w[n-2]
  # When a1 is negative, -a1*w[n-1] becomes positive (adding to w[n])
  # When a1 is positive, -a1*w[n-1] is negative (subtracting from w[n])
  if (a1 < 0) {
    # a1 is negative, so -a1*w[n-1] is positive: add |a1|*w[n-1] to w[n]
    state_reactions <- c(state_reactions,
      jn(mult_a1$species$output, ' -> ', w_current)
    )
  } else {
    # a1 is positive, so -a1*w[n-1] is negative: subtract a1*w[n-1] from w[n]
    state_reactions <- c(state_reactions,
      jn(mult_a1$species$output, ' + ', w_current, ' -> 0')
    )
  }
  
  if (a2 < 0) {
    # a2 is negative, so -a2*w[n-2] is positive: add |a2|*w[n-2] to w[n]
    state_reactions <- c(state_reactions,
      jn(mult_a2$species$output, ' -> ', w_current)
    )
  } else {
    # a2 is positive, so -a2*w[n-2] is negative: subtract a2*w[n-2] from w[n]
    state_reactions <- c(state_reactions,
      jn(mult_a2$species$output, ' + ', w_current, ' -> 0')
    )
  }
  
  # w[n] decay (to allow updates)
  state_reactions <- c(state_reactions,
    jn(w_current, ' -> 0')
  )
  
  state_ki <- rep(rate, length(state_reactions))
  
  # Reactions for output: y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2]
  output_reactions <- c()
  
  if (b0 >= 0) {
    output_reactions <- c(output_reactions,
      jn(mult_b0$species$output, ' -> ', nameOutput)
    )
  } else {
    # Negative coefficient (unusual for b0, but handle it)
    output_reactions <- c(output_reactions,
      jn(mult_b0$species$output, ' + ', nameOutput, ' -> 0')
    )
  }
  
  if (b1 >= 0) {
    output_reactions <- c(output_reactions,
      jn(mult_b1$species$output, ' -> ', nameOutput)
    )
  } else {
    output_reactions <- c(output_reactions,
      jn(mult_b1$species$output, ' + ', nameOutput, ' -> 0')
    )
  }
  
  if (b2 >= 0) {
    output_reactions <- c(output_reactions,
      jn(mult_b2$species$output, ' -> ', nameOutput)
    )
  } else {
    output_reactions <- c(output_reactions,
      jn(mult_b2$species$output, ' + ', nameOutput, ' -> 0')
    )
  }
  
  # Output decay
  output_reactions <- c(output_reactions,
    jn(nameOutput, ' -> 0')
  )
  
  output_ki <- rep(rate, length(output_reactions))
  
  # Combine all components
  all_species <- c(
    unlist(phase_gen$species),
    w_current, w_z1, w_z2,
    unlist(delay_w1$species),
    unlist(delay_w2$species),
    unlist(mult_a1$species),
    unlist(mult_a2$species),
    unlist(mult_b0$species),
    unlist(mult_b1$species),
    unlist(mult_b2$species),
    nameInput,
    nameOutput
  )
  
  all_species_unique <- unique(all_species)
  
  all_reactions <- c(
    phase_gen$reactions,
    delay_w1$reactions,
    delay_w2$reactions,
    mult_a1$reactions,
    mult_a2$reactions,
    mult_b0$reactions,
    mult_b1$reactions,
    mult_b2$reactions,
    state_reactions,
    output_reactions
  )
  
  all_ki <- c(
    phase_gen$ki,
    delay_w1$ki,
    delay_w2$ki,
    mult_a1$ki,
    mult_a2$ki,
    mult_b0$ki,
    mult_b1$ki,
    mult_b2$ki,
    state_ki,
    output_ki
  )
  
  all_ci <- c(
    phase_gen$ci,
    0, 0, 0,  # w_current, w_z1, w_z2
    delay_w1$ci,
    delay_w2$ci,
    mult_a1$ci,
    mult_a2$ci,
    mult_b0$ci,
    mult_b1$ci,
    mult_b2$ci,
    0,  # input
    0   # output
  )
  
  # Create a mapping to handle duplicate species
  ci_final <- c()
  for (sp in all_species_unique) {
    indices <- which(all_species == sp)
    concs <- all_ci[indices]
    final_conc <- ifelse(any(concs != 0), max(concs), 0)
    ci_final <- c(ci_final, final_conc)
  }
  
  biquad_module <- list(
    name = name,
    species = all_species_unique,
    species_nested = list(
      input = nameInput,
      output = nameOutput,
      w_current = w_current,
      w_z1 = w_z1,
      w_z2 = w_z2,
      phase_gen = phase_gen$species
    ),
    reactions = all_reactions,
    ci = ci_final,
    ki = all_ki,
    components = list(
      phase_gen = phase_gen,
      delay_w1 = delay_w1,
      delay_w2 = delay_w2,
      mult_a1 = mult_a1,
      mult_a2 = mult_a2,
      mult_b0 = mult_b0,
      mult_b1 = mult_b1,
      mult_b2 = mult_b2
    ),
    coefficients = list(
      b0 = b0, b1 = b1, b2 = b2,
      a1 = a1, a2 = a2
    )
  )
  
  return(biquad_module)
}

# =============================================================================
# Register Module for Pipelined Molecular Computation
# =============================================================================
# This module implements register components that hold values between clock
# phases, enabling non-blocking assignments and staged computation.
#
# Registers capture input data on a clock edge and hold it until the next
# clock cycle, implementing the behavior of blocking assignments (<=) in
# hardware description languages.
# =============================================================================

jn <- function(...) { paste(..., sep = '') }

#' Make Register Module
#'
#' Creates a register that captures and holds a value based on a clock enable signal.
#' The register implements non-blocking assignment behavior, where the output
#' is updated only when the clock enable is active.
#'
#' This enables pipelined operations like:
#'   a = b + c      (combinational)
#'   reg_a <= a     (captured on clock edge)
#'   c <= reg_a     (updated on next clock edge)
#'
#' @param name Register name prefix
#' @param input_data Name of input data signal
#' @param clock_enable Name of clock enable signal (e.g., 'clk_pos')
#' @param output_data Name of output data signal
#' @param initial_value Initial register value (default: 0)
#' @param rate Reaction rate for register operations
#'
#' @return A gate object compatible with circuit_add_gate
#' @export
#'
#' @examples
#' reg <- Make_Register(
#'   name = 'reg_a',
#'   input_data = 'data_in',
#'   clock_enable = 'clk_pos',
#'   output_data = 'data_out',
#'   initial_value = 0,
#'   rate = 1e-3
#' )
Make_Register <- function(name,
                          input_data,
                          clock_enable,
                          output_data,
                          initial_value = 0,
                          rate = 1e-3) {
  
  # Internal register state
  internal_state <- jn(name, '_state')
  capture_signal <- jn(name, '_capture')
  
  species <- list(
    input = input_data,
    clock = clock_enable,
    internal = internal_state,
    capture = capture_signal,
    output = output_data
  )
  
  ci <- c(
    0,              # input (driven by previous stage)
    0,              # clock (driven by clock system)
    initial_value,  # internal state (initial register value)
    0,              # capture signal
    initial_value   # output (same as initial state)
  )
  
  reactions <- c(
    # Capture: when clock is high and input is present, generate capture signal
    # input + clock -> input + clock + capture
    jn(input_data, ' + ', clock_enable, ' -> ', 
       input_data, ' + ', clock_enable, ' + ', capture_signal),
    
    # Update internal state: capture signal transfers input to internal state
    # input + capture -> internal + capture
    jn(input_data, ' + ', capture_signal, ' -> ', 
       internal_state, ' + ', capture_signal),
    
    # Output follows internal state
    # internal -> internal + output
    jn(internal_state, ' -> ', internal_state, ' + ', output_data),
    
    # Decay of capture signal
    jn(capture_signal, ' -> 0'),
    
    # Output decay (to allow updates)
    jn(output_data, ' -> 0')
  )
  
  ki <- c(
    rate * 5,   # Capture generation
    rate * 10,  # State update
    rate * 10,  # Output generation
    rate * 2,   # Capture decay
    rate * 2    # Output decay
  )
  
  gate <- list(
    name = name,
    species = species,
    ci = ci,
    reactions = reactions,
    ki = ki
  )
  
  return(gate)
}

#' Make Simple Latch
#'
#' Creates a simple latch that passes input to output when enable is high.
#' Unlike a register, a latch is level-sensitive (not edge-triggered).
#'
#' @param name Latch name prefix
#' @param input_data Name of input signal
#' @param enable Name of enable signal
#' @param output_data Name of output signal
#' @param initial_value Initial latch value
#' @param rate Reaction rate
#'
#' @return A gate object
#' @export
Make_Latch <- function(name,
                      input_data,
                      enable,
                      output_data,
                      initial_value = 0,
                      rate = 1e-3) {
  
  species <- list(
    input = input_data,
    enable = enable,
    output = output_data
  )
  
  ci <- c(0, 0, initial_value)
  
  reactions <- c(
    # When enabled, output follows input
    jn(input_data, ' + ', enable, ' -> ', input_data, ' + ', enable, ' + ', output_data),
    # Output decay when not driven
    jn(output_data, ' -> 0')
  )
  
  ki <- c(rate * 10, rate * 2)
  
  gate <- list(
    name = name,
    species = species,
    ci = ci,
    reactions = reactions,
    ki = ki
  )
  
  return(gate)
}

#' Make D Flip-Flop
#'
#' Creates a D flip-flop that captures data on the rising edge of the clock.
#' This is a proper edge-triggered register.
#'
#' @param name Flip-flop name prefix
#' @param d_input D input signal name
#' @param clock Clock signal name
#' @param q_output Q output signal name
#' @param initial_value Initial Q value
#' @param rate Reaction rate
#'
#' @return A gate object
#' @export
Make_D_FlipFlop <- function(name,
                            d_input,
                            clock,
                            q_output,
                            initial_value = 0,
                            rate = 1e-3) {
  
  # Use the register implementation
  Make_Register(
    name = name,
    input_data = d_input,
    clock_enable = clock,
    output_data = q_output,
    initial_value = initial_value,
    rate = rate
  )
}

#' Make Pipeline Register Bank
#'
#' Creates a bank of registers for holding multiple signals in a pipeline stage.
#' All registers share the same clock enable.
#'
#' @param name Bank name prefix
#' @param input_signals Vector of input signal names
#' @param clock_enable Clock enable signal name
#' @param output_signals Vector of output signal names (must match input_signals length)
#' @param initial_values Vector of initial values (optional)
#' @param rate Reaction rate
#'
#' @return A list of gate objects, one per register
#' @export
Make_Register_Bank <- function(name,
                               input_signals,
                               clock_enable,
                               output_signals,
                               initial_values = NULL,
                               rate = 1e-3) {
  
  if (length(input_signals) != length(output_signals)) {
    stop("input_signals and output_signals must have the same length")
  }
  
  if (is.null(initial_values)) {
    initial_values <- rep(0, length(input_signals))
  }
  
  if (length(initial_values) != length(input_signals)) {
    stop("initial_values must have the same length as input_signals")
  }
  
  registers <- list()
  
  for (i in seq_along(input_signals)) {
    reg_name <- jn(name, '_reg', i)
    
    registers[[i]] <- Make_Register(
      name = reg_name,
      input_data = input_signals[i],
      clock_enable = clock_enable,
      output_data = output_signals[i],
      initial_value = initial_values[i],
      rate = rate
    )
  }
  
  return(registers)
}

#' Zero-Order Hold (ZOH)
#' 
#' Implements a zero-order hold for signal reconstruction
#' Uses a two-phase pipeline: SAMPLE -> HOLD
#' 
#' In digital-to-analog conversion, ZOH maintains (holds) the sampled value
#' constant between sampling instants. This is the most common reconstruction
#' method used in DACs.
#' 
#' Pipeline stages:
#' 1. SAMPLE phase: Sample the input signal into internal storage and output
#' 2. HOLD phase: Maintain the output value (storage persists with slow decay)
#' 
#' The two-phase design ensures that:
#' - Input sampling is synchronized
#' - Value is held stable between samples
#' - Output is updated in controlled manner
#' 
#' @param name Module name
#' @param nameInput Input signal species name
#' @param nameOutput Output signal species name (held value)
#' @param sample_period Sampling period (time between samples)
#' @param cinput Initial input concentration
#' @param crange Concentration range
#' @param rate Base reaction rate
#' @return Zero-order hold module structure
#' @examples 
#' # Create ZOH with 50 time unit sampling period
#' zoh <- Make_ZeroOrderHold('zoh1', 'X_continuous', 'X_discrete', 
#'                           sample_period = 50, cinput = 0, 
#'                           crange = 10, rate = 1e-3)
Make_ZeroOrderHold <- function(name, nameInput, nameOutput, sample_period,
                               cinput = 0, crange = 10, rate = 1e-3) {
  
  # Helper function for string concatenation
  jn <- function(...) { paste(..., sep = '') }
  
  # Create two-phase generator using proven Dalchau oscillator
  phase_gen <- Make_Phase_Generator(
    name = jn(name, '_phase_gen'),
    period = sample_period,
    crange = crange,
    rate = rate
  )
  
  # Use the two phases from the oscillator
  # phase1 (osc_A) = SAMPLE phase
  # phase2 (osc_B) = HOLD phase
  phase_sample <- phase_gen$species$phase1  # This is osc_A
  phase_hold <- phase_gen$species$phase2    # This is osc_B
  
  # Internal storage species
  storage <- jn(name, '_storage')
  
  # Species list combining phase generator and ZOH-specific species
  species_list <- c(
    nameInput,
    storage,
    nameOutput,
    unlist(phase_gen$species)
  )
  
  species_unique <- unique(species_list)
  
  # Initial concentrations from phase generator
  ci_phase <- phase_gen$ci
  
  # Additional initial concentrations for ZOH-specific species
  ci_zoh <- c(
    cinput,    # input
    0,         # storage (empty initially)
    0          # output (empty initially)
  )
  
  # Combine and handle duplicates
  ci_all <- c(ci_zoh, ci_phase)
  
  # Create proper ci vector for unique species
  ci_final <- c()
  for (sp in species_unique) {
    # Find indices of this species in all species
    idx_in_list <- which(species_list == sp)
    # Take the first non-zero concentration if any
    concs <- ci_all[idx_in_list]
    ci_final <- c(ci_final, ifelse(any(concs != 0), max(concs), 0))
  }
  
  # =========================================================================
  # Reactions for two-phase sample-and-hold
  # =========================================================================
  
  # Phase generator reactions (oscillator)
  reactions_phase <- phase_gen$reactions
  ki_phase <- phase_gen$ki
  
  # ZOH-specific reactions
  reactions_zoh <- c()
  ki_zoh <- c()
  
  # SAMPLE phase: When phase_sample (osc_A) is high, sample input to storage and output
  # Input + PhaseSample -> Storage + Output + PhaseSample (catalytic)
  reactions_zoh <- c(reactions_zoh,
    jn(nameInput, ' + ', phase_sample, ' -> ', 
       storage, ' + ', nameOutput, ' + ', phase_sample)
  )
  ki_zoh <- c(ki_zoh, rate * 5)  # Fast sampling
  
  # HOLD phase: Storage maintains value with slow decay
  # Output decay (moderate - will be regenerated from storage)
  reactions_zoh <- c(reactions_zoh,
    jn(nameOutput, ' -> 0')
  )
  ki_zoh <- c(ki_zoh, rate * 0.2)  # Slow decay to hold value
  
  # Storage decay (very slow - maintains value)
  reactions_zoh <- c(reactions_zoh,
    jn(storage, ' -> 0')
  )
  ki_zoh <- c(ki_zoh, rate * 0.05)  # Very slow decay
  
  # Storage continuously produces output (to maintain held value)
  reactions_zoh <- c(reactions_zoh,
    jn(storage, ' -> ', storage, ' + ', nameOutput)
  )
  ki_zoh <- c(ki_zoh, rate * 2)  # Moderate production
  
  # Combine all reactions
  all_reactions <- c(reactions_phase, reactions_zoh)
  all_ki <- c(ki_phase, ki_zoh)
  
  zoh_module <- list(
    name = name,
    species = species_unique,
    species_nested = list(
      input = nameInput,
      storage = storage,
      output = nameOutput,
      phase_sample = phase_sample,
      phase_hold = phase_hold,
      phase_gen = phase_gen$species
    ),
    reactions = all_reactions,
    ci = ci_final,
    ki = all_ki,
    sample_period = sample_period,
    components = list(
      phase_gen = phase_gen
    )
  )
  
  return(zoh_module)
}
