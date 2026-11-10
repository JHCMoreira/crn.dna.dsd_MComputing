# DNAr is a program used to simulate formal Chemical Reaction Networks
# and the ones based on DNA.
# Copyright (C) 2017  Daniel Kneipp <danielv[at]dcc[dot]ufmg[dot]com[dot]br>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#
# 4-domain DNA approach
#

#' Check if the reaction is compatible with the ones
#' supported by \code{\link{react_4domain}()}.
#'
#' The reactions supported by \code{\link{react_4domain}()} are
#' only the uni or bimolecular. This function checks if the
#' \code{reaction} parameter meets this requirement, returning
#' \code{TRUE} if the reaction is supported, or \code{FALSE}
#' otherwise.
#'
#' @examples
#' DNAr:::check_reaction_4domain('A + B -> C')   # Should return TRUE
#' DNAr:::check_reaction_4domain('2A -> B')      # Should return TRUE
#' DNAr:::check_reaction_4domain('2A + B -> C')  # Should return FALSE
check_reaction_4domain <- function(reaction) {
    left_part <- get_first_part(reaction)
    sto <- get_stoichiometry_part(left_part)
    return(sto < 3)
}

#' Get buffer modules
#'
#' This function is used to add buffer modules according
#' to the theory described by Soloveichik D et al. `[1]`.
#' The parameters of this function follows the same
#' semantics of \code{\link{react_4domain}()}.
#'
#' @return \code{NULL} if no buffer modules were added. Otherwise
#' it returns a list with:
#'   - `lambda_1`      = lambda^{-1} value;
#'   - `new_species`   = vector with the new species added;
#'   - `new_cis`       = vector with the initial concentrations;
#'   - `new_reactions` = vector with the new reactions;
#'   - `new_ks`        = constant rate of the new reactions.
#'
#' @references
#'   - `[1]` \insertRef{soloveichik2010dna}{DNAr}
get_buff_modules <- function(reactions, ki, qmax, cmax) {
    sigmas <- list()
    bff_aux_species <- c('LS', 'HS', 'WS')
    #             LS   HS  WS
    bff_cis <- c(cmax, 0, cmax)

    # Calculate the sigma for each species.
    # Formation reactions (0 -> A) are ignored
    uni_count <- 0
    for(i in seq_along(reactions)) {
        reactants <- get_reactants(reactions[i])
        first_reactant <- reactants[[1]]

        if(is_bimolecular(reactions[[i]])) {
            if(is.null(sigmas[[first_reactant]])) {
                sigmas[first_reactant] <- ki[[i]]
            } else {
                sigmas[first_reactant] <- sigmas[[first_reactant]] + ki[[i]]
            }
        } else if(!is_formation(reactions[[i]])) {
            uni_count <- uni_count + 1
            if(is.null(sigmas[[first_reactant]])) {
                sigmas[first_reactant] <- 0
            }
        }
    }

    # There is no sigma (there is only formation or unimolecular reactions)
    if(length(sigmas) == 0) {
        return(NULL)
    }

    reaction_sigma <- max(unlist(sigmas))
    lambda_1 <- qmax / (qmax - reaction_sigma)

    new_ks <- c()
    new_bff_reactions <- c()
    new_species <- c()
    new_cis <- c()
    for(i in seq_along(sigmas)) {
        if(!(sigmas[[i]] == reaction_sigma)) {
            qs <- lambda_1 * (reaction_sigma - sigmas[[i]])
            aux_specs <- paste(bff_aux_species, as.character(i), sep = '')
            input_spec <- names(sigmas)[i]

            forward_reaction <- paste(input_spec, '+', aux_specs[1], '->',
                                      aux_specs[2], '+', aux_specs[3])
            backward_reaction <- paste(aux_specs[2], '+', aux_specs[3], '->',
                                       input_spec, '+', aux_specs[1])

            new_bff_reactions <- c(new_bff_reactions,
                                   forward_reaction,
                                   backward_reaction)
            new_ks <- c(new_ks, qs, qmax)
            new_species <- c(new_species, aux_specs[1],
                             aux_specs[2],
                             aux_specs[3])
            new_cis <- c(new_cis, bff_cis)
        }
    }

    # If even with bimolecular reactions, there is no buffer module
    # to add.
    if(length(new_ks) == 0) {
        return(NULL)
    } else {
        ret <- list(
            lambda_1 = lambda_1,
            new_species = new_species,
            new_cis = new_cis,
            new_reactions = new_bff_reactions,
            new_ks = new_ks
        )
        return(ret)
    }
}

normalize_dna_kinetics <- function(dna_kinetics, reactions) {
    defaults <- list(
        domain_lengths = data.frame(
            domain_id = paste0('domain_', seq_along(reactions)),
            num_bases = rep(6, length(reactions)),
            stringsAsFactors = FALSE
        ),
        sequence_composition = NULL,
        gamma = 0.5,
        k_step = 1.4e6,
        RT = 0.593,
        DG_per_BP = -1.2
    )

    if(is.null(dna_kinetics)) {
        return(defaults)
    }

    for(name in intersect(names(dna_kinetics), names(defaults))) {
        defaults[[name]] <- dna_kinetics[[name]]
    }

    return(defaults)
}

scale_dna_kinetics_rate <- function(reaction, reaction_index, base_rate, dna_kinetics, is_bimolecular_reaction) {
    domain_id <- paste0('domain_', reaction_index)
    domain_lengths <- dna_kinetics$domain_lengths
    L <- 6

    if(!is.null(domain_lengths) && all(c('domain_id', 'num_bases') %in% names(domain_lengths))) {
        if(domain_id %in% domain_lengths$domain_id) {
            L <- domain_lengths$num_bases[domain_lengths$domain_id == domain_id][1]
        }
    }

    if(is_bimolecular_reaction) {
        return(base_rate / (L^dna_kinetics$gamma))
    }

    if(grepl('migration', reaction, ignore.case = TRUE)) {
        return(dna_kinetics$k_step / (L^2))
    }

    sequence_composition <- dna_kinetics$sequence_composition
    if(!is.null(sequence_composition) && all(c('domain_id', 'gc') %in% names(sequence_composition))) {
        if(domain_id %in% sequence_composition$domain_id) {
            gc_count <- sequence_composition$gc[sequence_composition$domain_id == domain_id][1]
            at_count <- L - gc_count
            delta_G <- (gc_count * -2.0) + (at_count * -1.0)
        } else {
            delta_G <- L * dna_kinetics$DG_per_BP
        }
    } else {
        delta_G <- L * dna_kinetics$DG_per_BP
    }

    return(base_rate * exp(-delta_G / dna_kinetics$RT))
}

translate_4domain_crn <- function(
    species,
    ci,
    reactions,
    ki,
    qmax,
    cmax,
    alpha,
    beta,
    t,
    auto_buffer = TRUE,
    dna_kinetics = NULL,
    apply_dna_kinetics = FALSE
) {
    reactions <- check_crn(species, ci, reactions, ki, t)

    for(r in reactions) {
        if(!check_reaction_4domain(r)) {
            stop(paste('Failed to process reaction', r))
        }
    }

    dna_kinetics <- normalize_dna_kinetics(dna_kinetics, reactions)

    new_reactions <- c()
    new_species <- c(species)
    new_ks <- c()
    new_cis <- c(ci * beta)

    uni_aux_species <- c('G', 'O', 'T')
    bi_aux_species <- c('L', 'H', 'W', 'O', 'T')

    scaled_ks <- c()
    is_bimolecular_map <- c()
    for(i in seq_along(reactions)) {
        if(is_bimolecular(reactions[i])) {
            scaled_ks[i] <- ki[i] / (alpha * beta)
            is_bimolecular_map[i] <- TRUE
        } else {
            scaled_ks[i] <- ki[i] / alpha
            is_bimolecular_map[i] <- FALSE
        }
    }

    buffer_stuff <- NULL
    if(auto_buffer) {
        buffer_stuff <- get_buff_modules(reactions, scaled_ks, qmax, cmax)
    }

    if(!is.null(buffer_stuff)) {
        for(i in seq_along(new_cis)) {
            new_cis[i] <- new_cis[[i]] * buffer_stuff$lambda_1
        }
    }

    for(i in seq_along(reactions)) {
        new_reactions_for_i <- c()
        new_species_for_i <- c()
        new_ks_for_i <- c()
        new_cis_for_i <- c()

        if(is_bimolecular_map[i]) {
            aux <- paste(bi_aux_species, as.character(i), sep = '')

            left_part <- get_first_part(reactions[i])
            l_p_specs <- get_species(left_part)

            if(length(l_p_specs) == 1) {
                l_p_specs <- c(l_p_specs, l_p_specs)
            }

            right_part <- get_second_part(reactions[i])

            if(apply_dna_kinetics) {
                qi_with_buff <- scale_dna_kinetics_rate(
                    reaction = reactions[i],
                    reaction_index = i,
                    base_rate = scaled_ks[i],
                    dna_kinetics = dna_kinetics,
                    is_bimolecular_reaction = TRUE
                )
            } else {
                qi_with_buff <- scaled_ks[i]
            }

            new_reactions_for_i <- c(paste(l_p_specs[1], '+', aux[1], '->',
                                           aux[2], '+', aux[3]),
                                     paste(aux[2], '+', aux[3], '->',
                                           l_p_specs[1], '+', aux[1]),
                                     paste(l_p_specs[2], '+', aux[2], '->',
                                           aux[4]))

            new_species_for_i <- c(aux[1], aux[2], aux[3], aux[4])
            new_cis_for_i <- c(cmax, 0.0, cmax, 0.0)
            new_ks_for_i <- c(qi_with_buff, qmax, qmax)

            if(!isempty_part(right_part)) {
                new_reactions_for_i <- c(new_reactions_for_i,
                                         paste(aux[4], '+', aux[5], '->',
                                               right_part))
                new_species_for_i <- c(new_species_for_i, aux[5])
                new_cis_for_i <- c(new_cis_for_i, cmax)
                new_ks_for_i <- c(new_ks_for_i, qmax)
            }
        } else {
            aux <- paste(uni_aux_species, as.character(i), sep = '')

            left_part <- get_first_part(reactions[i])
            l_p_specs <- get_species(left_part)
            right_part <- get_second_part(reactions[i])

            if(is_formation(reactions[[i]])) {
                new_reactions_for_i <- c(paste(aux[1], '->', aux[2]))
            } else {
                new_reactions_for_i <- c(paste(l_p_specs, '+', aux[1], '->',
                                               aux[2]))
            }
            new_species_for_i <- c(aux[1], aux[2])
            new_cis_for_i <- c(cmax, 0.0)

            if(apply_dna_kinetics) {
                qi_with_buff <- scale_dna_kinetics_rate(
                    reaction = reactions[i],
                    reaction_index = i,
                    base_rate = scaled_ks[i] / cmax,
                    dna_kinetics = dna_kinetics,
                    is_bimolecular_reaction = FALSE
                )
            } else {
                qi_with_buff <- scaled_ks[i] / cmax
            }

            new_ks_for_i <- c(qi_with_buff)

            if(!isempty_part(right_part)) {
                new_reactions_for_i <- c(new_reactions_for_i,
                                         paste(aux[2], '+', aux[3], '->',
                                               right_part))
                new_species_for_i <- c(new_species_for_i, aux[3])
                new_cis_for_i <- c(new_cis_for_i, cmax)
                new_ks_for_i <- c(new_ks_for_i, qmax)
            }
        }

        new_species <- c(new_species, new_species_for_i)
        new_reactions <- c(new_reactions, new_reactions_for_i)
        new_ks <- c(new_ks, new_ks_for_i)
        new_cis <- c(new_cis, new_cis_for_i)
    }

    if(!is.null(buffer_stuff)) {
        new_species <- c(new_species, buffer_stuff$new_species)
        new_reactions <- c(new_reactions, buffer_stuff$new_reactions)
        new_ks <- c(new_ks, buffer_stuff$new_ks)
        new_cis <- c(new_cis, buffer_stuff$new_cis)
    }

    return(list(
        species = new_species,
        ci = new_cis,
        reactions = new_reactions,
        ki = new_ks
    ))
}

#' Translate a formal CRN into a set of DNA based reactions according to the
#' approach described by Soloveichik D et al. `[1]`
#'
#' This function is used to simulate a chemical circuit based on DNA made
#' to behavior as expected from a CRN specified by the parameters. In another
#' words, given the CRN^* passed through the parameters, another \eqn{CRN_2} is
#' created based on reactions between strands of DNA. CRN_2 is simulated using
#' \code{\link{react}()}. The matrix behavior of CRN_2 is returned but only
#' of the species specified in the \code{species} parameter, the behavior of
#' the  auxiliary ones are not returned. The parameters of these functions
#' follows the same pattern of \code{\link{react}()}, only with some additions
#' required by this approach `[1]` (here named as 4-domain).
#'
#' @section Known limitations:
#'   - It only support uni or bimolecular reactions;
#'   - Because of \code{\link{react}()} known limitation, this function also
#'   doesn't support bidirectional reactions;
#'   - The species names `L`, `H`, `W`, `O`, `T`, `G`, `LS`,
#'   `HS`, `WS` with or without numbers after it are not supported because
#'   these are the reserved for the auxiliary ones. Ex.: `L2` and `LS2`
#'   are not supported but `LT` and `LT2` are.
#'
#' @param species      A vector with the species of the reaction. The order of
#'                     this vector is important because it will define the
#'                     column order of the returned behavior. The species names
#'                     `L[0-9]*`, `H[0-9]*`, `W[0-9]*`, `O[0-9]*`,
#'                     `T[0-9]*`, `G[0-9]*`, `LS[0-9]*`, `HS[0-9]*`,
#'                     `WS[0-9]*` are not supported. For more information
#'                     about this, see the Section of **Known limitations**.
#' @param ci           A vector specifying the initial concentrations of the
#'                     \code{species} specified, in order.
#' @param reactions    A vector with the reactions of the CRN^*.
#' @param ki           A vector defining the constant rate of each reaction
#'                     in \code{reactions}, in order.
#' @param qmax         Maximum rate constant for the auxiliary reactions.
#' @param cmax         Maximum initial concentration for the auxiliary species.
#' @param alpha,beta   Rescaling parameters.
#' @param t            A vector specifying the time interval. Each value
#'                     would be a specific time point.
#' @param auto_buffer  With the default value of `TRUE`, this specifies if
#'                     buffer modules should be generated automatically.
#' @param verbose      Be verbose and print information about the integration
#'                     process with `deSolve::diagnostics.deSolve()`. Default
#'                     value is `FALSE`
#' @param ...          Parameters passed to `deSolve::ode()`.
#'
#' @return A list with the attributes `behavior`, `species`, `ci`, `reactions`
#' and `ki`. These attributes are:
#'   - `behavior`: A matrix with each line being a specific point in the time
#'                 and each column but the first being the concentration of a
#'                 species. The first column is the time interval. The
#'                 behavior of the auxiliary species are also returned.
#'  - `species`  : A vector with all the species used in the reactions.
#'  - `ci`       : The initial concentration of each species.
#'  - `reactions`: All the reactions computed, including the ones generated
#'                 according to the 4-domain approach.
#'  - `ki`       : The rate constants of the reactions.
#'
#' @export
#'
#' @references
#'   - `[1]` \insertRef{soloveichik2010dna}{DNAr}
#'
#' @example demo/main_4domain.R
react_4domain <- function(
    species,
    ci,
    reactions,
    ki,
    qmax,
    cmax,
    alpha,
    beta,
    t,
    auto_buffer = TRUE,
    dna_kinetics = NULL,
    verbose = FALSE,
    forced_concentrations = NULL,
    engine = 'desolve',
    ...
) {
    translated <- translate_4domain_crn(
        species = species,
        ci = ci,
        reactions = reactions,
        ki = ki,
        qmax = qmax,
        cmax = cmax,
        alpha = alpha,
        beta = beta,
        t = t,
        auto_buffer = auto_buffer,
        dna_kinetics = dna_kinetics,
        apply_dna_kinetics = FALSE
    )

    # Run the reaction
    b <- react4(
        species   = translated$species,
        ci        = translated$ci,
        reactions = translated$reactions,
        ki        = translated$ki,
        t         = t,
        verbose   = verbose,
        forced_concentrations = forced_concentrations,
        engine = engine,
        ...
    )

    # Arrange the data to be returned in a list
    result <- list(
        behavior  = b,
        species   = translated$species,
        ci        = translated$ci,
        reactions = translated$reactions,
        ki        = translated$ki
    )

    # Return the behavior of all species (including the auxiliary ones),
    # initial concentrations, reactions and rate constants.
    return(result)
}

#' Translate a CRN to 4-domain DNA and simulate it stochastically
#'
#' This is the stochastic counterpart to \\code{\\link{react_4domain}()}.
#' The CRN is translated to the 4-domain DNA reaction network and then
#' simulated with \\code{\\link{react_stochastic}()} using DNA-specific
#' kinetic modifiers.
#'
#' @param dna_kinetics Optional list with DNA kinetic modifiers. Supported
#'   entries are `domain_lengths`, `sequence_composition`, `gamma`, `k_step`,
#'   `RT`, and `DG_per_BP`.
#'
#' @export
react_4domain_stochastic <- function(
    species,
    ci,
    reactions,
    ki,
    qmax,
    cmax,
    alpha,
    beta,
    t,
    auto_buffer = TRUE,
    dna_kinetics = NULL,
    seed = NULL,
    verbose = FALSE,
    ...
) {
    translated <- translate_4domain_crn(
        species = species,
        ci = ci,
        reactions = reactions,
        ki = ki,
        qmax = qmax,
        cmax = cmax,
        alpha = alpha,
        beta = beta,
        t = t,
        auto_buffer = auto_buffer,
        dna_kinetics = dna_kinetics,
        apply_dna_kinetics = TRUE
    )

    behavior <- react_stochastic(
        species = translated$species,
        ci = translated$ci,
        reactions = translated$reactions,
        ki = translated$ki,
        t = t,
        seed = seed,
        verbose = verbose,
        ...
    )

    result <- list(
        behavior = behavior,
        species = translated$species,
        ci = translated$ci,
        reactions = translated$reactions,
        ki = translated$ki
    )

    return(result)
}

#' Automatically modify a CRN to make it compatible with the DNA simulation
#'
#' This function automatically set the `qmax`, `cmax`, `alpha`, `beta`, and
#' optionally the parameters `method` and `t` (duration time).
#'
#' @param crn     The CRN to be modified.
#' @param t       Optional parameter used to redefine the duration time.
#' @param method  Optional parameter which sets the solver method.
#'
#' @return  The modified CRN with the nwe parameters, ready to be used with
#'          the function \code{\link{react_4domain}()}.
#'
#' @export
update_crn_4domain <- function(crn, t = NULL, method = NULL) {
    crn$qmax <- max(crn$ki) * 1e4
    crn$cmax <- max(crn$ci) * 1e4
    crn$alpha <- 1
    crn$beta <- 1
    if(!is.null(method)){
        crn$method <- method
    }
    if(!is.null(t)) {
        crn$t <- t
    }

    return(crn)
}
