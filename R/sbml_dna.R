# Papel D: ponte entre o modelo SBML importado e a tradução para DNA.
#
# Fluxo: read_sbml_to_dnar() -> prepare_model_for_dna() -> dna_compatibility()
#        -> dna_crn_from_model() -> react_4domain().
#
# Depende de R/sbml.R (papel B) e R/sbml_kinetics.R (papel C):
# usa model$kinetic_info, sbml_split_reversible(), sbml_kinetic_report(),
# model$conversion_factors (dx_scale) e model$fixed_species.

# Nomes usados pelas espécies auxiliares de translate_4domain_crn() e
# get_buff_modules(). Espécies SBML com esses nomes são renomeadas.
DNA_RESERVED_PATTERN <- "^(L|H|W|O|T|G|LS|HS|WS)[0-9]*$"

#' Split a DNAr reaction string into reactant and product terms.
#'
#' @param reaction A reaction such as `"2A + B -> C"` or `"0 -> A"`.
#' @return A list with `reactants` and `products`, each a data.frame with
#'   columns `species` and `coef`.
#' @keywords internal
dna_parse_reaction <- function(reaction) {
    sides <- strsplit(reaction, "->", fixed = TRUE)[[1]]
    if (length(sides) != 2L) {
        stop("reação sem '->' ou com mais de uma seta: ", reaction)
    }
    parse_side <- function(side) {
        side <- trimws(side)
        if (side %in% c("", "0")) {
            return(data.frame(species = character(), coef = numeric(),
                              stringsAsFactors = FALSE))
        }
        terms <- trimws(strsplit(side, "+", fixed = TRUE)[[1]])
        coef_text <- regmatches(terms, regexpr("^[0-9.]*", terms))
        species <- trimws(substring(terms, nchar(coef_text) + 1L))
        coef <- ifelse(nzchar(coef_text), suppressWarnings(as.numeric(coef_text)), 1)
        data.frame(species = species, coef = coef, stringsAsFactors = FALSE)
    }
    list(reactants = parse_side(sides[[1]]), products = parse_side(sides[[2]]))
}

#' Format the output of [dna_parse_reaction()] back into a reaction string.
#'
#' @param parsed Output of [dna_parse_reaction()].
#' @param rename Optional named character vector `old -> new` species names.
#' @return A reaction string.
#' @keywords internal
dna_format_reaction <- function(parsed, rename = NULL) {
    format_side <- function(side) {
        if (nrow(side) == 0L) {
            return("0")
        }
        names_out <- side$species
        if (!is.null(rename)) {
            hit <- names_out %in% names(rename)
            names_out[hit] <- rename[names_out[hit]]
        }
        coef <- ifelse(side$coef == 1, "",
                       format(side$coef, trim = TRUE, scientific = FALSE))
        paste(paste0(coef, names_out), collapse = " + ")
    }
    paste(format_side(parsed$reactants), "->", format_side(parsed$products))
}

#' Prepare an imported SBML model for the DNA translation.
#'
#' Splits reversible mass-action reactions into forward and backward steps
#' and attaches the support table of the papel C detector, recomputed on the
#' split model so that each elementary step is checked.
#'
#' @param model Object returned by `read_sbml_to_dnar()`.
#' @return The split model with an extra `dna_report` data.frame.
#' @keywords internal
prepare_model_for_dna <- function(model) {
    if (is.null(model$kinetic_info)) {
        stop("modelo sem kinetic_info: importe com read_sbml_to_dnar() ",
             "carregando R/sbml_kinetics.R")
    }
    split <- sbml_split_reversible(model)
    split$dna_report <- sbml_kinetic_report(split)
    split
}

#' Common derivative scale of the non-fixed species, if there is one.
#'
#' SBML kinetic laws are amount/time and `react4()` converts them with a
#' per-species `dx_scale`. A single mass-action constant can only absorb that
#' scale when it is the same for every species that changes.
#'
#' @param model A model from [prepare_model_for_dna()].
#' @return The common scale, or `NA` when species use different scales.
#' @keywords internal
dna_uniform_scale <- function(model) {
    scale <- model$conversion_factors
    if (is.null(scale) || length(scale) == 0L) {
        return(1)
    }
    free <- setdiff(names(scale), model$fixed_species)
    if (length(free) == 0L) {
        return(1)
    }
    scale <- as.numeric(scale[free])
    reference <- scale[[1]]
    if (!is.finite(reference) || reference <= 0) {
        return(NA_real_)
    }
    if (all(abs(scale - reference) <= 1e-12 * abs(reference))) reference else NA_real_
}

#' Check whether a prepared SBML model can be translated to DNA.
#'
#' A model is compatible when every reaction is mass action with at most two
#' reactants (papel C detector), all stoichiometries are integers, no reaction
#' is a formation `0 -> X` (known failure in the 4-domain buffers), and all
#' changing species share the same `dx_scale`.
#'
#' @param model A model from [prepare_model_for_dna()].
#' @return A list with `ok` (logical) and `reasons` (character).
#' @keywords internal
dna_compatibility <- function(model) {
    reasons <- character()
    if (length(model$species) == 0L || length(model$reactions) == 0L) {
        return(list(ok = FALSE, reasons = "modelo sem reações"))
    }

    report <- model$dna_report
    if (is.null(report)) {
        stop("modelo sem dna_report: use prepare_model_for_dna() antes")
    }
    bad <- which(!report$dna_ok)
    if (length(bad) > 0L) {
        reasons <- c(reasons, sprintf("reação %s: %s", report$id[bad], report$reason[bad]))
    }

    for (i in seq_along(model$reactions)) {
        parsed <- dna_parse_reaction(model$reactions[[i]])
        coef <- c(parsed$reactants$coef, parsed$products$coef)
        if (any(!is.finite(coef)) || any(abs(coef - round(coef)) > 1e-12)) {
            reasons <- c(reasons, sprintf(
                "reação %s: estequiometria não inteira", report$id[[i]]))
        }
        if (nrow(parsed$reactants) == 0L) {
            reasons <- c(reasons, sprintf(
                "reação %s: reação de formação (0 -> X) ainda não suportada na tradução 4-domain",
                report$id[[i]]))
        }
    }

    if (is.na(dna_uniform_scale(model))) {
        reasons <- c(reasons, paste(
            "dx_scale diferente entre espécies (volumes ou fatores de",
            "conversão distintos): nenhum ki único reproduz a dinâmica"))
    }

    reasons <- unique(reasons)
    list(ok = length(reasons) == 0L, reasons = reasons)
}

#' Build the mass-action CRN that goes into the DNA translation.
#'
#' Uses the constant `k` found by the papel C detector (not `model$ki`),
#' multiplies it by the common `dx_scale`, and renames species whose names
#' collide with the auxiliary species of the translation.
#'
#' @param model A compatible model from [prepare_model_for_dna()].
#' @return A list with `species`, `ci`, `reactions`, `ki`, `fixed_species`
#'   (renamed) and `name_map` (new name -> original SBML id).
#' @keywords internal
dna_crn_from_model <- function(model) {
    compat <- dna_compatibility(model)
    if (!compat$ok) {
        stop("modelo incompatível com DNA: ", paste(compat$reasons, collapse = "; "))
    }
    k <- vapply(model$kinetic_info, function(info) as.numeric(info$k), numeric(1))
    ki <- k * dna_uniform_scale(model)

    original <- model$species
    renamed <- ifelse(grepl(DNA_RESERVED_PATTERN, original),
                      paste0("sp_", original), original)
    rename <- setNames(renamed, original)
    reactions <- vapply(model$reactions, function(reaction) {
        dna_format_reaction(dna_parse_reaction(reaction), rename)
    }, character(1), USE.NAMES = FALSE)

    fixed <- model$fixed_species
    if (is.null(fixed)) {
        fixed <- character()
    }
    list(
        species = unname(renamed),
        ci = unname(as.numeric(model$ci[original])),
        reactions = reactions,
        ki = unname(ki),
        fixed_species = unname(rename[fixed]),
        name_map = setNames(original, renamed)
    )
}
