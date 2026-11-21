#' Check whether the libSBML bridge is available.
#'
#' The Docker development image installs libSBML through its Python bindings.
#' Keeping the check here makes failures explicit for users running DNAr
#' outside that image.
#'
#' @return A logical scalar.
#' @export
sbml_available <- function() {
    nzchar(Sys.which("python3")) &&
        isTRUE(system2(
            "python3",
            c("-c", shQuote("import libsbml")),
            stdout = FALSE,
            stderr = FALSE
        ) == 0)
}

sbml_python_script <- function() {
    script <- tempfile(fileext = ".py")
    writeLines(c(
        "import json",
        "import math",
        "import sys",
        "import libsbml",
        "",
        "path = sys.argv[1]",
        "reader = libsbml.SBMLReader()",
        "document = reader.readSBML(path)",
        "model = document.getModel()",
        "if model is None:",
        "    raise RuntimeError('SBML document does not contain a model')",
        "",
        "def number(value):",
        "    if value is None:",
        "        return None",
        "    value = float(value)",
        "    if math.isnan(value):",
        "        return 'NaN'",
        "    if math.isinf(value):",
        "        return 'Inf' if value > 0 else '-Inf'",
        "    return value",
        "",
        "def optional_id(obj):",
        "    return obj.getId() if obj.isSetId() else None",
        "",
        "def math_formula(obj):",
        "    math = obj.getMath()",
        "    return None if math is None else libsbml.formulaToL3String(math)",
        "",
        "compartments = {}",
        "for i in range(model.getNumCompartments()):",
        "    compartment = model.getCompartment(i)",
        "    size = compartment.getSize() if compartment.isSetSize() else None",
        "    compartments[compartment.getId()] = {",
        "        'size': number(size),",
        "        'spatialDimensions': (",
        "            number(compartment.getSpatialDimensions())",
        "            if compartment.isSetSpatialDimensions() else None",
        "        )",
        "    }",
        "",
        "parameters = {}",
        "for i in range(model.getNumParameters()):",
        "    parameter = model.getParameter(i)",
        "    if parameter.isSetValue():",
        "        parameters[parameter.getId()] = number(parameter.getValue())",
        "",
        "species = []",
        "for i in range(model.getNumSpecies()):",
        "    item = model.getSpecies(i)",
        "    compartment_id = item.getCompartment()",
        "    # Keep raw amount/concentration; R converts around initialAssignment",
        "    # so compartment size from IA is honored and hosu*concentration",
        "    # is multiplied by volume.",
        "    initial = None",
        "    initial_type = None",
        "    if item.isSetInitialConcentration():",
        "        initial = item.getInitialConcentration()",
        "        initial_type = 'concentration'",
        "    elif item.isSetInitialAmount():",
        "        initial = item.getInitialAmount()",
        "        initial_type = 'amount'",
        "    if initial is None:",
        "        initial = 0.0",
        "    species.append({",
        "        'id': item.getId(),",
        "        'initial': number(initial),",
        "        'initialType': initial_type,",
        "        'compartment': compartment_id,",
        "        'boundaryCondition': item.getBoundaryCondition(),",
        "        'constant': item.getConstant(),",
        "        'hasOnlySubstanceUnits': item.getHasOnlySubstanceUnits(),",
        "        'conversionFactor': (",
        "            item.getConversionFactor() if item.isSetConversionFactor() else None",
        "        )",
        "    })",
        "",
        "def species_reference(ref):",
        "    return {",
        "        'id': optional_id(ref),",
        "        'species': ref.getSpecies(),",
        "        'stoich': number(ref.getStoichiometry())",
        "    }",
        "",
        "reactions = []",
        "for i in range(model.getNumReactions()):",
        "    reaction = model.getReaction(i)",
        "    reactants = [species_reference(reaction.getReactant(j))",
        "                 for j in range(reaction.getNumReactants())]",
        "    products = [species_reference(reaction.getProduct(j))",
        "                for j in range(reaction.getNumProducts())]",
        "    formula = None",
        "    local_parameters = {}",
        "    law = reaction.getKineticLaw()",
        "    if law is not None:",
        "        formula = libsbml.formulaToL3String(law.getMath())",
        "        for j in range(law.getNumLocalParameters()):",
        "            parameter = law.getLocalParameter(j)",
        "            if parameter.isSetValue():",
        "                local_parameters[parameter.getId()] = number(parameter.getValue())",
        "        for j in range(law.getNumParameters()):",
        "            parameter = law.getParameter(j)",
        "            if parameter.isSetValue():",
        "                local_parameters[parameter.getId()] = number(parameter.getValue())",
        "    reactions.append({'id': reaction.getId(), 'reversible': reaction.getReversible(),",
        "                      'reactants': reactants,",
        "                      'products': products, 'formula': formula,",
        "                      'parameters': local_parameters})",
        "",
        "initial_assignments = []",
        "for i in range(model.getNumInitialAssignments()):",
        "    assignment = model.getInitialAssignment(i)",
        "    initial_assignments.append({",
        "        'symbol': assignment.getSymbol(),",
        "        'formula': math_formula(assignment)",
        "    })",
        "",
        "rules = []",
        "for i in range(model.getNumRules()):",
        "    rule = model.getRule(i)",
        "    type_code = rule.getTypeCode()",
        "    if type_code == libsbml.SBML_ASSIGNMENT_RULE:",
        "        rule_type = 'assignment'",
        "        variable = rule.getVariable()",
        "    elif type_code == libsbml.SBML_RATE_RULE:",
        "        rule_type = 'rate'",
        "        variable = rule.getVariable()",
        "    elif type_code == libsbml.SBML_ALGEBRAIC_RULE:",
        "        rule_type = 'algebraic'",
        "        variable = None",
        "    else:",
        "        rule_type = 'unknown'",
        "        variable = None",
        "    rules.append({",
        "        'type': rule_type,",
        "        'variable': variable,",
        "        'formula': math_formula(rule)",
        "    })",
        "",
        "errors = []",
        "warnings = []",
        "for i in range(document.getNumErrors()):",
        "    error = document.getError(i)",
        "    entry = {'line': error.getLine(), 'message': error.getMessage()}",
        "    if error.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:",
        "        errors.append(entry)",
        "    else:",
        "        warnings.append(entry)",
        "",
        "print(json.dumps({",
        "    'compartments': compartments,",
        "    'species': species,",
        "    'parameters': parameters,",
        "    'reactions': reactions,",
        "    'initialAssignments': initial_assignments,",
        "    'rules': rules,",
        "    'conversionFactor': (",
        "        model.getConversionFactor() if model.isSetConversionFactor() else None",
        "    ),",
        "    'errors': errors,",
        "    'warnings': warnings",
        "}))"
    ), script)
    script
}

sbml_decode_number <- function(value) {
    if (is.null(value) || length(value) != 1L) {
        return(value)
    }
    if (is.character(value)) {
        if (identical(value, "Inf")) {
            return(Inf)
        }
        if (identical(value, "-Inf")) {
            return(-Inf)
        }
        if (identical(value, "NaN")) {
            return(NaN)
        }
    }
    value
}

sbml_decode_numbers <- function(node) {
    if (is.list(node)) {
        return(lapply(node, sbml_decode_numbers))
    }
    sbml_decode_number(node)
}

sbml_run_bridge <- function(xml_file_path) {
    if (!is.character(xml_file_path) || length(xml_file_path) != 1L ||
        !nzchar(xml_file_path)) {
        stop("'xml_file_path' must be a single non-empty path.")
    }
    if (!file.exists(xml_file_path)) {
        stop("SBML file does not exist: ", xml_file_path)
    }
    if (!sbml_available()) {
        stop(
            "libSBML is unavailable. Install python-libsbml and ensure ",
            "'python3' is on PATH."
        )
    }

    script <- sbml_python_script()
    on.exit(unlink(script), add = TRUE)
    output <- system2(
        "python3",
        c(script, normalizePath(xml_file_path, mustWork = TRUE)),
        stdout = TRUE,
        stderr = TRUE
    )
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
        stop("libSBML could not read the file: ", paste(output, collapse = "\n"))
    }
    parsed <- jsonlite::fromJSON(paste(output, collapse = ""), simplifyVector = FALSE)
    sbml_decode_numbers(parsed)
}

#' Validate an SBML document using libSBML.
#'
#' @param xml_file_path Path to an SBML document.
#' @return A list with `valid`, `errors`, and `warnings`.
#' @export
validate_sbml <- function(xml_file_path) {
    parsed <- sbml_run_bridge(xml_file_path)
    list(
        valid = length(parsed$errors) == 0L,
        errors = parsed$errors,
        warnings = parsed$warnings
    )
}

sbml_reaction_side <- function(references) {
    if (length(references) == 0L) {
        return("0")
    }
    parts <- vapply(references, function(reference) {
        stoich <- as.numeric(reference$stoich)
        coefficient <- if (isTRUE(all.equal(stoich, 1))) {
            ""
        } else {
            format(stoich, trim = TRUE, scientific = FALSE)
        }
        paste0(coefficient, reference$species)
    }, character(1))
    paste(parts, collapse = " + ")
}

#' Build an evaluation environment for SBML MathML/L3 formulas.
#'
#' Maps common `formulaToL3String` identifiers onto R functions so kinetic
#' laws can be evaluated with `eval(parse(text = formula), envir = env)`
#' without textual rewriting.
#'
#' @param values Named list or vector of symbol values.
#' @return An environment suitable for `eval()`.
#' @export
sbml_math_env <- function(values) {
    env <- new.env(parent = baseenv())
    env$ceil <- ceiling
    env$ln <- log
    env$exponentiale <- exp(1)
    env$true <- TRUE
    env$false <- FALSE
    env$piecewise <- function(...) {
        args <- list(...)
        n <- length(args)
        if (n == 0L) {
            return(NA_real_)
        }
        # SBML piecewise(value1, cond1, ..., valueN, condN[, otherwise]).
        # A lone otherwise argument is valid and must not call seq(1, 0, by=2).
        if (n == 1L) {
            return(args[[1L]])
        }
        last_piece <- if (n %% 2L == 0L) n else (n - 1L)
        for (i in seq.int(1L, last_piece, by = 2L)) {
            if (isTRUE(args[[i + 1L]])) {
                return(args[[i]])
            }
        }
        if (n %% 2L == 1L) args[[n]] else NA_real_
    }
    list2env(as.list(values), envir = env)
    env
}

sbml_compartment_size <- function(parsed, compartment_id) {
    if (is.null(compartment_id) || !nzchar(compartment_id)) {
        return(NA_real_)
    }
    compartment <- parsed$compartments[[compartment_id]]
    if (is.null(compartment) || is.null(compartment$size)) {
        return(NA_real_)
    }
    as.numeric(compartment$size)
}

#' Collect SBML name -> value bindings used when evaluating kinetic laws.
#'
#' Includes compartment sizes, parameters, species raw initials (amount or
#' concentration as declared; unit conversion is applied around
#' `initialAssignment` via `sbml_apply_species_initial_units()`), and
#' stoichiometries keyed by `speciesReference` id.
#'
#' @param parsed Bridge output from `sbml_run_bridge()`.
#' @return A named list of numeric symbol values.
sbml_build_symbols <- function(parsed) {
    symbols <- list()

    for (name in names(parsed$compartments)) {
        size <- parsed$compartments[[name]]$size
        if (!is.null(size)) {
            symbols[[name]] <- as.numeric(size)
        }
    }

    for (name in names(parsed$parameters)) {
        symbols[[name]] <- as.numeric(parsed$parameters[[name]])
    }

    for (item in parsed$species) {
        # Raw initialAmount / initialConcentration from the bridge.
        symbols[[item$id]] <- as.numeric(item$initial)
    }

    for (reaction in parsed$reactions) {
        for (ref in c(reaction$reactants, reaction$products)) {
            if (!is.null(ref$id) && nzchar(ref$id)) {
                symbols[[ref$id]] <- as.numeric(ref$stoich)
            }
        }
    }

    symbols
}

#' Convert species attribute initials into simulated units.
#'
#' Always converts from the bridge raw `initial` (amount or concentration as
#' declared), never from a previously converted symbol value. Call once with
#' the declared compartment size before resolving `initialAssignment` (so
#' formulas see concentrations), then again afterwards for species that were
#' not overwritten, using the resolved size.
#'
#' Species whose value was set by an `initialAssignment` are left unchanged
#' (already in the species units).
#'
#' Simulated units follow DNAr/SBML practice:
#' - `hasOnlySubstanceUnits` or 0D compartment: amount
#' - otherwise: concentration (= amount / volume)
#'
#' @param parsed Bridge output from `sbml_run_bridge()`.
#' @param symbols Symbol table (before or after `initialAssignment`).
#' @param assigned_species Species ids overwritten by initialAssignment.
#' @param strict If `TRUE`, missing/invalid volume is an error; if `FALSE`,
#'   skip conversion (needed when size exists only via `initialAssignment`).
#' @return Updated symbol table.
sbml_apply_species_initial_units <- function(parsed, symbols,
                                             assigned_species = character(),
                                             strict = TRUE) {
    assigned_species <- as.character(assigned_species)
    for (item in parsed$species) {
        id <- item$id
        if (id %in% assigned_species) {
            next
        }
        initial_type <- item$initialType
        if (is.null(initial_type) || !nzchar(initial_type)) {
            next
        }
        compartment <- parsed$compartments[[item$compartment]]
        dims <- if (is.null(compartment)) {
            NA_real_
        } else {
            as.numeric(compartment$spatialDimensions)
        }
        amount_state <- isTRUE(item$hasOnlySubstanceUnits) ||
            (!is.na(dims) && dims == 0)
        # Always from raw attribute so a second pass can redo with resolved V.
        value <- as.numeric(item$initial)
        if (length(value) != 1L || !is.finite(value)) {
            next
        }

        need_volume <- (identical(initial_type, "amount") && !amount_state) ||
            (identical(initial_type, "concentration") && amount_state)
        if (!need_volume) {
            symbols[[id]] <- value
            next
        }

        volume <- symbols[[item$compartment]]
        if (length(volume) != 1L) {
            volume <- if (is.null(compartment)) NULL else compartment$size
        }
        volume <- as.numeric(volume)
        # NULL -> numeric(0); length check avoids `!is.finite(numeric(0))`.
        volume_ok <- length(volume) == 1L && is.finite(volume) &&
            !(identical(initial_type, "amount") && !amount_state && volume == 0)
        if (!volume_ok) {
            if (isTRUE(strict)) {
                stop(
                    "volume invalido para especie '", id,
                    "' ao converter initial",
                    if (identical(initial_type, "amount")) "Amount" else "Concentration"
                )
            }
            next
        }

        if (identical(initial_type, "amount") && !amount_state) {
            symbols[[id]] <- value / volume
        } else {
            symbols[[id]] <- value * volume
        }
    }
    symbols
}

#' Per-species SBML metadata, including compartment volume for Papel C.
#'
#' @param parsed Bridge output from `sbml_run_bridge()`.
#' @return A named list of species info lists (`volume`, flags, ...).
sbml_build_species_info <- function(parsed) {
    info <- lapply(parsed$species, function(item) {
        list(
            compartment = item$compartment,
            volume = sbml_compartment_size(parsed, item$compartment),
            spatialDimensions = {
                compartment <- parsed$compartments[[item$compartment]]
                if (is.null(compartment)) NA_real_ else as.numeric(compartment$spatialDimensions)
            },
            boundaryCondition = isTRUE(item$boundaryCondition),
            constant = isTRUE(item$constant),
            hasOnlySubstanceUnits = isTRUE(item$hasOnlySubstanceUnits),
            conversionFactor = item$conversionFactor,
            initialType = item$initialType
        )
    })
    setNames(info, vapply(parsed$species, `[[`, character(1), "id"))
}

#' Species that must stay at their initial value (`boundaryCondition`/`constant`).
#'
#' @param parsed Bridge output from `sbml_run_bridge()`.
#' @return Character vector of fixed species ids.
sbml_fixed_species <- function(parsed) {
    if (length(parsed$species) == 0L) {
        return(character())
    }
    ids <- vapply(parsed$species, `[[`, character(1), "id")
    keep <- vapply(parsed$species, function(item) {
        isTRUE(item$boundaryCondition) || isTRUE(item$constant)
    }, logical(1))
    ids[keep]
}

#' Constant non-species symbols exposed as simulation outputs.
#'
#' Compartment sizes and parameters (including those filled by
#' `initialAssignment`) are returned; species and `speciesReference` ids
#' are excluded.
#'
#' @param parsed Bridge output from `sbml_run_bridge()`.
#' @param symbols Resolved symbol table.
#' @param species Character vector of species ids.
#' @return Named list of numeric constants.
sbml_build_outputs <- function(parsed, symbols, species) {
    ref_ids <- character()
    for (reaction in parsed$reactions) {
        for (ref in c(reaction$reactants, reaction$products)) {
            if (!is.null(ref$id) && nzchar(ref$id)) {
                ref_ids <- c(ref_ids, ref$id)
            }
        }
    }
    keep <- setdiff(names(symbols), c(species, ref_ids))
    if (length(keep) == 0L) {
        return(list())
    }
    lapply(symbols[keep], as.numeric)
}

#' Per-species derivative scales for SBML kinetic laws.
#'
#' SBML kinetic laws are kept as amount/time. Each species then converts to
#' its simulated unit via `react4(dx_scale=...)`:
#' - concentration in a 1D/2D/3D compartment: `conversionFactor / V`
#' - `hasOnlySubstanceUnits` or 0D compartment: `conversionFactor`
#'
#' Species-level `conversionFactor` overrides the model default; missing
#' factors are treated as 1. Compartment size prefers the resolved symbol
#' table (so `initialAssignment` to size is honored).
#'
#' @param parsed Bridge output from `sbml_run_bridge()`.
#' @param symbols Resolved symbol table.
#' @return Named numeric vector aligned with species ids.
sbml_species_conversion_factors <- function(parsed, symbols) {
    if (length(parsed$species) == 0L) {
        return(numeric())
    }
    model_cf <- parsed$conversionFactor
    setNames(vapply(parsed$species, function(item) {
        cf_name <- item$conversionFactor
        if (is.null(cf_name) || !nzchar(cf_name)) {
            cf_name <- model_cf
        }
        cf <- 1
        if (!is.null(cf_name) && nzchar(cf_name)) {
            value <- symbols[[cf_name]]
            if (is.null(value) || !is.finite(as.numeric(value))) {
                stop("conversionFactor '", cf_name, "' sem valor numerico")
            }
            cf <- as.numeric(value)
        }
        compartment <- parsed$compartments[[item$compartment]]
        dims <- if (is.null(compartment)) {
            NA_real_
        } else {
            as.numeric(compartment$spatialDimensions)
        }
        amount_state <- isTRUE(item$hasOnlySubstanceUnits) ||
            (!is.na(dims) && dims == 0)
        if (amount_state) {
            return(cf)
        }
        volume <- symbols[[item$compartment]]
        if (is.null(volume)) {
            volume <- if (is.null(compartment)) NULL else compartment$size
        }
        volume <- as.numeric(volume)
        if (!is.finite(volume) || volume <= 0) {
            stop(
                "volume invalido para especie '", item$id,
                "' ao montar dx_scale"
            )
        }
        cf / volume
    }, numeric(1)), vapply(parsed$species, `[[`, character(1), "id"))
}

#' Build `react4` forcing functions that hold fixed species at `ci`.
#'
#' With the forced value equal to the initial state, the
#' `(value - y) / tau` term in `react4` stays zero and the species does
#' not drift.
#'
#' @param model DNAr SBML model from `read_sbml_to_dnar()`.
#' @return Named list of constant functions, or `NULL` when none.
#' @export
sbml_forced_concentrations <- function(model) {
    fixed <- model$fixed_species
    if (is.null(fixed) || length(fixed) == 0L) {
        return(NULL)
    }
    lapply(model$ci[fixed], function(v) {
        force(v)
        function(t) v
    })
}

sbml_merge_values <- function(symbols, local_parameters) {
    values <- symbols
    for (name in names(local_parameters)) {
        values[[name]] <- as.numeric(local_parameters[[name]])
    }
    values
}

#' Resolve SBML `initialAssignment` formulas in dependency order (Kahn).
#'
#' Only dependencies among assignment targets are ordered; other symbols are
#' read from `values`. Cyclic target graphs raise an error.
#'
#' @param assignments List of lists with `symbol` and `formula`.
#' @param values Named list or vector of already-known symbol values.
#' @return Named list of values with assignment targets applied.
sbml_resolve_initial_assignments <- function(assignments, values) {
    if (length(assignments) == 0L) {
        return(values)
    }
    targets <- vapply(assignments, `[[`, character(1), "symbol")
    exprs <- setNames(
        lapply(assignments, function(a) parse(text = a$formula)[[1]]),
        targets
    )
    deps <- lapply(exprs, function(e) intersect(all.vars(e), targets))
    indeg <- vapply(deps, length, integer(1))
    queue <- targets[indeg == 0L]
    ord <- character()
    while (length(queue)) {
        n <- queue[[1]]
        queue <- queue[-1]
        ord <- c(ord, n)
        for (t in targets) {
            if (n %in% deps[[t]]) {
                indeg[[t]] <- indeg[[t]] - 1L
                if (indeg[[t]] == 0L) {
                    queue <- c(queue, t)
                }
            }
        }
    }
    if (length(ord) < length(targets)) {
        stop(
            "initialAssignment cíclico: ",
            paste(setdiff(targets, ord), collapse = ", ")
        )
    }
    env <- sbml_math_env(values)
    for (t in ord) {
        assign(t, eval(exprs[[t]], envir = env), envir = env)
    }
    as.list(unlist(mget(union(names(values), targets), envir = env)))
}

#' Read an SBML model into a DNAr CRN.
#'
#' Species attribute initials are converted into simulated units (concentration
#' unless `hasOnlySubstanceUnits` or 0D) twice around `initialAssignment`:
#' first with the declared compartment size (`strict = FALSE`) so IA formulas
#' read concentrations; then again from the raw attribute with the resolved
#' size for species that were not overwritten. With hosu +
#' `initialConcentration`, the value is multiplied by volume.
#' `initialAssignment` formulas are resolved (Kahn topological order) before
#' `ci` and `ki` are computed. Kinetic laws stay in amount/time; per-species
#' `conversion_factors` (`CF/V`, or `CF` for hosu/0D) are meant for
#' `react4(dx_scale=...)`. Constant compartments and parameters are listed
#' in `outputs` for non-species reference columns. Reactions with mass-action
#' kinetic laws expose their constant in `ki`. Non-mass-action laws are
#' reported in `kinetic_laws` as full-propensity functions. Species with
#' `boundaryCondition` or `constant` are listed in `fixed_species` for use
#' with `forced_concentrations`.
#'
#' @param xml_file_path Path to an SBML document.
#' @return A DNAr-compatible list containing `species`, `ci`, `reactions`,
#'   `ki`, `symbols`, `outputs`, `conversion_factors`, `species_info`,
#'   `fixed_species`, `rules` (metadata only; not applied), and SBML
#'   kinetic-law metadata.
#' @export
read_sbml_to_dnar <- function(xml_file_path) {
    parsed <- sbml_run_bridge(xml_file_path)
    if (length(parsed$errors) > 0L) {
        messages <- vapply(parsed$errors, `[[`, character(1), "message")
        stop("Invalid SBML: ", paste(messages, collapse = "; "))
    }
    if (length(parsed$compartments) > 1L) {
        stop("não suportado: múltiplos compartimentos")
    }

    assigned_symbols <- if (length(parsed$initialAssignments) == 0L) {
        character()
    } else {
        vapply(parsed$initialAssignments, `[[`, character(1), "symbol")
    }
    # Declared size first so IA math sees concentrations; skip if V is IA-only.
    symbols <- sbml_apply_species_initial_units(
        parsed,
        sbml_build_symbols(parsed),
        assigned_species = character(),
        strict = FALSE
    )
    symbols <- sbml_resolve_initial_assignments(
        parsed$initialAssignments,
        symbols
    )
    # Redo from raw with resolved V; keep IA-assigned species as-is.
    symbols <- sbml_apply_species_initial_units(
        parsed,
        symbols,
        assigned_species = assigned_symbols,
        strict = TRUE
    )
    species <- vapply(parsed$species, `[[`, character(1), "id")
    ci <- setNames(
        vapply(species, function(id) as.numeric(symbols[[id]]), numeric(1)),
        species
    )
    global_parameters <- parsed$parameters
    species_info <- sbml_build_species_info(parsed)
    fixed_species <- sbml_fixed_species(parsed)

    # DNAr parses a single net rate; SBML `reversible` with one kinetic law
    # is still that net flux. Emitting `<->` breaks the `->` splitter.
    reactions <- vapply(parsed$reactions, function(reaction) {
        paste(
            sbml_reaction_side(reaction$reactants),
            "->",
            sbml_reaction_side(reaction$products)
        )
    }, character(1))
    conversion_factors <- sbml_species_conversion_factors(parsed, symbols)
    outputs <- sbml_build_outputs(parsed, symbols, species)

    sbml_full_rate <- function(value_or_fun) {
        if (is.function(value_or_fun)) {
            attr(value_or_fun, "sbml_full_propensity") <- TRUE
            return(value_or_fun)
        }
        rate <- as.numeric(value_or_fun)
        fun <- function(t, y) rate
        attr(fun, "sbml_full_propensity") <- TRUE
        fun
    }

    # A kinetic-law formula is mass action when it consists of a constant
    # parameter multiplied by exactly the reactant powers. Otherwise ki is
    # a function that already includes the full propensity (no extra MA term).
    ki <- lapply(parsed$reactions, function(reaction) {
        formula <- reaction$formula
        values <- sbml_merge_values(symbols, reaction$parameters)
        if (is.null(formula) || !nzchar(formula)) {
            warning("Reaction '", reaction$id, "' has no kinetic law; using 1.")
            return(1)
        }
        constant_values <- values
        for (name in species) {
            constant_values[[name]] <- NULL
        }
        parameter_names <- names(constant_values)
        candidate <- formula
        for (name in parameter_names) {
            candidate <- gsub(
                paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                format(constant_values[[name]], scientific = FALSE, trim = TRUE),
                candidate,
                perl = TRUE
            )
        }
        candidate <- gsub("\\s+", "", candidate)
        candidate <- gsub("\\^1(?:\\.0+)?", "", candidate)
        residual <- formula
        for (name in c(parameter_names, species)) {
            residual <- gsub(
                paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                "",
                residual,
                perl = TRUE
            )
        }
        residual <- gsub("\\s+", "", residual)
        species_in_formula <- setNames(vapply(
            species,
            function(name) {
                grepl(
                    paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                    formula,
                    perl = TRUE
                )
            },
            logical(1)
        ), species)
        references_only_reactants <- all(vapply(
            species,
            function(name) {
                !isTRUE(species_in_formula[[name]]) || any(vapply(
                    reaction$reactants,
                    function(ref) identical(ref$species, name),
                    logical(1)
                ))
            },
            logical(1)
        ))
        if (!grepl("[A-Za-z_]", residual) &&
            !grepl("\\+", residual) &&
            references_only_reactants) {
            at_one <- candidate
            for (name in species) {
                at_one <- gsub(
                    paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                    "1",
                    at_one,
                    perl = TRUE
                )
            }
            value <- suppressWarnings(eval(parse(text = at_one), envir = baseenv()))
            if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
                if (!any(species_in_formula)) {
                    # Constant law (e.g. speciesReference / parameter only):
                    # stoichiometry still consumes/produces, but rate is bare.
                    return(sbml_full_rate(value))
                }
                return(as.numeric(value))
            }
        }
        expr <- tryCatch(parse(text = formula)[[1]], error = function(e) NULL)
        if (is.null(expr)) {
            warning("Could not reduce kinetic law for reaction '", reaction$id,
                    "'; using 1.")
            return(1)
        }
        constants <- constant_values
        fun <- function(t, y) {
            env <- sbml_math_env(c(constants, as.list(y)))
            val <- suppressWarnings(eval(expr, envir = env))
            if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
                stop(
                    "Non-finite kinetic law rate at t=", t,
                    " for species state: ",
                    paste(names(y), "=", as.numeric(y), collapse = ", ")
                )
            }
            as.numeric(val)
        }
        sbml_full_rate(fun)
    })

    # Arbitro por arvore sintatica (R/sbml_kinetics.R): uma lei so e acao de
    # massas se for constante * produto dos reagentes com expoente = estequiometria.
    # Se o detector diz "general" e a heuristica acima devolveu um numero, ela
    # era um falso positivo (ex.: lei `k*A` para a reacao A+B->C): usa a lei
    # completa como funcao de propensao em vez de k.
    kinetic_info <- sbml_classify_model(parsed, symbols, species)
    for (i in seq_along(ki)) {
        if (kinetic_info[[i]]$class == "general" && is.numeric(ki[[i]])) {
            reaction <- parsed$reactions[[i]]
            values <- sbml_merge_values(symbols, reaction$parameters)
            ki[[i]] <- sbml_general_rate_function(
                reaction$formula, values[setdiff(names(values), species)]
            )
        }
    }
    if (length(ki) > 0L &&
        all(vapply(ki, function(x) is.numeric(x) && length(x) == 1L, logical(1)))) {
        ki <- unlist(ki, use.names = FALSE)
    }

    structure(
        list(
            species = species,
            ci = ci,
            reactions = reactions,
            ki = ki,
            symbols = symbols,
            outputs = outputs,
            conversion_factors = conversion_factors,
            species_info = species_info,
            fixed_species = fixed_species,
            rules = if (is.null(parsed$rules)) list() else parsed$rules,
            kinetic_info = kinetic_info,
            kinetic_laws = lapply(parsed$reactions, function(reaction) {
                list(
                    id = reaction$id,
                    formula = reaction$formula,
                    parameters = sbml_merge_values(symbols, reaction$parameters)
                )
            })
        ),
        class = c("dnar_sbml", "list")
    )
}

#' @rdname read_sbml_to_dnar
#' @export
sbml_to_dnar <- read_sbml_to_dnar

convert__to_latex <- function(mathml_str) {
    if (!requireNamespace("V8", quietly = TRUE)) {
        stop("The 'V8' package is required. Run install.packages('V8') first.")
    }
    if (!requireNamespace("xml2", quietly = TRUE)) {
        stop("The 'xml2' package is required. Run install.packages('xml2') first.")
    }

    clean_mathml <- as.character(xml2::read_xml(mathml_str))
    clean_mathml <- gsub("\n", "", clean_mathml)
    clean_mathml <- gsub("'", "\\'", clean_mathml, fixed = TRUE)
    ctx <- V8::v8()
    js_url <- "https://unpkg.com"
    tryCatch({
        ctx$source(js_url)
    }, error = function(error) {
        stop("Failed to download the JavaScript converter library. Check your internet connection.")
    })
    ctx$assign("mml", clean_mathml)
    tryCatch({
        ctx$eval("MathMLToLaTeX.convert(mml)")
    }, error = function(error) {
        ctx$eval("typeof MathMLToLaTeX !== 'undefined' ? MathMLToLaTeX.convert(mml) : MathMLToLaTeX.default.convert(mml)")
    })
}

solve_latex_equation <- function(latex_str, ...) {
    if (!requireNamespace("latex2r", quietly = TRUE)) {
        stop("The 'latex2r' package is required.")
    }
    clean_latex <- gsub("^\\$|\\$$|^\\\\\\[|\\\\\\]$", "", trimws(latex_str))
    if (grepl("=", clean_latex)) {
        clean_latex <- unlist(strsplit(clean_latex, "="))[[2L]]
    }
    r_math_function <- latex2r::latex2fun(clean_latex)
    tryCatch(
        r_math_function(...),
        error = function(error) {
            stop(
                paste(
                    "Evaluation failed. Make sure you passed all necessary variable names correctly.",
                    "\nError:", error$message
                )
            )
        }
    )
}
