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
        "compartment_sizes = {}",
        "for i in range(model.getNumCompartments()):",
        "    compartment = model.getCompartment(i)",
        "    size = compartment.getSize() if compartment.isSetSize() else None",
        "    compartment_sizes[compartment.getId()] = size",
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
        "    initial = None",
        "    initial_type = None",
        "    if item.isSetInitialConcentration():",
        "        initial = item.getInitialConcentration()",
        "        initial_type = 'concentration'",
        "    elif item.isSetInitialAmount():",
        "        initial = item.getInitialAmount()",
        "        initial_type = 'amount'",
        "        size = compartment_sizes.get(item.getCompartment())",
        "        if not item.getHasOnlySubstanceUnits() and size not in (None, 0):",
        "            initial = initial / size",
        "    if initial is None:",
        "        initial = 0.0",
        "    species.append({",
        "        'id': item.getId(),",
        "        'initial': number(initial),",
        "        'initialType': initial_type,",
        "        'compartment': item.getCompartment(),",
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

#' Read an SBML model into a DNAr CRN.
#'
#' Species initial amounts are converted to concentrations using their
#' compartment size. Reactions with mass-action kinetic laws expose their
#' constant in `ki`. Non-mass-action laws are reported in `kinetic_laws` and
#' use their value at the initial state as a compatibility fallback.
#'
#' @param xml_file_path Path to an SBML document.
#' @return A DNAr-compatible list containing `species`, `ci`, `reactions`,
#'   `ki`, and SBML kinetic-law metadata.
#' @export
read_sbml_to_dnar <- function(xml_file_path) {
    parsed <- sbml_run_bridge(xml_file_path)
    if (length(parsed$errors) > 0L) {
        messages <- vapply(parsed$errors, `[[`, character(1), "message")
        stop("Invalid SBML: ", paste(messages, collapse = "; "))
    }

    species <- vapply(parsed$species, `[[`, character(1), "id")
    ci <- setNames(
        vapply(parsed$species, function(item) as.numeric(item$initial), numeric(1)),
        species
    )
    global_parameters <- parsed$parameters

    reactions <- vapply(parsed$reactions, function(reaction) {
        paste(
            sbml_reaction_side(reaction$reactants),
            if (isTRUE(reaction$reversible)) "<->" else "->",
            sbml_reaction_side(reaction$products)
        )
    }, character(1))

    # A kinetic-law formula is mass action when it consists of a constant
    # parameter multiplied by exactly the reactant powers.
    ki <- vapply(parsed$reactions, function(reaction) {
        formula <- reaction$formula
        values <- c(global_parameters, reaction$parameters)
        if (is.null(formula) || !nzchar(formula)) {
            warning("Reaction '", reaction$id, "' has no kinetic law; using 1.")
            return(1)
        }
        parameter_names <- names(values)
        reactant_term <- paste(
            vapply(reaction$reactants, function(ref) {
                paste0(ref$species, "^", ref$stoich)
            }, character(1)),
            collapse = " * "
        )
        candidate <- formula
        for (name in parameter_names) {
            candidate <- gsub(
                paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                format(values[[name]], scientific = FALSE, trim = TRUE),
                candidate,
                perl = TRUE
            )
        }
        candidate <- gsub("\\s+", "", candidate)
        candidate <- gsub("\\^1(?:\\.0+)?", "", candidate)
        species_names <- vapply(parsed$species, `[[`, character(1), "id")
        residual <- formula
        for (name in c(parameter_names, species_names)) {
            residual <- gsub(
                paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                "",
                residual,
                perl = TRUE
            )
        }
        residual <- gsub("\\s+", "", residual)
        references_only_reactants <- all(vapply(
            species_names,
            function(name) {
                !grepl(
                    paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                    formula,
                    perl = TRUE
                ) || any(vapply(
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
            for (name in species_names) {
                at_one <- gsub(
                    paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                    "1",
                    at_one,
                    perl = TRUE
                )
            }
            value <- suppressWarnings(eval(parse(text = at_one), envir = baseenv()))
            if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
                return(as.numeric(value))
            }
        }
        initial_values <- c(as.list(ci), as.list(values))
        for (name in names(initial_values)) {
            candidate <- gsub(
                paste0("(?<![A-Za-z0-9_])", name, "(?![A-Za-z0-9_])"),
                format(initial_values[[name]], scientific = FALSE, trim = TRUE),
                candidate,
                perl = TRUE
            )
        }
        candidate <- gsub("\\bceil\\(", "ceiling(", candidate)
        value <- suppressWarnings(eval(parse(text = candidate), envir = baseenv()))
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
            warning("Could not reduce kinetic law for reaction '", reaction$id,
                    "'; using 1.")
            return(1)
        }
        warning("Reaction '", reaction$id,
                "' is not mass action; ki uses its initial-state rate.")
        as.numeric(value)
    }, numeric(1))

    structure(
        list(
            species = species,
            ci = ci,
            reactions = reactions,
            ki = ki,
            kinetic_laws = lapply(parsed$reactions, function(reaction) {
                list(id = reaction$id, formula = reaction$formula,
                     parameters = c(global_parameters, reaction$parameters))
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
