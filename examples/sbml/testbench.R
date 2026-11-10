# Reproducible SBML Test Suite runner for DNAr.
#
# From the repository root:
#   source("examples/sbml/testbench.R")
#   suite <- run_sbml_cases("external/sbml-test-suite/cases/semantic",
#                           case_ids = c("00001", "00002"))
#
# The command-line entry point is scripts/run_sbml_baseline.R.

sbml_source_paths <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
}, character(1))
sbml_source_paths <- sbml_source_paths[!is.na(sbml_source_paths)]
sbml_testbench_file <- if (length(sbml_source_paths)) {
    tail(sbml_source_paths, 1L)
} else {
    file.path("examples", "sbml", "testbench.R")
}
sbml_project_root <- normalizePath(
    file.path(dirname(sbml_testbench_file), "..", ".."),
    mustWork = TRUE
)

sbml_previous_wd <- setwd(sbml_project_root)
tryCatch({
    source(file.path("R", "parser.R"), local = environment())
    source(file.path("R", "util_functions.R"), local = environment())
    source(file.path("R", "crn_reactor.R"), local = environment())
    source(file.path("R", "sbml.R"), local = environment())
}, finally = {
    setwd(sbml_previous_wd)
})

`%||%` <- function(value, fallback) {
    if (is.null(value) || length(value) == 0L) fallback else value
}

sbml_split_values <- function(value) {
    if (is.null(value) || !nzchar(trimws(value))) {
        return(character())
    }
    values <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
    values[nzchar(values)]
}

sbml_read_fields <- function(filepath) {
    if (!file.exists(filepath)) {
        stop("Arquivo nao encontrado: ", filepath)
    }
    lines <- readLines(filepath, warn = FALSE, encoding = "UTF-8")
    matches <- regexec(
        "^\\s*([A-Za-z][A-Za-z0-9]*):\\s*(.*?)\\s*$", lines,
        perl = TRUE
    )
    captures <- regmatches(lines, matches)
    captures <- captures[lengths(captures) == 3L]
    if (length(captures) == 0L) {
        return(list())
    }
    keys <- vapply(captures, `[[`, character(1), 2L)
    values <- vapply(captures, `[[`, character(1), 3L)
    as.list(stats::setNames(values, keys))
}

parse_sbml_settings <- function(filepath) {
    fields <- sbml_read_fields(filepath)
    read_number <- function(key, integer = FALSE, allow_blank = FALSE) {
        value <- fields[[key]]
        if (is.null(value) || !nzchar(value)) {
            if (allow_blank) return(if (integer) NA_integer_ else NA_real_)
            stop("Chave '", key, "' ausente em ", filepath)
        }
        parsed <- suppressWarnings(as.numeric(value))
        if (length(parsed) != 1L || !is.finite(parsed)) {
            stop("Valor invalido para '", key, "' em ", filepath)
        }
        if (integer) as.integer(parsed) else parsed
    }

    list(
        start = read_number("start", allow_blank = TRUE),
        duration = read_number("duration", allow_blank = TRUE),
        steps = read_number("steps", integer = TRUE, allow_blank = TRUE),
        variables = sbml_split_values(fields$variables),
        absolute = read_number("absolute"),
        relative = read_number("relative"),
        amount = sbml_split_values(fields$amount),
        concentration = sbml_split_values(fields$concentration)
    )
}

parse_sbml_metadata <- function(filepath) {
    fields <- sbml_read_fields(filepath)
    list(
        category = fields$category %||% "",
        synopsis = fields$synopsis %||% "",
        component_tags = sbml_split_values(fields$componentTags),
        test_tags = sbml_split_values(fields$testTags),
        test_type = fields$testType %||% "",
        levels = sbml_split_values(fields$levels),
        generated_by = fields$generatedBy %||% "",
        packages_present = sbml_split_values(fields$packagesPresent)
    )
}

discover_sbml_cases <- function(base_dir, case_ids = NULL) {
    if (!dir.exists(base_dir)) {
        stop("Diretorio da SBML Test Suite nao encontrado: ", base_dir)
    }
    settings <- list.files(
        base_dir,
        pattern = "-settings\\.txt$",
        recursive = TRUE,
        full.names = TRUE
    )
    if (length(settings) == 0L) {
        stop("Nenhum caso semantico encontrado em: ", base_dir)
    }
    ids <- sub("-settings\\.txt$", "", basename(settings))
    cases <- data.frame(
        case_id = ids,
        case_dir = dirname(settings),
        stringsAsFactors = FALSE
    )
    if (!is.null(case_ids)) {
        numeric_ids <- suppressWarnings(as.integer(case_ids))
        requested <- as.character(case_ids)
        requested[!is.na(numeric_ids)] <- sprintf("%05d", numeric_ids[!is.na(numeric_ids)])
        cases <- cases[cases$case_id %in% requested, , drop = FALSE]
        missing <- setdiff(requested, cases$case_id)
        if (length(missing) > 0L) {
            stop("Casos nao encontrados: ", paste(missing, collapse = ", "))
        }
    }
    cases[order(cases$case_id), , drop = FALSE]
}

find_sbml_case_files <- function(case_id, case_dir,
                                 preferred_level = "l3v2") {
    prefix <- paste0(case_id, "-sbml-")
    files <- list.files(
        case_dir,
        pattern = paste0("^", case_id, "-sbml-l[0-9]+v[0-9]+.*\\.xml$"),
        full.names = TRUE
    )
    files <- files[!grepl("-sedml\\.xml$", files)]
    if (length(files) == 0L) {
        fallback <- file.path(case_dir, paste0(case_id, ".xml"))
        files <- fallback[file.exists(fallback)]
    }
    if (length(files) == 0L) {
        stop("Arquivo SBML nao encontrado para o caso ", case_id)
    }

    variants <- sub("\\.xml$", "", substring(basename(files), nchar(prefix) + 1L))
    priorities <- unique(c(
        preferred_level,
        "l3v2", "l3v1", "l2v5", "l2v4", "l2v3", "l2v2", "l2v1", "l1v2"
    ))
    rank <- match(variants, priorities)
    missing_rank <- which(is.na(rank))
    rank[missing_rank] <- length(priorities) + seq_along(missing_rank)
    chosen <- which.min(rank)
    level_version <- regexec("^l([0-9]+)v([0-9]+)", variants[chosen])
    parsed <- regmatches(variants[chosen], level_version)[[1L]]
    list(
        path = files[chosen],
        variant = variants[chosen],
        level = if (length(parsed) == 3L) as.integer(parsed[2L]) else NA_integer_,
        version = if (length(parsed) == 3L) as.integer(parsed[3L]) else NA_integer_
    )
}

sbml_known_unsupported <- function(metadata) {
    if (!identical(metadata$test_type, "TimeCourse")) {
        return(paste0("testType nao suportado: ", metadata$test_type))
    }
    if (length(metadata$packages_present) > 0L) {
        return(paste0(
            "pacote SBML Level 3 nao suportado: ",
            paste(metadata$packages_present, collapse = ", ")
        ))
    }
    unsupported_components <- c(
        "AssignmentRule", "RateRule", "AlgebraicRule", "InitialAssignment",
        "EventNoDelay", "EventWithDelay", "EventPriority",
        "FunctionDefinition", "CSymbolDelay", "CSymbolRateOf"
    )
    found <- intersect(metadata$component_tags, unsupported_components)
    if (length(found) > 0L) {
        return(paste0(
            "componente ainda nao representado pelo DNAr: ",
            paste(found, collapse = ", ")
        ))
    }
    ""
}

compare_sbml_results <- function(result, reference, variables = NULL,
                                 absolute_tolerance = 1e-6,
                                 relative_tolerance = 1e-4) {
    result <- as.data.frame(result, check.names = FALSE)
    reference <- as.data.frame(reference, check.names = FALSE)
    if (!"time" %in% names(result) || !"time" %in% names(reference)) {
        stop("Resultado e referencia devem conter a coluna 'time'.")
    }
    if (is.null(variables)) {
        variables <- setdiff(names(reference), "time")
    }
    missing_result <- setdiff(variables, names(result))
    missing_reference <- setdiff(variables, names(reference))
    if (length(missing_result) > 0L || length(missing_reference) > 0L) {
        stop(
            "Variaveis ausentes. Resultado: ",
            paste(missing_result, collapse = ", "),
            "; referencia: ", paste(missing_reference, collapse = ", ")
        )
    }
    if (nrow(result) != nrow(reference)) {
        stop("Resultado e referencia possuem numero de linhas diferente.")
    }
    if (!length(variables) || !nrow(reference)) stop("Comparacao sem dados.")
    if (any(!is.finite(result$time)) || any(!is.finite(reference$time)) ||
        any(abs(result$time - reference$time) > 1e-10 * pmax(1, abs(reference$time)))) {
        stop("Tempos do resultado e da referencia diferem.")
    }

    values <- do.call(rbind, lapply(variables, function(name) {
        actual <- result[[name]]
        expected <- reference[[name]]
        both_nan <- is.nan(actual) & is.nan(expected)
        same_infinity <- is.infinite(actual) & is.infinite(expected) &
            sign(actual) == sign(expected)
        same_infinity[is.na(same_infinity)] <- FALSE
        finite <- is.finite(actual) & is.finite(expected)
        difference <- rep(Inf, length(expected))
        difference[finite] <- abs(actual[finite] - expected[finite])
        difference[both_nan | same_infinity] <- 0
        allowed <- absolute_tolerance + relative_tolerance * abs(expected)
        point_passed <- (finite & difference <= allowed) |
            both_nan | same_infinity
        relative <- difference / pmax(abs(expected), .Machine$double.eps)
        relative[!is.finite(relative) & difference == 0] <- 0
        data.frame(
            variable = name,
            max_absolute_error = max(difference),
            max_relative_error = max(relative),
            failing_points = sum(!point_passed),
            stringsAsFactors = FALSE
        )
    }))
    list(
        passed = all(values$failing_points == 0L),
        by_variable = values,
        max_absolute_error = max(values$max_absolute_error),
        max_relative_error = max(values$max_relative_error),
        failing_points = sum(values$failing_points),
        absolute_tolerance = absolute_tolerance,
        relative_tolerance = relative_tolerance
    )
}

sbml_failure_cause <- function(metadata, reason = "", status = "falhou") {
    if (identical(status, "passou")) return("nenhuma")
    components <- metadata$component_tags %||% character()
    tests <- metadata$test_tags %||% character()
    packages <- metadata$packages_present %||% character()
    if (length(packages) || any(grepl("^(comp|fbc):", components))) {
        return("pacote_level_3")
    }
    if (any(grepl("^Event", components))) return("eventos")
    if (any(c("AssignmentRule", "RateRule", "AlgebraicRule") %in% components)) {
        return("regras")
    }
    if ("InitialAssignment" %in% components) return("initial_assignment")
    if (any(c("FunctionDefinition", "CSymbolDelay", "CSymbolRateOf") %in% components) ||
        "UncommonMathML" %in% tests) return("mathml_ou_funcao")
    if (any(c("NonUnityCompartment", "MultiCompartment", "0D-Compartment",
              "VolumeConcentrationRates") %in% tests)) return("compartimentos_ou_unidades")
    if (grepl("compartment", reason, ignore.case = TRUE)) {
        return("compartimentos_ou_unidades")
    }
    if (any(c("BoundaryCondition", "ConstantSpecies", "HasOnlySubstanceUnits") %in% tests)) {
        return("semantica_de_especies")
    }
    if (any(c("NonUnityStoichiometry", "AssignedConstantStoichiometry",
              "AssignedVariableStoichiometry") %in% tests)) return("estequiometria")
    if ("ReversibleReaction" %in% tests) return("reacao_reversivel")
    if (grepl("variaveis.*ausentes|nao simulada", reason, ignore.case = TRUE)) {
        return("variavel_de_saida")
    }
    if (grepl("valid|libSBML|SBML", reason, ignore.case = TRUE)) return("importacao_ou_validacao")
    if (identical(status, "falhou") && !nzchar(reason)) return("divergencia_numerica")
    if (identical(status, "falhou")) return("simulacao_ou_comparacao")
    "nenhuma"
}

empty_sbml_metadata <- function() {
    list(
        component_tags = character(), test_tags = character(),
        test_type = "", packages_present = character()
    )
}

sbml_project_output <- function(simulation, settings, xml_file) {
    script <- file.path(sbml_project_root, "scripts", "sbml_output_units.py")
    output <- system2("python3", c(shQuote(script), shQuote(normalizePath(xml_file))),
                      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(output, "status")) && attr(output, "status") != 0L) {
        stop("Falha ao obter unidades de saida: ", paste(output, collapse = "\n"))
    }
    units <- jsonlite::fromJSON(paste(output, collapse = ""), simplifyVector = FALSE)
    result <- as.data.frame(simulation, check.names = FALSE)
    for (name in settings$variables) {
        info <- units[[name]]
        if (is.null(info)) next
        want_amount <- name %in% settings$amount
        want_concentration <- name %in% settings$concentration
        if (!want_amount && !want_concentration) {
            stop("Unidade de saida nao especificada para ", name)
        }
        if ((want_amount && !isTRUE(info$amount_state)) ||
            (want_concentration && isTRUE(info$amount_state))) {
            volume <- info$volume
            if (is.null(volume) || !is.finite(volume) || volume <= 0) {
                stop("Volume ausente ou invalido para converter ", name)
            }
            result[[name]] <- if (want_amount) result[[name]] * volume else result[[name]] / volume
        }
    }
    result
}

run_sbml_case <- function(case_id, base_dir,
                          absolute_tolerance = 1e-6,
                          relative_tolerance = 1e-4,
                          preferred_level = "l3v2",
                          engine = "desolve",
                          preflight = TRUE,
                          case_dir = NULL) {
    if (is.null(case_dir)) {
        discovered <- discover_sbml_cases(base_dir, case_ids = case_id)
        case_dir <- discovered$case_dir[1L]
    }
    metadata_file <- file.path(case_dir, paste0(case_id, "-model.m"))
    settings_file <- file.path(case_dir, paste0(case_id, "-settings.txt"))
    reference_file <- file.path(case_dir, paste0(case_id, "-results.csv"))
    metadata <- if (file.exists(metadata_file)) {
        parse_sbml_metadata(metadata_file)
    } else empty_sbml_metadata()
    settings <- parse_sbml_settings(settings_file)
    model_file <- find_sbml_case_files(case_id, case_dir, preferred_level)

    base_result <- list(
        case_id = case_id,
        case_dir = case_dir,
        xml_file = model_file$path,
        sbml_level = model_file$level,
        sbml_version = model_file$version,
        sbml_variant = model_file$variant,
        metadata = metadata,
        settings = settings,
        validation = NULL,
        comparison = NULL,
        warnings = character(),
        status = "falhou",
        cause = "simulacao_ou_comparacao",
        reason = ""
    )

    unsupported <- if (preflight) sbml_known_unsupported(metadata) else ""
    if (nzchar(unsupported)) {
        base_result$status <- "nao suportado"
        base_result$reason <- unsupported
        base_result$cause <- sbml_failure_cause(
            metadata, unsupported, base_result$status
        )
        return(base_result)
    }
    if (any(!is.finite(c(settings$start, settings$duration, settings$steps)))) {
        base_result$reason <- "start, duration ou steps ausente para TimeCourse"
        base_result$cause <- "infraestrutura_do_caso"
        return(base_result)
    }
    if (!file.exists(reference_file)) {
        base_result$status <- "nao suportado"
        base_result$reason <- "arquivo de referencia ausente"
        base_result$cause <- "infraestrutura_do_caso"
        return(base_result)
    }
    if (!sbml_available()) {
        base_result$reason <- "libSBML nao esta disponivel neste ambiente."
        base_result$cause <- "infraestrutura_do_runner"
        return(base_result)
    }

    warnings <- character()
    evaluated <- tryCatch(
        withCallingHandlers((function() {
            validation <- validate_sbml(model_file$path)
            base_result$validation <<- validation
            if (!isTRUE(validation$valid)) {
                messages <- vapply(
                    validation$errors, `[[`, character(1), "message"
                )
                stop("SBML invalido: ", paste(messages, collapse = "; "))
            }
            model <- read_sbml_to_dnar(model_file$path)
            missing_variables <- setdiff(settings$variables, model$species)
            if (length(missing_variables) > 0L) {
                return(list(
                    unsupported = paste0(
                        "variavel de referencia nao simulada pelo DNAr: ",
                        paste(missing_variables, collapse = ", ")
                    ),
                    validation = validation,
                    model = model
                ))
            }
            times <- seq(
                settings$start,
                settings$start + settings$duration,
                length.out = settings$steps + 1L
            )
            integration_times <- sort(unique(c(0, times)))
            simulation <- react(
                species = model$species,
                ci = model$ci,
                reactions = model$reactions,
                ki = model$ki,
                t = integration_times,
                engine = engine,
                rtol = 1e-10, atol = 1e-12
            )
            simulation <- as.data.frame(simulation, check.names = FALSE)
            simulation <- simulation[match(times, simulation$time), , drop = FALSE]
            simulation <- sbml_project_output(simulation, settings, model_file$path)
            reference <- read.csv(
                reference_file, check.names = FALSE,
                na.strings = c("NaN", "nan")
            )
            comparison <- compare_sbml_results(
                simulation,
                reference,
                variables = settings$variables,
                absolute_tolerance = absolute_tolerance %||% settings$absolute,
                relative_tolerance = relative_tolerance %||% settings$relative
            )
            list(
                validation = validation,
                model = model,
                simulation = simulation,
                reference = reference,
                comparison = comparison
            )
        })(), warning = function(warning) {
            warnings <<- c(warnings, conditionMessage(warning))
            invokeRestart("muffleWarning")
        }),
        error = function(error) error
    )

    base_result$warnings <- unique(warnings)
    if (inherits(evaluated, "error")) {
        base_result$reason <- conditionMessage(evaluated)
        base_result$cause <- sbml_failure_cause(
            metadata, base_result$reason, base_result$status
        )
        return(base_result)
    }
    base_result$validation <- evaluated$validation
    if (!is.null(evaluated$unsupported)) {
        base_result$status <- "nao suportado"
        base_result$reason <- evaluated$unsupported
        base_result$cause <- sbml_failure_cause(
            metadata, base_result$reason, base_result$status
        )
        return(base_result)
    }
    base_result$model <- evaluated$model
    base_result$simulation <- evaluated$simulation
    base_result$reference <- evaluated$reference
    base_result$comparison <- evaluated$comparison
    base_result$status <- if (evaluated$comparison$passed) "passou" else "falhou"
    base_result$cause <- sbml_failure_cause(
        metadata, "", base_result$status
    )
    base_result
}

sbml_result_row <- function(result) {
    comparison <- result$comparison
    validation <- result$validation
    suite_comparison <- if (!is.null(comparison)) {
        compare_sbml_results(result$simulation, result$reference,
                             variables = result$settings$variables,
                             absolute_tolerance = result$settings$absolute,
                             relative_tolerance = result$settings$relative)
    } else NULL
    data.frame(
        case_id = result$case_id,
        sbml_level = result$sbml_level,
        sbml_version = result$sbml_version,
        sbml_variant = result$sbml_variant,
        status = result$status,
        passed_suite_tolerances = if (is.null(suite_comparison)) NA else suite_comparison$passed,
        cause = result$cause,
        validation_valid = if (is.null(validation)) NA else validation$valid,
        max_absolute_error = if (is.null(comparison)) NA_real_ else comparison$max_absolute_error,
        max_relative_error = if (is.null(comparison)) NA_real_ else comparison$max_relative_error,
        failing_points = if (is.null(comparison)) NA_integer_ else comparison$failing_points,
        used_absolute_tolerance = if (is.null(comparison)) NA_real_ else comparison$absolute_tolerance,
        used_relative_tolerance = if (is.null(comparison)) NA_real_ else comparison$relative_tolerance,
        suite_absolute_tolerance = result$settings$absolute,
        suite_relative_tolerance = result$settings$relative,
        component_tags = paste(result$metadata$component_tags, collapse = ";"),
        test_tags = paste(result$metadata$test_tags, collapse = ";"),
        packages_present = paste(result$metadata$packages_present, collapse = ";"),
        warning_count = length(result$warnings),
        warnings = paste(result$warnings, collapse = " | "),
        reason = result$reason,
        xml_file = result$xml_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

summarize_sbml_tags <- function(results, field) {
    rows <- do.call(rbind, lapply(results, function(result) {
        tags <- result$metadata[[field]]
        if (length(tags) == 0L) return(NULL)
        data.frame(
            tag = tags,
            status = result$status,
            stringsAsFactors = FALSE
        )
    }))
    if (is.null(rows) || nrow(rows) == 0L) {
        return(data.frame(tag = character(), status = character(), count = integer()))
    }
    counts <- stats::aggregate(
        rep(1L, nrow(rows)), rows[c("tag", "status")], sum
    )
    names(counts)[3L] <- "count"
    counts[order(counts$tag, counts$status), , drop = FALSE]
}

run_sbml_cases <- function(base_dir, case_ids = NULL,
                           absolute_tolerance = 1e-6,
                           relative_tolerance = 1e-4,
                           preferred_level = "l3v2",
                           engine = "desolve",
                           preflight = TRUE,
                           progress = interactive(),
                           checkpoint_dir = NULL) {
    cases <- discover_sbml_cases(base_dir, case_ids)
    results <- vector("list", nrow(cases))
    for (index in seq_len(nrow(cases))) {
        case_id <- cases$case_id[index]
        if (progress) {
            message(sprintf("[%d/%d] %s", index, nrow(cases), case_id))
        }
        results[[index]] <- tryCatch(
            run_sbml_case(
                case_id = case_id,
                base_dir = base_dir,
                absolute_tolerance = absolute_tolerance,
                relative_tolerance = relative_tolerance,
                preferred_level = preferred_level,
                engine = engine,
                preflight = preflight,
                case_dir = cases$case_dir[index]
            ),
            error = function(error) {
                metadata <- empty_sbml_metadata()
                list(
                    case_id = case_id, sbml_level = NA_integer_,
                    sbml_version = NA_integer_, sbml_variant = "",
                    metadata = metadata,
                    settings = list(absolute = NA_real_, relative = NA_real_),
                    validation = NULL, comparison = NULL,
                    warnings = character(), status = "falhou",
                    cause = "infraestrutura_do_runner",
                    reason = conditionMessage(error), xml_file = ""
                )
            }
        )
        if (!is.null(checkpoint_dir)) {
            dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
            saveRDS(results[[index]], file.path(checkpoint_dir, paste0(case_id, ".rds")))
        }
    }
    summarize_sbml_cases(results)
}

summarize_sbml_cases <- function(results) {
    if (!length(results)) stop("Nenhum resultado para consolidar.")
    summary <- do.call(rbind, lapply(results, sbml_result_row))
    status_counts <- as.data.frame(table(summary$status), stringsAsFactors = FALSE)
    names(status_counts) <- c("status", "count")
    compatibility_groups <- summary[c("sbml_level", "sbml_version", "status")]
    for (field in c("sbml_level", "sbml_version")) {
        compatibility_groups[[field]] <- ifelse(is.na(compatibility_groups[[field]]),
                                                "unknown", as.character(compatibility_groups[[field]]))
    }
    compatibility <- stats::aggregate(
        rep(1L, nrow(summary)),
        compatibility_groups,
        sum
    )
    names(compatibility)[4L] <- "count"
    list(
        results = results,
        summary = summary,
        status_counts = status_counts,
        compatibility = compatibility,
        by_component_tag = summarize_sbml_tags(results, "component_tags"),
        by_test_tag = summarize_sbml_tags(results, "test_tags")
    )
}

write_sbml_baseline <- function(suite, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    causes <- stats::aggregate(rep(1L, nrow(suite$summary)),
                               suite$summary[c("cause", "status")], sum)
    names(causes)[3L] <- "count"
    artifacts <- list(
        cases = suite$summary,
        status = suite$status_counts,
        compatibility = suite$compatibility,
        causes = causes,
        component_tags = suite$by_component_tag,
        test_tags = suite$by_test_tag
    )
    paths <- file.path(output_dir, paste0(names(artifacts), ".csv"))
    Map(function(data, path) {
        utils::write.csv(data, path, row.names = FALSE, na = "")
        path
    }, artifacts, paths)
}

if (sys.nframe() == 0L) {
    message("Source this file or use scripts/run_sbml_baseline.R.")
}
