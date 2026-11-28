# Papel D: roda casos da SBML Test Suite pelo circuito de DNA (4-domain).
#
# Na raiz do repositório:
#   source("examples/sbml/testbench_dna.R")
#   r <- run_sbml_case_dna("00062", "external/sbml-test-suite/cases/semantic")
#   r$status_dna; r$metrics
#
# Reaproveita o testbench do grupo (examples/sbml/testbench.R): descoberta dos
# casos, settings, filtro de recursos, conversão de unidades de saída e
# comparação com a referência.

dna_testbench_dir <- local({
    paths <- vapply(sys.frames(), function(frame) {
        if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
    }, character(1))
    paths <- paths[!is.na(paths)]
    if (length(paths)) dirname(tail(paths, 1L)) else file.path("examples", "sbml")
})
source(file.path(dna_testbench_dir, "testbench.R"))
local({
    previous <- setwd(sbml_project_root)
    on.exit(setwd(previous))
    source(file.path("R", "4domain_reactor.R"), local = globalenv())
    source(file.path("R", "sbml_dna.R"), local = globalenv())
})

# Erros do circuito DNA em relação à referência, além do passou/falhou.
dna_error_metrics <- function(simulation, reference, variables) {
    diffs <- unlist(lapply(variables, function(v) simulation[[v]] - reference[[v]]))
    scale <- max(abs(unlist(reference[variables])), na.rm = TRUE)
    mse <- mean(diffs^2, na.rm = TRUE)
    list(
        mse = mse,
        nrmse = if (is.finite(scale) && scale > 0) sqrt(mse) / scale else NA_real_
    )
}

#' Load one SBML Test Suite case and prepare it for the DNA circuit.
#'
#' Does the expensive part once (reading the SBML with libSBML, checking DNA
#' compatibility, simulating the CRN and finding the output units), so that
#' [simulate_case_dna()] can then try many `qmax`/`cmax` values cheaply.
#'
#' @param case_id Case id, e.g. "00062".
#' @param base_dir Folder with the semantic cases.
#' @return A list with `result` (the partial result) and `context` (`NULL`
#'   when the case stops before the DNA simulation).
load_sbml_case_dna <- function(case_id, base_dir,
                               absolute_tolerance = 1e-6,
                               relative_tolerance = 1e-4,
                               preferred_level = "l3v2",
                               case_dir = NULL) {
    if (is.null(case_dir)) {
        case_dir <- discover_sbml_cases(base_dir, case_ids = case_id)$case_dir[1L]
    }
    metadata_file <- file.path(case_dir, paste0(case_id, "-model.m"))
    metadata <- if (file.exists(metadata_file)) {
        parse_sbml_metadata(metadata_file)
    } else empty_sbml_metadata()
    settings <- parse_sbml_settings(file.path(case_dir, paste0(case_id, "-settings.txt")))
    reference_file <- file.path(case_dir, paste0(case_id, "-results.csv"))
    model_file <- find_sbml_case_files(case_id, case_dir, preferred_level)

    result <- list(
        case_id = case_id, xml_file = model_file$path,
        level_version = model_file$variant,
        status_dna = "falhou", status_crn = NA_character_, reason = "",
        metrics = NULL, params = NULL,
        n_dna_species = NA_integer_, n_dna_reactions = NA_integer_,
        sim_seconds = NA_real_, warnings = character()
    )
    stop_here <- function(reason, status = "falhou") {
        result$status_dna <- status
        result$reason <- reason
        list(result = result, context = NULL)
    }

    preflight <- sbml_known_unsupported(metadata)
    if (nzchar(preflight)) {
        return(stop_here(preflight, "nao suportado"))
    }
    if (any(!is.finite(c(settings$start, settings$duration, settings$steps)))) {
        return(stop_here("start, duration ou steps ausente para TimeCourse"))
    }
    if (!file.exists(reference_file)) {
        return(stop_here("arquivo de referencia ausente", "nao suportado"))
    }

    warnings <- character()
    keep_warning <- function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }

    model <- tryCatch(
        withCallingHandlers(read_sbml_to_dnar(model_file$path), warning = keep_warning),
        error = function(e) e
    )
    result$warnings <- unique(warnings)
    if (inherits(model, "error")) {
        message <- conditionMessage(model)
        # O importador recusa de propósito o que está fora do escopo
        # (ex.: "não suportado: múltiplos compartimentos").
        if (grepl("^n(ã|a)o suportado", message)) {
            return(stop_here(paste("importacao:", message), "nao suportado"))
        }
        return(stop_here(paste("importacao:", message)))
    }
    output_names <- names(model$outputs)
    missing_variables <- setdiff(settings$variables, c(model$species, output_names))
    if (length(missing_variables) > 0L) {
        return(stop_here(paste0("variavel de referencia nao simulada pelo DNAr: ",
                                paste(missing_variables, collapse = ", ")),
                         "nao suportado"))
    }

    prepared <- prepare_model_for_dna(model)
    compat <- dna_compatibility(prepared)
    if (!compat$ok) {
        return(stop_here(paste("DNA:", paste(compat$reasons, collapse = "; ")),
                         "nao suportado"))
    }
    dna_crn <- dna_crn_from_model(prepared)

    times <- seq(settings$start, settings$start + settings$duration,
                 length.out = settings$steps + 1L)
    integration_times <- sort(unique(c(0, times)))
    reference <- read.csv(reference_file, check.names = FALSE,
                          na.strings = c("NaN", "nan"))
    time_col <- match("time", tolower(names(reference)))
    if (!is.na(time_col)) names(reference)[[time_col]] <- "time"

    prepared_run <- tryCatch(withCallingHandlers({
        # CRN de referência: mesmo caminho do testbench do grupo.
        forced_crn <- NULL
        if (length(model$fixed_species) > 0L) {
            forced_crn <- lapply(model$ci[model$fixed_species], function(v) {
                force(v); function(t) v
            })
        }
        crn <- as.data.frame(react4(
            species = model$species, ci = model$ci, reactions = model$reactions,
            ki = model$ki, t = integration_times,
            forced_concentrations = forced_crn, dx_scale = model$conversion_factors,
            rtol = 1e-10, atol = 1e-12
        ), check.names = FALSE)
        crn <- crn[match(times, crn$time), , drop = FALSE]

        # Conversão de unidades de saída (quantidade x concentração): é um
        # fator constante por variável, obtido uma vez com o testbench do grupo.
        ones <- as.data.frame(as.list(setNames(rep(1, length(model$species) + 1L),
                                               c("time", model$species))),
                              check.names = FALSE)
        # model = model: usa o volume resolvido por initialAssignment (como o testbench).
        unit_factor <- unlist(sbml_project_output(ones, settings, model_file$path,
                                                  model = model)[model$species])

        list(crn = crn, unit_factor = unit_factor)
    }, warning = keep_warning), error = function(e) e)
    result$warnings <- unique(warnings)
    if (inherits(prepared_run, "error")) {
        return(stop_here(paste("simulacao:", conditionMessage(prepared_run))))
    }

    context <- list(
        settings = settings, model = model, dna_crn = dna_crn,
        times = times, integration_times = integration_times,
        reference = reference, unit_factor = prepared_run$unit_factor,
        absolute_tolerance = absolute_tolerance,
        relative_tolerance = relative_tolerance
    )
    context$crn <- dna_finish_output(prepared_run$crn, context)
    context$comparison_crn <- dna_compare(context$crn, context)
    list(result = result, context = context)
}

# Saída no formato da referência: unidades e variáveis que não são espécies.
dna_finish_output <- function(simulation, context) {
    for (name in names(context$unit_factor)) {
        simulation[[name]] <- simulation[[name]] * context$unit_factor[[name]]
    }
    outputs <- context$model$outputs
    for (name in intersect(names(outputs), context$settings$variables)) {
        simulation[[name]] <- outputs[[name]]
    }
    simulation
}

dna_compare <- function(simulation, context) {
    compare_sbml_results(simulation, context$reference,
                         variables = context$settings$variables,
                         absolute_tolerance = context$absolute_tolerance,
                         relative_tolerance = context$relative_tolerance)
}

#' Simulate the DNA circuit of a loaded case and compare it with the reference.
#'
#' Default factors come from the calibration (examples/sbml/dna_calibration.R):
#' qmax and cmax at 1e5 gave the most approved cases. With `stiff_retry`, a
#' stiff run (the solver does not finish) is repeated with qmax ten times
#' smaller, down to `min_qmax_factor`.
#'
#' @param loaded Output of [load_sbml_case_dna()].
#' @param qmax_factor,cmax_factor `qmax = qmax_factor * max(ki)` and
#'   `cmax = cmax_factor * beta * dna_material_bound()` (initial amounts plus
#'   what formation reactions can produce in the simulated time).
#' @param alpha,beta Time and concentration scales of the DNA circuit; the
#'   result is mapped back by simulating at `alpha * t` and dividing by `beta`.
#' @param stiff_retry Retry stiff runs with a smaller qmax.
#' @return The result list (see [run_sbml_case_dna()]); `params$retries`
#'   counts the retries and `params$qmax_factor` is the factor actually used.
simulate_case_dna <- function(loaded, qmax_factor = 1e5, cmax_factor = 1e5,
                              alpha = 1, beta = 1,
                              stiff_retry = FALSE, min_qmax_factor = 1e2) {
    result <- simulate_case_dna_once(loaded, qmax_factor, cmax_factor, alpha, beta)
    retries <- 0L
    while (stiff_retry && grepl("sistema rigido", result$reason) &&
           qmax_factor / 10 >= min_qmax_factor) {
        qmax_factor <- qmax_factor / 10
        retries <- retries + 1L
        result <- simulate_case_dna_once(loaded, qmax_factor, cmax_factor, alpha, beta)
    }
    if (!is.null(result$params)) {
        result$params$retries <- retries
    }
    result
}

simulate_case_dna_once <- function(loaded, qmax_factor, cmax_factor, alpha, beta) {
    result <- loaded$result
    context <- loaded$context
    if (is.null(context)) {
        return(result)
    }
    dna_crn <- context$dna_crn
    integration_times <- context$integration_times

    warnings <- character()
    simulated <- tryCatch(withCallingHandlers({
        qmax <- qmax_factor * max(dna_crn$ki)
        cmax <- cmax_factor * beta *
            dna_material_bound(dna_crn, max(integration_times) - min(integration_times))
        forced_dna <- NULL
        if (length(dna_crn$fixed_species) > 0L) {
            fixed_ci <- dna_crn$ci[match(dna_crn$fixed_species, dna_crn$species)]
            forced_dna <- setNames(lapply(fixed_ci, function(v) {
                v <- v * beta; force(v); function(t) v
            }), dna_crn$fixed_species)
        }
        started <- Sys.time()
        dna_run <- react_4domain(
            species = dna_crn$species, ci = dna_crn$ci,
            reactions = dna_crn$reactions, ki = dna_crn$ki,
            qmax = qmax, cmax = cmax, alpha = alpha, beta = beta,
            t = integration_times * alpha,
            forced_concentrations = forced_dna,
            rtol = 1e-10, atol = 1e-12
        )
        elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

        dna <- as.data.frame(dna_run$behavior, check.names = FALSE)
        if (nrow(dna) < length(integration_times)) {
            stop("integracao interrompida: sistema rigido ",
                 "(qmax x cmax alto demais para o resolvedor)")
        }
        # Só as espécies formais, com tempo e concentração revertidos.
        dna <- dna[, c("time", dna_crn$species), drop = FALSE]
        dna$time <- integration_times
        dna[dna_crn$species] <- dna[dna_crn$species] / beta
        names(dna) <- c("time", unname(dna_crn$name_map[dna_crn$species]))
        dna <- dna[match(context$times, dna$time), , drop = FALSE]

        list(dna = dna_finish_output(dna, context),
             qmax = qmax, cmax = cmax, elapsed = elapsed,
             n_species = length(dna_run$species),
             n_reactions = length(dna_run$reactions))
    }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) e)

    result$warnings <- unique(c(result$warnings, warnings))
    result$status_crn <- if (context$comparison_crn$passed) "passou" else "falhou"
    result$params <- list(qmax_factor = qmax_factor, cmax_factor = cmax_factor,
                          alpha = alpha, beta = beta)
    if (inherits(simulated, "error")) {
        result$reason <- paste("simulacao:", conditionMessage(simulated))
        return(result)
    }

    comparison_dna <- dna_compare(simulated$dna, context)
    errors <- dna_error_metrics(simulated$dna, context$reference,
                                context$settings$variables)
    result$status_dna <- if (comparison_dna$passed) "passou" else "falhou"
    result$comparison_dna <- comparison_dna
    result$metrics <- list(
        max_abs_err_dna = comparison_dna$max_absolute_error,
        max_rel_err_dna = comparison_dna$max_relative_error,
        mse_dna = errors$mse, nrmse_dna = errors$nrmse,
        max_rel_err_crn = context$comparison_crn$max_relative_error
    )
    result$params$qmax <- simulated$qmax
    result$params$cmax <- simulated$cmax
    result$n_dna_species <- simulated$n_species
    result$n_dna_reactions <- simulated$n_reactions
    result$sim_seconds <- simulated$elapsed
    result$reference <- context$reference
    result$crn <- context$crn
    result$dna <- simulated$dna
    result
}

#' Run one SBML Test Suite case through the CRN and the DNA circuit.
#'
#' Same as [load_sbml_case_dna()] followed by [simulate_case_dna()], with the
#' calibrated defaults (qmax and cmax factors 1e5, retry when stiff).
#'
#' @return A list with `status_dna` ("passou", "falhou", "nao suportado"),
#'   `reason`, `status_crn`, error metrics, parameters and the series
#'   `reference`, `crn` and `dna`.
run_sbml_case_dna <- function(case_id, base_dir,
                              qmax_factor = 1e5, cmax_factor = 1e5,
                              alpha = 1, beta = 1, stiff_retry = TRUE,
                              absolute_tolerance = 1e-6,
                              relative_tolerance = 1e-4,
                              preferred_level = "l3v2",
                              case_dir = NULL) {
    loaded <- load_sbml_case_dna(case_id, base_dir,
                                 absolute_tolerance = absolute_tolerance,
                                 relative_tolerance = relative_tolerance,
                                 preferred_level = preferred_level,
                                 case_dir = case_dir)
    simulate_case_dna(loaded, qmax_factor = qmax_factor, cmax_factor = cmax_factor,
                      alpha = alpha, beta = beta, stiff_retry = stiff_retry)
}

# Uma linha da tabela de resumo para um resultado de run_sbml_case_dna().
dna_case_row <- function(result) {
    pick <- function(x) if (is.null(x) || length(x) == 0L) NA else x
    metric <- function(name) pick(result$metrics[[name]])
    param <- function(name) pick(result$params[[name]])
    data.frame(
        case_id = result$case_id,
        level_version = pick(result$level_version),
        status_crn = pick(result$status_crn),
        status_dna = result$status_dna,
        max_rel_err_dna = metric("max_rel_err_dna"),
        max_abs_err_dna = metric("max_abs_err_dna"),
        mse_dna = metric("mse_dna"),
        nrmse_dna = metric("nrmse_dna"),
        max_rel_err_crn = metric("max_rel_err_crn"),
        qmax_factor = param("qmax_factor"),
        cmax_factor = param("cmax_factor"),
        retries = param("retries"),
        qmax = param("qmax"),
        cmax = param("cmax"),
        alpha = param("alpha"),
        beta = param("beta"),
        n_dna_species = pick(result$n_dna_species),
        n_dna_reactions = pick(result$n_dna_reactions),
        sim_seconds = pick(result$sim_seconds),
        reason = pick(result$reason),
        stringsAsFactors = FALSE
    )
}

#' Run several SBML Test Suite cases through the DNA circuit.
#'
#' Each case runs inside its own `tryCatch`, so an unexpected error is
#' recorded as "falhou" with the message and the batch goes on.
#'
#' @param base_dir Folder with the semantic cases.
#' @param case_ids Case ids to run (all cases when `NULL`).
#' @param out_csv Optional path of the summary CSV.
#' @param verbose Print one line per case.
#' @param ... Passed to [run_sbml_case_dna()] (e.g. `qmax_factor`).
#' @return A list with `results` (one list per case) and `summary`
#'   (data.frame, one row per case).
run_sbml_cases_dna <- function(base_dir, case_ids = NULL, out_csv = NULL,
                               verbose = TRUE, ...) {
    cases <- discover_sbml_cases(base_dir, case_ids)
    results <- vector("list", nrow(cases))
    for (i in seq_len(nrow(cases))) {
        id <- cases$case_id[[i]]
        results[[i]] <- tryCatch(
            run_sbml_case_dna(id, base_dir, case_dir = cases$case_dir[[i]], ...),
            error = function(e) list(case_id = id, status_dna = "falhou",
                                     reason = paste("erro inesperado:", conditionMessage(e)))
        )
        if (verbose) {
            cat(sprintf("[%d/%d] %s DNA: %s\n", i, nrow(cases), id,
                        results[[i]]$status_dna))
        }
    }
    summary <- do.call(rbind, lapply(results, dna_case_row))
    if (!is.null(out_csv)) {
        dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
        write.csv(summary, out_csv, row.names = FALSE)
    }
    list(results = results, summary = summary)
}
