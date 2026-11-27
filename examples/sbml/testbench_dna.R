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

#' Run one SBML Test Suite case through the CRN and the DNA circuit.
#'
#' @param case_id Case id, e.g. "00062".
#' @param base_dir Folder with the semantic cases.
#' @param qmax_factor,cmax_factor `qmax = qmax_factor * max(ki)` and
#'   `cmax = cmax_factor * beta * dna_material_bound()` (initial amounts plus
#'   what formation reactions can produce in the simulated time).
#' @param alpha,beta Time and concentration scales of the DNA circuit; the
#'   result is mapped back by simulating at `alpha * t` and dividing by `beta`.
#' @return A list with `status_dna` ("passou", "falhou", "nao suportado"),
#'   `reason`, `status_crn`, error metrics, parameters and the series
#'   `reference`, `crn` and `dna`.
run_sbml_case_dna <- function(case_id, base_dir,
                              qmax_factor = 1e4, cmax_factor = 1e4,
                              alpha = 1, beta = 1,
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
    unsupported <- function(reason) {
        result$status_dna <- "nao suportado"
        result$reason <- reason
        result
    }

    preflight <- sbml_known_unsupported(metadata)
    if (nzchar(preflight)) {
        return(unsupported(preflight))
    }
    if (any(!is.finite(c(settings$start, settings$duration, settings$steps)))) {
        result$reason <- "start, duration ou steps ausente para TimeCourse"
        return(result)
    }
    if (!file.exists(reference_file)) {
        return(unsupported("arquivo de referencia ausente"))
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
    if (inherits(model, "error")) {
        result$warnings <- unique(warnings)
        message <- conditionMessage(model)
        # O importador recusa de propósito o que está fora do escopo
        # (ex.: "não suportado: múltiplos compartimentos").
        if (grepl("^n(ã|a)o suportado", message)) {
            return(unsupported(paste("importacao:", message)))
        }
        result$reason <- paste("importacao:", message)
        return(result)
    }
    output_names <- names(model$outputs)
    missing_variables <- setdiff(settings$variables, c(model$species, output_names))
    if (length(missing_variables) > 0L) {
        return(unsupported(paste0("variavel de referencia nao simulada pelo DNAr: ",
                                  paste(missing_variables, collapse = ", "))))
    }

    prepared <- prepare_model_for_dna(model)
    compat <- dna_compatibility(prepared)
    if (!compat$ok) {
        return(unsupported(paste("DNA:", paste(compat$reasons, collapse = "; "))))
    }
    dna_crn <- dna_crn_from_model(prepared)

    times <- seq(settings$start, settings$start + settings$duration,
                 length.out = settings$steps + 1L)
    integration_times <- sort(unique(c(0, times)))
    reference <- read.csv(reference_file, check.names = FALSE,
                          na.strings = c("NaN", "nan"))
    time_col <- match("time", tolower(names(reference)))
    if (!is.na(time_col)) names(reference)[[time_col]] <- "time"

    # Monta a saída no formato da referência (unidades e variáveis extras).
    finish_output <- function(simulation) {
        simulation <- simulation[match(times, simulation$time), , drop = FALSE]
        simulation <- sbml_project_output(simulation, settings, model_file$path)
        for (name in intersect(output_names, settings$variables)) {
            simulation[[name]] <- model$outputs[[name]]
        }
        simulation
    }
    compare <- function(simulation) {
        compare_sbml_results(simulation, reference, variables = settings$variables,
                             absolute_tolerance = absolute_tolerance,
                             relative_tolerance = relative_tolerance)
    }

    simulated <- tryCatch(withCallingHandlers({
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

        # Circuito DNA.
        # cmax cobre o material inicial e o produzido por formações (0 -> X),
        # na escala de concentração do circuito (beta).
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

        # Só as espécies formais, com tempo e concentração revertidos.
        dna <- as.data.frame(dna_run$behavior, check.names = FALSE)
        dna <- dna[, c("time", dna_crn$species), drop = FALSE]
        dna$time <- integration_times
        dna[dna_crn$species] <- dna[dna_crn$species] / beta
        names(dna) <- c("time", unname(dna_crn$name_map[dna_crn$species]))

        list(crn = finish_output(crn), dna = finish_output(dna),
             qmax = qmax, cmax = cmax, elapsed = elapsed,
             n_species = length(dna_run$species),
             n_reactions = length(dna_run$reactions))
    }, warning = keep_warning), error = function(e) e)

    result$warnings <- unique(warnings)
    if (inherits(simulated, "error")) {
        result$reason <- paste("simulacao:", conditionMessage(simulated))
        return(result)
    }

    comparison_crn <- compare(simulated$crn)
    comparison_dna <- compare(simulated$dna)
    errors <- dna_error_metrics(simulated$dna, reference, settings$variables)

    result$status_crn <- if (comparison_crn$passed) "passou" else "falhou"
    result$status_dna <- if (comparison_dna$passed) "passou" else "falhou"
    result$comparison_dna <- comparison_dna
    result$metrics <- list(
        max_abs_err_dna = comparison_dna$max_absolute_error,
        max_rel_err_dna = comparison_dna$max_relative_error,
        mse_dna = errors$mse, nrmse_dna = errors$nrmse,
        max_rel_err_crn = comparison_crn$max_relative_error
    )
    result$params <- list(qmax = simulated$qmax, cmax = simulated$cmax,
                          alpha = alpha, beta = beta)
    result$n_dna_species <- simulated$n_species
    result$n_dna_reactions <- simulated$n_reactions
    result$sim_seconds <- simulated$elapsed
    result$reference <- reference
    result$crn <- simulated$crn
    result$dna <- simulated$dna
    result
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
