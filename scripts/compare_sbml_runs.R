#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("Uso: Rscript scripts/compare_sbml_runs.R ANTES DEPOIS SAIDA")
}

read_run <- function(path) {
    list(
        manifest = jsonlite::read_json(file.path(path, "manifest.json")),
        cases = read.csv(file.path(path, "cases.csv"),
                         colClasses = c(case_id = "character", reason = "character"))
    )
}

before <- read_run(args[[1L]])
after <- read_run(args[[2L]])
settings <- c("suite_commit", "libsbml", "solver", "rtol", "atol",
              "absolute", "relative", "preferred_level", "preflight")
settings <- c(settings, intersect(c("r", "deSolve"),
                                  union(names(before$manifest), names(after$manifest))))
for (setting in settings) {
    a <- before$manifest[[setting]]
    b <- after$manifest[[setting]]
    if (is.null(a) || !identical(a, b)) stop("Configuracao diferente: ", setting)
}

a <- before$cases
b <- after$cases
if (!nrow(a) || anyNA(a$case_id) || anyNA(b$case_id) ||
    anyDuplicated(a$case_id) || anyDuplicated(b$case_id) ||
    !setequal(a$case_id, b$case_id)) {
    stop("As rodadas devem ter os mesmos IDs, sem duplicatas.")
}
b <- b[match(a$case_id, b$case_id), , drop = FALSE]
if (!identical(a$sbml_variant, b$sbml_variant)) {
    stop("As rodadas selecionaram variantes SBML diferentes.")
}

statuses <- c("passou", "falhou", "nao suportado")
if (!all(c(a$status, b$status) %in% statuses)) stop("Status de caso desconhecido.")
count <- function(status) as.integer(table(factor(status, levels = statuses)))
totals <- data.frame(status = statuses, antes = count(a$status), depois = count(b$status))
totals$delta <- totals$depois - totals$antes
transitions <- as.data.frame(table(
    antes = factor(a$status, levels = statuses),
    depois = factor(b$status, levels = statuses)
), responseName = "casos")

compared_before <- !is.na(a$passed_suite_tolerances)
compared_after <- !is.na(b$passed_suite_tolerances)
suite_changed <- compared_before & compared_after &
    a$passed_suite_tolerances != b$passed_suite_tolerances
changed <- which(a$status != b$status | compared_before != compared_after | suite_changed)
changes <- data.frame(
    case_id = a$case_id[changed],
    sbml_variant = a$sbml_variant[changed],
    antes = a$status[changed],
    depois = b$status[changed],
    comparado_antes = compared_before[changed],
    comparado_depois = compared_after[changed],
    passou_tol_suite_antes = a$passed_suite_tolerances[changed],
    passou_tol_suite_depois = b$passed_suite_tolerances[changed],
    erro_absoluto_antes = a$max_absolute_error[changed],
    erro_absoluto_depois = b$max_absolute_error[changed],
    motivo_antes = a$reason[changed],
    motivo_depois = b$reason[changed]
)

resource_counts <- function(status, tags, type) {
    tags[is.na(tags)] <- ""
    tags <- lapply(strsplit(tags, ";", fixed = TRUE), unique)
    resource <- unlist(tags, use.names = FALSE)
    counts <- table(resource, factor(rep(status, lengths(tags)), levels = statuses))
    result <- data.frame(
        tipo = rep(type, nrow(counts)),
        recurso = rownames(counts),
        casos = as.integer(rowSums(counts)),
        passou = as.integer(counts[, "passou"]),
        falhou = as.integer(counts[, "falhou"]),
        nao_suportado = as.integer(counts[, "nao suportado"])
    )
    result[nzchar(result$recurso), , drop = FALSE]
}

resource_matrix <- function(cases, stage) {
    group <- rbind(
        resource_counts(cases$status, cases$component_tags, "componentTags"),
        resource_counts(cases$status, cases$test_tags, "testTags")
    )
    data.frame(rodada = rep(stage, nrow(group)), group)
}

variants <- ifelse(is.na(b$sbml_variant) | !nzchar(b$sbml_variant), "desconhecida", b$sbml_variant)
versions <- as.data.frame(table(variante = rep(variants, 2L), rodada = rep(c("antes", "depois"), each = nrow(b)),
                                status = factor(c(a$status, b$status), levels = statuses)),
                          responseName = "casos")
outputs <- list(
    status = totals, transicoes = transitions, mudancas = changes,
    "matriz-recursos" = rbind(resource_matrix(a, "antes"), resource_matrix(b, "depois")),
    "matriz-nivel-versao" = versions
)
dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
for (name in names(outputs)) {
    write.csv(outputs[[name]], file.path(args[[3L]], paste0(name, ".csv")), row.names = FALSE, na = "")
}
print(totals, row.names = FALSE)
print(transitions[transitions$casos > 0L, ], row.names = FALSE)
