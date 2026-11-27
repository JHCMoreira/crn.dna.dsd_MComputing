#!/usr/bin/env Rscript

# Triagem das falhas de leitura do baseline A (papel B).
# Agrupa cases.csv por regex na mensagem e por tags SBML.
#
# Uso (na raiz do repositorio):
#   Rscript scripts/papel_b_triagem.R
#   Rscript scripts/papel_b_triagem.R --cases=results/sbml-acdd21a/cases.csv
#   Rscript scripts/papel_b_triagem.R --output=reports/papel-b/triagem.csv

args <- commandArgs(trailingOnly = TRUE)
read_arg <- function(name, default = NULL) {
    prefix <- paste0("--", name, "=")
    match <- args[startsWith(args, prefix)]
    if (length(match) == 0L) default else substring(match[1L], nchar(prefix) + 1L)
}

cases_path <- read_arg("cases", file.path("reports", "sbml-baseline", "cases.csv"))
output_path <- read_arg("output", file.path("reports", "papel-b", "triagem.csv"))

if (!file.exists(cases_path)) {
    stop("cases.csv nao encontrado: ", cases_path)
}

cases <- utils::read.csv(cases_path, stringsAsFactors = FALSE)
required <- c("reason", "test_tags", "status", "cause")
missing <- setdiff(required, names(cases))
if (length(missing) > 0L) {
    stop("Colunas ausentes em ", cases_path, ": ", paste(missing, collapse = ", "))
}

cases$reason[is.na(cases$reason)] <- ""
cases$test_tags[is.na(cases$test_tags)] <- ""
cases$status[is.na(cases$status)] <- ""
cases$cause[is.na(cases$cause)] <- ""

split_tags <- function(text) {
    if (!nzchar(text)) {
        return(character(0))
    }
    tags <- trimws(strsplit(text, ";", fixed = TRUE)[[1L]])
    tags[nzchar(tags)]
}

has_tags <- function(row_tags, required) {
    tags <- split_tags(row_tags)
    all(required %in% tags)
}

# Numeros do checkpoint do papel B: 224, 12, 3, 1, 3, 1 e 205.
# A linha de InitialAssignment no checkpoint (205) corresponde a
# cause==initial_assignment (197) + BoundaryCondition falhou sem mensagem (8).
groups <- data.frame(
    grupo = c(
        "b1_compartimento",
        "b2_contorno",
        "b5_species_reference_stoich",
        "b4_infinito",
        "b6_conversion_factor",
        "b5_id_na_formula",
        "b3_initial_assignment"
    ),
    criterio = c(
        "regex: object '(compartment|C)' not found",
        "tags: BoundaryCondition & HasOnlySubstanceUnits",
        "regex: _stoich",
        "regex: Infinity",
        "tag: ConversionFactors & status=falhou",
        "regex: object 'J0' not found",
        "cause: initial_assignment + BoundaryCondition falhou sem reason"
    ),
    stringsAsFactors = FALSE
)

hit <- list(
    grepl("object '(compartment|C)' not found", cases$reason, perl = TRUE),
    vapply(
        cases$test_tags, has_tags, logical(1),
        required = c("BoundaryCondition", "HasOnlySubstanceUnits")
    ),
    grepl("_stoich", cases$reason, fixed = TRUE),
    grepl("Infinity", cases$reason, fixed = TRUE),
    cases$status == "falhou" & vapply(
        cases$test_tags, has_tags, logical(1), required = "ConversionFactors"
    ),
    grepl("object 'J0' not found", cases$reason, fixed = TRUE),
    cases$cause == "initial_assignment" |
        (cases$status == "falhou" &
            !nzchar(cases$reason) &
            vapply(cases$test_tags, has_tags, logical(1), required = "BoundaryCondition"))
)

groups$n <- vapply(hit, sum, integer(1))
esperado <- c(224L, 12L, 3L, 1L, 3L, 1L, 205L)
if (!identical(as.integer(groups$n), esperado)) {
    warning(
        "Contagens diferem do checkpoint (224, 12, 3, 1, 3, 1, 205): ",
        paste(groups$n, collapse = ", ")
    )
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(groups, output_path, row.names = FALSE)

message("Triagem escrita em ", output_path)
print(groups, row.names = FALSE)
