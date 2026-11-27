# Roda a SBML Test Suite pelo circuito de DNA (papel D) e grava o resumo.
#
# Uso, na raiz do repositório:
#   Rscript examples/sbml/run_dna_suite.R                     # todos os casos
#   Rscript examples/sbml/run_dna_suite.R --cases=00001,00062 # só alguns
#   Rscript examples/sbml/run_dna_suite.R --factor=1e5        # qmax e cmax
#
# Saídas em results/ (fora do Git):
#   dna_summary.csv  uma linha por caso (status CRN e DNA, erros, parâmetros)
#   dna_results.rds  resultados completos, com as séries para os gráficos

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
suite <- option("suite", file.path("external", "sbml-test-suite", "cases", "semantic"))
out_dir <- option("out", "results")
factor <- as.numeric(option("factor", "1e4"))
cases <- option("cases", "")
case_ids <- if (nzchar(cases)) strsplit(cases, ",", fixed = TRUE)[[1]] else NULL

source(file.path("examples", "sbml", "testbench_dna.R"))

started <- Sys.time()
out <- run_sbml_cases_dna(
    suite, case_ids,
    out_csv = file.path(out_dir, "dna_summary.csv"),
    qmax_factor = factor, cmax_factor = factor
)
saveRDS(out$results, file.path(out_dir, "dna_results.rds"))

s <- out$summary
cat("\nCasos:", nrow(s), "| fator qmax/cmax:", factor,
    "| tempo:", format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
cat("\nStatus no DNA:\n")
print(table(s$status_dna))
cat("\nCRN x DNA (casos simulados):\n")
print(table(CRN = s$status_crn, DNA = s$status_dna, useNA = "no"))
