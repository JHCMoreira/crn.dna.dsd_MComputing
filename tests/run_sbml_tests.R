# Roda arquivos de teste do SBML sem instalar o pacote DNAr.
# Uso (na raiz do repositorio, dentro do container):
#   Rscript tests/run_sbml_tests.R                      # test_sbml_kinetics.R
#   Rscript tests/run_sbml_tests.R test_sbml.R test_sbml_testbench.R
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) args <- "test_sbml_kinetics.R"
suppressMessages({library(testthat); library(dplyr); library(ggplot2)})
for (f in c("parser.R", "util_functions.R", "crn_reactor.R",
            "sbml.R", "sbml_kinetics.R")) {
    source(file.path("R", f))
}
setwd(file.path("tests", "testthat"))
falhou <- FALSE
for (arq in args) {
    cat("\n=====", arq, "=====\n")
    linhas <- readLines(arq)
    linhas <- linhas[!grepl("^library\\(DNAr\\)", linhas)]
    tmp <- tempfile(pattern = "tmp_run_", tmpdir = ".", fileext = ".R")
    writeLines(linhas, tmp)
    res <- as.data.frame(test_file(tmp, reporter = "summary", env = globalenv()))
    unlink(tmp)
    cat(sprintf("%s -> testes: %d, falhas: %d, erros: %d, pulados: %d\n",
                arq, nrow(res), sum(res$failed), sum(res$error), sum(res$skipped)))
    if (sum(res$failed) + sum(res$error) > 0) falhou <- TRUE
}
quit(status = if (falhou) 1 else 0)
