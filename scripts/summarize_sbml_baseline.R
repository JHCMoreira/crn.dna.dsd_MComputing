#!/usr/bin/env Rscript
# Rebuild CSV summaries from saved per-case results without resimulation.
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1L] else file.path("results", "sbml")
source(file.path("examples", "sbml", "testbench.R"))
files <- list.files(file.path(output_dir, "cases-rds"), "\\.rds$", full.names = TRUE)
results <- lapply(files, readRDS)
for (index in seq_along(results)) {
    item <- results[[index]]
    item$cause <- sbml_failure_cause(item$metadata, item$reason, item$status)
    results[[index]] <- item
}
suite <- summarize_sbml_cases(results)
write_sbml_baseline(suite, output_dir)
print(suite$status_counts, row.names = FALSE)
