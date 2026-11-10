#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
read_arg <- function(name, default = NULL) {
    prefix <- paste0("--", name, "=")
    match <- args[startsWith(args, prefix)]
    if (length(match) == 0L) default else substring(match[1L], nchar(prefix) + 1L)
}

suite_dir <- read_arg(
    "suite-dir", file.path("external", "sbml-test-suite", "cases", "semantic")
)
output_dir <- read_arg("output-dir", file.path("results", "sbml"))
case_arg <- read_arg("cases")
case_ids <- if (is.null(case_arg)) NULL else trimws(strsplit(case_arg, ",")[[1L]])
absolute <- as.numeric(read_arg("absolute", "1e-6"))
relative <- as.numeric(read_arg("relative", "1e-4"))
preferred_level <- read_arg("preferred-level", "l3v2")
preflight <- !"--no-preflight" %in% args
if (length(absolute) != 1L || !is.finite(absolute) || absolute < 0 ||
    length(relative) != 1L || !is.finite(relative) || relative < 0) {
    stop("Tolerancias devem ser numeros finitos nao negativos.")
}

source(file.path("examples", "sbml", "testbench.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
manifest <- list(
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    git_head = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    worktree_dirty = length(system2("git", c("status", "--porcelain"), stdout = TRUE)) > 0L,
    suite_commit = system2("git", c("-C", shQuote(suite_dir), "rev-parse", "HEAD"), stdout = TRUE),
    libsbml = system2("python3", c("-c", shQuote("import libsbml; print(libsbml.getLibSBMLDottedVersion())")), stdout = TRUE),
    solver = "deSolve/lsoda", rtol = 1e-10, atol = 1e-12,
    absolute = absolute, relative = relative,
    preferred_level = preferred_level, preflight = preflight
)
jsonlite::write_json(manifest, file.path(output_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
suite <- run_sbml_cases(
    base_dir = suite_dir,
    case_ids = case_ids,
    absolute_tolerance = absolute,
    relative_tolerance = relative,
    preferred_level = preferred_level,
    preflight = preflight,
    progress = TRUE,
    checkpoint_dir = file.path(output_dir, "cases-rds")
)
paths <- write_sbml_baseline(suite, output_dir)

print(suite$status_counts, row.names = FALSE)
message("Artifacts:")
for (path in unlist(paths, use.names = FALSE)) message("  ", path)
