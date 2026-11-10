#!/usr/bin/env Rscript

# Download only the semantic cases. Sparse checkout avoids Windows path-length
# failures caused by unrelated Java dependencies in the full repository.

args <- commandArgs(trailingOnly = TRUE)
read_arg <- function(name, default) {
    prefix <- paste0("--", name, "=")
    match <- args[startsWith(args, prefix)]
    if (length(match) == 0L) default else substring(match[1L], nchar(prefix) + 1L)
}

target <- read_arg("target", file.path("external", "sbml-test-suite"))
version <- read_arg("version", "3.5.0")
repository <- "https://github.com/sbmlteam/sbml-test-suite.git"

run_git <- function(arguments, stdout = TRUE) {
    status <- system2("git", vapply(arguments, shQuote, character(1)), stdout = stdout, stderr = stdout)
    code <- if (is.character(status)) attr(status, "status") else status
    if (!is.null(code) && code != 0L) {
        stop("git failed: git ", paste(arguments, collapse = " "))
    }
    status
}

if (dir.exists(file.path(target, ".git"))) {
    current <- trimws(run_git(c("-C", target, "describe", "--tags", "--exact-match")))
    if (!identical(current, version)) {
        stop(
            "A suite already exists at ", target, " (", current,
            "); expected tag ", version, "."
        )
    }
    message("SBML Test Suite ", version, " already available at ", target)
    quit(status = 0L)
}
if (file.exists(target) || dir.exists(target)) {
    stop("Target exists but is not a git checkout: ", target)
}
dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)

run_git(c(
    "clone", "--filter=blob:none", "--no-checkout", "--depth", "1",
    "--branch", version, repository, target
), stdout = "")
run_git(c("-C", target, "sparse-checkout", "init", "--cone"), stdout = "")
run_git(c(
    "-C", target, "sparse-checkout", "set", "cases/semantic"
), stdout = "")
run_git(c("-C", target, "checkout", version), stdout = "")

case_count <- length(list.dirs(
    file.path(target, "cases", "semantic"), recursive = FALSE
))
message(
    "Downloaded SBML Test Suite ", version, " (", case_count,
    " semantic cases) to ", target
)
