# Reproduz os achados do papel B sobre falhas de leitura SBML.
# Rode na raiz do repositorio:
#   Rscript examples/sbml/fixtures/repro_b.R
#   Rscript examples/sbml/fixtures/repro_b.R examples/sbml/fixtures
#
# A saida e gravada em examples/sbml/fixtures/baseline_acdd21a_b.txt

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1L) args[1L] else "examples/sbml/fixtures"
output <- file.path(dir, "baseline_acdd21a_b.txt")

source(file.path("R", "sbml.R"))
source(file.path("R", "sbml_kinetics.R"))

esperado <- list(
    b1_compartimento = "erro object 'compartment' not found",
    b2_contorno = "S1 e S2 importadas; S1 de contorno tratada como dinamica",
    b3_initial_assignment = "IAs ignoradas; S1 inicial 0 (cadeia p1->p2->S1 perdida)",
    b4_infinito = "erro de parse JSON por Infinity no bridge",
    b5_species_reference = "erro object 'S1_stoich' not found",
    b6_conversion_factor = "conversionFactor ignorado; S1 ci sem fator"
)

arquivos <- sort(list.files(dir, pattern = "^b[0-9].*\\.xml$", full.names = TRUE))
if (length(arquivos) == 0L) {
    stop("Nenhum b*.xml encontrado em ", dir)
}

git_head <- Sys.getenv("GIT_HEAD", unset = "")
if (!nzchar(git_head)) {
    git_head <- tryCatch(
        {
            out <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
            status <- attr(out, "status")
            if (!is.null(status) && status != 0L) {
                "unavailable"
            } else {
                paste(out, collapse = "")
            }
        },
        error = function(e) "unavailable"
    )
}
if (!nzchar(git_head)) {
    git_head <- "unavailable"
}

lines <- c(
    paste0("Baseline papel B (leitura SBML), gerado em ",
           format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("git HEAD: ", git_head),
    ""
)

for (path in arquivos) {
    id <- sub("\\.xml$", "", basename(path))
    block <- c(paste0("== ", id, " =="), paste0("   esperado: ", esperado[[id]]))

    avisos <- character(0)
    resultado <- tryCatch(
        withCallingHandlers(
            read_sbml_to_dnar(path),
            warning = function(w) {
                avisos <<- c(avisos, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) e
    )

    if (inherits(resultado, "error")) {
        block <- c(block, paste0("   obtido: ERRO -> ", conditionMessage(resultado)))
    } else {
        block <- c(
            block,
            paste0(
                "   obtido: species = ", paste(resultado$species, collapse = ", "),
                "; ci = ", paste(paste(names(resultado$ci), resultado$ci, sep = "="), collapse = ", ")
            ),
            paste0(
                "           reactions = ", paste(resultado$reactions, collapse = " | "),
                "; ki = ", paste(resultado$ki, collapse = ", ")
            )
        )
    }
    for (aviso in avisos) {
        block <- c(block, paste0("   aviso: ", aviso))
    }
    block <- c(block, "")
    lines <- c(lines, block)
    cat(paste(block, collapse = "\n"), "\n", sep = "")
}

writeLines(lines, output)
message("Baseline gravado em ", output)
