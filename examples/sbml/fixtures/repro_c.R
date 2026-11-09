# Reproduz os achados do papel C sobre o suporte atual a SBML.
# Rode na raiz do repositorio, depois de copiar os c*.xml para a pasta indicada:
#   Rscript repro_c.R examples/sbml/fixtures
#
# ATENCAO: este script nao foi executado por quem o escreveu (sem R disponivel).
# Se der erro de sintaxe ou de pacote, corrija e anote o que mudou.

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1L) args[1L] else "examples/sbml/fixtures"

source(file.path("R", "parser.R"))
source(file.path("R", "util_functions.R"))
source(file.path("R", "crn_reactor.R"))
source(file.path("R", "sbml.R"))

esperado <- list(
    c0_controle = "reacao 'A + B -> C', ki = 0.5, aceito como mass-action",
    c1_volume = "reacao 'A + B -> C', ki = 0.5 (coeficiente 1 / V=2)",
    c2_parametro_local = "reacao 'A -> B', ki = 0.5 (o local vence o global)",
    c3_reversivel = "duas reacoes: 'A -> B' ki = 1 e 'B -> A' ki = 0.5",
    c4_ordem_errada = "marcado como nao suportado (lei k*A nao e mass-action de A + B)"
)

arquivos <- sort(list.files(dir, pattern = "^c[0-9].*\\.xml$", full.names = TRUE))
if (length(arquivos) == 0L) stop("Nenhum c*.xml encontrado em ", dir)

for (path in arquivos) {
    id <- sub("\\.xml$", "", basename(path))
    cat("\n==", id, "\n")
    cat("   esperado:", esperado[[id]], "\n")

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
        cat("   obtido: ERRO ->", conditionMessage(resultado), "\n")
    } else {
        cat("   obtido: reactions =", paste(resultado$reactions, collapse = " | "), "\n")
        cat("           ki        =", paste(resultado$ki, collapse = ", "), "\n")
    }
    for (aviso in avisos) cat("   aviso:", aviso, "\n")
}
