# Testes de regressão da tradução 4-domain e do react3().
# Carrega os arquivos de R/ diretamente (o pacote DNAr não instala hoje
# por causa de R/SIGNALS_N_SYSTEMS.R).
library(testthat)

# Alguns arquivos de R/ fazem source() com caminho relativo à raiz do projeto.
root <- if (file.exists(file.path("R", "4domain_reactor.R"))) "." else file.path("..", "..")
old_wd <- setwd(root)
for (f in c("parser.R", "util_functions.R", "crn_reactor.R", "4domain_reactor.R")) {
    source(file.path("R", f))
}
setwd(old_wd)

max_rel_error <- function(dna, crn, species) {
    max(vapply(species, function(s) {
        max(abs(dna[[s]] - crn[[s]])) / max(abs(crn[[s]]))
    }, numeric(1)))
}

test_that("translate_4domain_crn gera reações com a seta '->' aceita pelo parser", {
    crn <- translate_4domain_crn(
        species = c("A", "B", "C"), ci = c(1e3, 1e3, 0),
        reactions = c("A + B -> C", "C -> A"), ki = c(1e-7, 1e-4),
        qmax = 1e-3, cmax = 1e5, alpha = 1, beta = 1, t = seq(0, 10, 1)
    )
    expect_false(any(grepl("-->", crn$reactions, fixed = TRUE)))
    expect_true(all(grepl(" -> ", crn$reactions, fixed = TRUE)))
    # todas as espécies auxiliares precisam aparecer nas reações traduzidas
    sto <- get_M(check_crn(crn$species, crn$ci, crn$reactions, crn$ki, 0), crn$species)
    expect_true(all(colSums(abs(sto$react) + abs(sto$prod)) > 0))
})

test_that("react_4domain reproduz A -> B (unimolecular)", {
    sp <- c("A", "B"); ci <- c(1e4, 0); rx <- "A -> B"; k <- 5e-5
    t <- seq(0, 20000, 100)
    crn <- as.data.frame(react4(sp, ci, rx, k, t))
    dna <- as.data.frame(react_4domain(sp, ci, rx, k, qmax = k * 1e4,
                                       cmax = 1e8, alpha = 1, beta = 1,
                                       t = t)$behavior)
    expect_lt(max_rel_error(dna, crn, sp), 1e-2)
})

test_that("react_4domain reproduz A + B -> C (bimolecular)", {
    sp <- c("A", "B", "C"); ci <- c(1e3, 1e3, 0); rx <- "A + B -> C"; k <- 1e-7
    t <- seq(0, 72000, 100)
    crn <- as.data.frame(react4(sp, ci, rx, k, t))
    dna <- as.data.frame(react_4domain(sp, ci, rx, k, qmax = 1e-3, cmax = 1e5,
                                       alpha = 1, beta = 1, t = t)$behavior)
    expect_lt(max_rel_error(dna, crn, sp), 1e-2)
})

test_that("react3 aceita argumentos extras do deSolve (atol, rtol)", {
    expect_no_error(
        react3(c("A", "B"), c(1, 0), "A -> B", 1, seq(0, 1, 0.1),
               atol = 1e-12, rtol = 1e-10)
    )
})
