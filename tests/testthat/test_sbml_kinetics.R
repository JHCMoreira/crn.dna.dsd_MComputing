library(testthat)
library(DNAr)

# Testes do detector de acao de massas (papel C).
# Os de "sbml_classify_kinetic_law" sao puros (sem libSBML); os que leem
# as fixtures c0-c4 precisam do libSBML.

refs <- function(...) {
    lapply(list(...), function(x) list(species = x[[1]], stoich = x[[2]]))
}
classify <- function(formula, reactants, products = list(),
                     constants = list(k = 2, k1 = 1, k2 = 0.5, V = 2)) {
    sbml_classify_kinetic_law(formula, reactants, products, constants,
                              species = c("A", "B", "C"))
}

test_that("acao de massas: ordem dos fatores, parenteses e potencias nao importam", {
    ab <- refs(list("A", 1), list("B", 1))
    expect_equal(classify("k * A * B", ab)$k, 2)
    expect_equal(classify("B * (A * k)", ab)$k, 2)
    expect_equal(classify("V * k * A * B", ab)$k, 4)  # volume entra em k
    expect_equal(classify("k * A^2", refs(list("A", 2)))$k, 2)
    expect_equal(classify("k * A * A", refs(list("A", 2)))$class, "mass_action")
    expect_equal(classify("k * A / V", refs(list("A", 1)))$k, 1)
    expect_equal(classify("k", list())$class, "mass_action")  # fonte de ordem 0
})

test_that("lei que nao bate com a estequiometria vira 'general' com razao", {
    ab <- refs(list("A", 1), list("B", 1))
    r <- classify("k * A", ab)
    expect_equal(r$class, "general")
    expect_match(r$reason, "falta o reagente B")
    r <- classify("k * A^2 * B", ab)
    expect_match(r$reason, "expoente de A")
    r <- classify("k * A * B * C", ab)
    expect_match(r$reason, "C aparece na lei mas nao e reagente")
    expect_match(classify("k * A / (k2 + A)", refs(list("A", 1)))$reason, "soma")
    expect_match(classify("k * factorial(A)", refs(list("A", 1)))$reason, "factorial")
    expect_match(classify("k * A^B", refs(list("A", 1), list("B", 1)))$reason, "expoente")
    expect_equal(classify("k", refs(list("A", 1)))$class, "general")  # ordem 0 num reagente
})

test_that("reversivel f - g e separado em direto e inverso", {
    r <- classify("k1 * A - k2 * B", refs(list("A", 1)), refs(list("B", 1)))
    expect_equal(r$class, "reversible_mass_action")
    expect_equal(c(r$k_forward, r$k_backward), c(1, 0.5))
    r <- classify("-k2 * B + k1 * A", refs(list("A", 1)), refs(list("B", 1)))
    expect_equal(c(r$k_forward, r$k_backward), c(1, 0.5))
    r <- classify("k1 * A - k2 * A", refs(list("A", 1)), refs(list("B", 1)))
    expect_equal(r$class, "general")
    expect_match(r$reason, "termo inverso")
})

fixture <- function(name) {
    file.path("..", "..", "examples", "sbml", "fixtures", name)
}

test_that("fixtures c0-c4: classificacao e ki", {
    skip_if_not(sbml_available())
    load <- function(n) suppressWarnings(read_sbml_to_dnar(fixture(n)))

    c0 <- load("c0_controle.xml")
    expect_equal(c0$kinetic_info[[1]]$class, "mass_action")
    expect_equal(c0$ki, 0.5)

    c1 <- load("c1_volume.xml")  # k * V, por desenho (dx_scale divide por V)
    expect_equal(c1$kinetic_info[[1]]$class, "mass_action")
    expect_equal(c1$ki, 1)

    c2 <- load("c2_parametro_local.xml")  # local vence o global
    expect_equal(c2$ki, 0.5)

    c4 <- load("c4_ordem_errada.xml")  # lei k*A para A + B -> C
    expect_equal(c4$kinetic_info[[1]]$class, "general")
    expect_match(c4$kinetic_info[[1]]$reason, "falta o reagente B")
    expect_true(is.function(c4$ki[[1]]))
    expect_true(isTRUE(attr(c4$ki[[1]], "sbml_full_propensity")))
    expect_equal(c4$ki[[1]](0, c(A = 3, B = 7, C = 0)), 6)  # k*A, sem B
})

test_that("reversivel c3: fluxo liquido por padrao, dois passos com o split", {
    skip_if_not(sbml_available())
    c3 <- suppressWarnings(read_sbml_to_dnar(fixture("c3_reversivel.xml")))
    expect_equal(c3$reactions, "A -> B")
    expect_equal(c3$kinetic_info[[1]]$class, "reversible_mass_action")

    s <- sbml_split_reversible(c3)
    expect_equal(s$reactions, c("A -> B", "B -> A"))
    expect_equal(s$ki, c(1, 0.5))
    expect_equal(s$split_from, c(1L, 1L))
    expect_true(all(sbml_kinetic_report(s)$dna_ok))

    t <- seq(0, 5, by = 0.5)
    net <- react4(c3$species, c3$ci, c3$reactions, c3$ki, t)
    two <- react4(s$species, s$ci, s$reactions, s$ki, t)
    expect_equal(net$A, two$A, tolerance = 1e-6)
    expect_equal(net$B, two$B, tolerance = 1e-6)
})

test_that("tabela de suporte: crn_ok / dna_ok / razao", {
    skip_if_not(sbml_available())
    rep <- sbml_kinetic_report(
        suppressWarnings(read_sbml_to_dnar(fixture("c4_ordem_errada.xml")))
    )
    expect_equal(rep$kinetic_class, "general")
    expect_true(rep$crn_ok)
    expect_false(rep$dna_ok)
    expect_match(rep$reason, "B")
})
