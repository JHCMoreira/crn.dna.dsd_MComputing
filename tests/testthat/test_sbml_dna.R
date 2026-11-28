# Testes do papel D: compatibilidade SBML -> DNA e montagem da CRN traduzível.
# Carrega os arquivos de R/ com source(), como test_4domain_translation.R.
library(testthat)

root <- if (file.exists(file.path("R", "sbml_dna.R"))) "." else file.path("..", "..")
skip_if_not(file.exists(file.path(root, "R", "sbml_kinetics.R")),
            "R/sbml_kinetics.R (papel C) ainda não está neste branch")
old_wd <- setwd(root)
for (f in c("parser.R", "util_functions.R", "crn_reactor.R", "4domain_reactor.R",
            "sbml.R", "sbml_kinetics.R", "sbml_dna.R")) {
    source(file.path("R", f))
}
setwd(old_wd)

# Modelo "de brinquedo" no formato devolvido por read_sbml_to_dnar().
toy_model <- function(species, ci, reactions, info, conversion_factors = NULL,
                      fixed_species = character()) {
    info <- lapply(seq_along(info), function(i) {
        x <- info[[i]]
        defaults <- list(class = "mass_action", k = NA_real_, k_forward = NA_real_,
                         k_backward = NA_real_, reason = "", id = paste0("r", i),
                         reversible = FALSE, n_reactants = 1)
        modifyList(defaults, x)
    })
    if (is.null(conversion_factors)) {
        conversion_factors <- setNames(rep(1, length(species)), species)
    }
    list(
        species = species,
        ci = setNames(ci, species),
        reactions = reactions,
        ki = vapply(info, function(x) if (is.na(x$k)) 1 else x$k, numeric(1)),
        kinetic_info = info,
        kinetic_laws = lapply(info, function(x) list(id = x$id)),
        conversion_factors = conversion_factors,
        fixed_species = fixed_species
    )
}

reasons_of <- function(model) dna_compatibility(prepare_model_for_dna(model))$reasons

test_that("leitura e escrita de reações preservam a estequiometria", {
    parsed <- dna_parse_reaction("2A + B -> C")
    expect_equal(parsed$reactants$species, c("A", "B"))
    expect_equal(parsed$reactants$coef, c(2, 1))
    expect_equal(nrow(dna_parse_reaction("0 -> A")$reactants), 0L)
    expect_equal(dna_format_reaction(parsed), "2A + B -> C")
    expect_equal(dna_format_reaction(dna_parse_reaction("A -> 0")), "A -> 0")
})

test_that("modelo de ação de massas é compatível e usa o k do detector", {
    m <- toy_model(c("A", "B", "C"), c(1, 1, 0), "A + B -> C",
                   list(list(k = 0.5, n_reactants = 2)))
    p <- prepare_model_for_dna(m)
    expect_true(dna_compatibility(p)$ok)
    crn <- dna_crn_from_model(p)
    expect_equal(crn$reactions, "A + B -> C")
    expect_equal(crn$ki, 0.5)
})

test_that("reversível é separado em ida e volta antes da checagem", {
    m <- toy_model(c("A", "B"), c(1, 0), "A -> B",
                   list(list(class = "reversible_mass_action", k_forward = 1,
                             k_backward = 0.5, n_reactants = 1)))
    p <- prepare_model_for_dna(m)
    expect_true(dna_compatibility(p)$ok)
    crn <- dna_crn_from_model(p)
    expect_equal(crn$reactions, c("A -> B", "B -> A"))
    expect_equal(crn$ki, c(1, 0.5))
})

test_that("reversível com 3 produtos é recusado no passo inverso", {
    m <- toy_model(c("A", "B", "C", "D", "E"), c(1, 1, 0, 0, 0),
                   "A + B -> C + D + E",
                   list(list(class = "reversible_mass_action", k_forward = 1,
                             k_backward = 2, n_reactants = 2)))
    expect_match(reasons_of(m), "r1_rev", all = FALSE)
})

test_that("lei que não é ação de massas vira não suportado com o motivo do C", {
    m <- toy_model(c("A", "B"), c(1, 0), "A -> B",
                   list(list(class = "general", reason = "funcao 'exp' aplicada a uma especie")))
    r <- reasons_of(m)
    expect_length(r, 1L)
    expect_match(r, "reação r1: funcao 'exp'")
})

test_that("mais de 2 reagentes vira não suportado", {
    m <- toy_model(c("A", "B", "C"), c(1, 1, 0), "2A + B -> C",
                   list(list(k = 1, n_reactants = 3)))
    expect_match(reasons_of(m), "mais de 2 reagentes", all = FALSE)
})

test_that("reação de formação 0 -> X é aceita e entra no limite de material", {
    m <- toy_model(c("A", "B"), c(0, 3), c("0 -> A", "A -> B"),
                   list(list(k = 2, n_reactants = 0), list(k = 1)))
    p <- prepare_model_for_dna(m)
    expect_true(dna_compatibility(p)$ok)
    crn <- dna_crn_from_model(p)
    expect_equal(crn$reactions, c("0 -> A", "A -> B"))
    # material inicial (3) + formação 2 * duração 10
    expect_equal(dna_material_bound(crn, 10), 23)
    expect_equal(dna_material_bound(list(ci = 0, reactions = "A -> B", ki = 1), 10), 1)
})

test_that("estequiometria não inteira vira não suportado", {
    m <- toy_model(c("A", "B"), c(1, 0), "0.5A -> B",
                   list(list(k = 1, n_reactants = 0.5)))
    expect_match(reasons_of(m), "não inteira", all = FALSE)
})

test_that("dx_scale uniforme multiplica o ki; diferente entre espécies é recusado", {
    m <- toy_model(c("A", "B"), c(1, 0), "A -> B", list(list(k = 3)),
                   conversion_factors = c(A = 0.5, B = 0.5))
    expect_equal(dna_crn_from_model(prepare_model_for_dna(m))$ki, 1.5)

    m2 <- modifyList(m, list(conversion_factors = c(A = 0.5, B = 1)))
    expect_match(reasons_of(m2), "dx_scale", all = FALSE)

    # espécie fixa (contorno) não conta para a uniformidade
    m3 <- toy_model(c("A", "B"), c(1, 0), "A -> B", list(list(k = 3)),
                    conversion_factors = c(A = 1, B = 0.5), fixed_species = "A")
    expect_true(dna_compatibility(prepare_model_for_dna(m3))$ok)
})

test_that("espécie com nome reservado é renomeada e o mapa volta ao original", {
    m <- toy_model(c("T1", "B", "G"), c(1, 0, 1), c("2T1 -> B", "G + B -> T1"),
                   list(list(k = 1, n_reactants = 2), list(k = 2, n_reactants = 2)),
                   fixed_species = "G")
    crn <- dna_crn_from_model(prepare_model_for_dna(m))
    expect_equal(crn$species, c("sp_T1", "B", "sp_G"))
    expect_equal(crn$reactions, c("2sp_T1 -> B", "sp_G + B -> sp_T1"))
    expect_equal(crn$fixed_species, "sp_G")
    expect_equal(unname(crn$name_map[c("sp_T1", "B")]), c("T1", "B"))
})

test_that("modelos mínimos do papel C (c0-c4) são classificados para DNA", {
    skip_if_not(sbml_available())
    fixtures <- file.path(root, "examples", "sbml", "fixtures")
    skip_if_not(file.exists(file.path(fixtures, "c0_controle.xml")))
    check <- function(name) {
        dna_compatibility(prepare_model_for_dna(
            read_sbml_to_dnar(file.path(fixtures, name))))
    }
    expect_true(check("c0_controle.xml")$ok)
    expect_true(check("c1_volume.xml")$ok)
    expect_true(check("c2_parametro_local.xml")$ok)
    expect_true(check("c3_reversivel.xml")$ok)
    expect_false(check("c4_ordem_errada.xml")$ok)
})
