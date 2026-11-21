# Testes do papel D com casos reais da SBML Test Suite pelo circuito DNA.
# Precisam da suíte baixada em external/sbml-test-suite (scripts/download_sbml_test_suite.R).
library(testthat)

root <- normalizePath(if (file.exists(file.path("R", "sbml_dna.R"))) "." else file.path("..", ".."))
suite <- file.path(root, "external", "sbml-test-suite", "cases", "semantic")
skip_if_not(file.exists(file.path(root, "R", "sbml_kinetics.R")),
            "R/sbml_kinetics.R (papel C) ainda não está neste branch")
skip_if_not(dir.exists(suite), "SBML Test Suite não baixada em external/")
source(file.path(root, "examples", "sbml", "testbench_dna.R"), chdir = FALSE)
skip_if_not(sbml_available())

test_that("caso bimolecular reversível (00062) passa pelo DNA", {
    r <- run_sbml_case_dna("00062", suite)
    expect_equal(r$status_dna, "passou")
    expect_equal(r$status_crn, "passou")
    expect_lt(r$metrics$nrmse_dna, 1e-3)
    expect_gt(r$n_dna_reactions, 2L)
})

test_that("reverter alpha e beta devolve a mesma curva (00062)", {
    r <- run_sbml_case_dna("00062", suite, alpha = 2, beta = 10)
    expect_equal(r$status_dna, "passou")
})

test_that("espécie de contorno (00007) passa pelo DNA", {
    r <- run_sbml_case_dna("00007", suite)
    expect_equal(r$status_dna, "passou")
})

test_that("compartimento com volume diferente de 1 (00021) passa pelo DNA", {
    r <- run_sbml_case_dna("00021", suite)
    expect_equal(r$status_dna, "passou")
})

test_that("lei que não é ação de massas (00028) vira não suportado com motivo", {
    r <- run_sbml_case_dna("00028", suite)
    expect_equal(r$status_dna, "nao suportado")
    expect_match(r$reason, "^DNA: ")
})

test_that("falha de importação (01066) não interrompe o fluxo", {
    expect_no_error(r <- run_sbml_case_dna("01066", suite))
    expect_true(r$status_dna %in% c("passou", "falhou", "nao suportado"))
})
