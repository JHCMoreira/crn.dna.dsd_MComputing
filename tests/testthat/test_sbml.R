library(testthat)
library(DNAr)

test_that("libSBML validates and imports the example case", {
    skip_if_not(sbml_available())

    path <- file.path("..", "..", "examples", "sbml", "cases", "00028.xml")
    validation <- validate_sbml(path)
    expect_true(validation$valid)
    expect_length(validation$errors, 0)

    crn <- read_sbml_to_dnar(path)
    expect_equal(crn$species, c("S1", "S2"))
    expect_equal(unname(crn$ci), c(1, 0))
    expect_equal(crn$reactions, "S1 -> S2")
    expect_equal(crn$ki, 0.96)
    expect_equal(crn$kinetic_laws[[1]]$formula,
                 "factorial(ceil(p1 * S1)) * p2^-1")
})

test_that("missing SBML files fail explicitly", {
    skip_if_not(sbml_available())
    expect_error(read_sbml_to_dnar("does-not-exist.xml"), "does not exist")
})

test_that("bridge exports extended SBML metadata", {
    skip_if_not(sbml_available())

    path <- file.path("..", "..", "examples", "sbml", "fixtures",
                      "b6_conversion_factor.xml")
    parsed <- sbml_run_bridge(path)

    expect_equal(parsed$conversionFactor, "modelconv")
    expect_equal(parsed$compartments$compartment$spatialDimensions, 3)
    expect_equal(parsed$compartments$compartment$size, 1)
    expect_equal(parsed$species[[1]]$compartment, "compartment")
    expect_false(parsed$species[[1]]$boundaryCondition)
    expect_false(parsed$species[[1]]$constant)
    expect_false(parsed$species[[1]]$hasOnlySubstanceUnits)
    expect_equal(parsed$species[[1]]$conversionFactor, "s1conv")
    expect_equal(parsed$species[[1]]$initialType, "amount")
    expect_null(parsed$species[[2]]$conversionFactor)
    expect_null(parsed$reactions[[1]]$reactants[[1]]$id)
})

test_that("symbols table unlocks compartment kinetics and exposes volume", {
    skip_if_not(sbml_available())

    path <- file.path("..", "..", "examples", "sbml", "fixtures",
                      "b1_compartimento.xml")
    crn <- suppressWarnings(read_sbml_to_dnar(path))
    expect_equal(crn$species, c("S1", "S2"))
    expect_equal(crn$symbols$compartment, 2)
    expect_equal(crn$symbols$k, 1)
    expect_equal(crn$symbols$S1, 1)
    expect_equal(crn$species_info$S1$volume, 2)
    expect_equal(crn$species_info$S2$volume, 2)
    # Amount-rate law compartment*k*S1 still reduces; Papel C divides by V.
    expect_equal(crn$ki, 2)

    stoich <- file.path("..", "..", "examples", "sbml", "fixtures",
                        "b5_species_reference.xml")
    crn_stoich <- suppressWarnings(read_sbml_to_dnar(stoich))
    expect_equal(crn_stoich$symbols$S1_stoich, 1)
    expect_equal(crn_stoich$ki, 0.5)

    env <- sbml_math_env(list(x = 1.5, flag = TRUE))
    expect_equal(eval(parse(text = "ceil(x)"), envir = env), 2)
    expect_equal(eval(parse(text = "piecewise(10, flag, 0)"), envir = env), 10)
})

test_that("bridge preserves Inf/NaN and speciesReference ids", {
    skip_if_not(sbml_available())

    infinite <- file.path("..", "..", "examples", "sbml", "fixtures",
                          "b4_infinito.xml")
    parsed_inf <- sbml_run_bridge(infinite)
    expect_true(is.infinite(parsed_inf$parameters$k) && parsed_inf$parameters$k > 0)

    case_00951 <- file.path(
        "..", "..", "external", "sbml-test-suite", "cases", "semantic",
        "00951", "00951-sbml-l3v2.xml"
    )
    skip_if_not(file.exists(case_00951))
    parsed_case <- sbml_run_bridge(case_00951)
    expect_true(is.infinite(parsed_case$parameters$P) && parsed_case$parameters$P > 0)
    expect_true(is.infinite(parsed_case$parameters$Q) && parsed_case$parameters$Q < 0)
    expect_true(is.nan(parsed_case$parameters$R))
    crn <- read_sbml_to_dnar(case_00951)
    expect_equal(crn$species, character())

    stoich <- file.path("..", "..", "examples", "sbml", "fixtures",
                        "b5_species_reference.xml")
    parsed_ref <- sbml_run_bridge(stoich)
    expect_equal(parsed_ref$reactions[[1]]$reactants[[1]]$id, "S1_stoich")

    ia <- file.path("..", "..", "examples", "sbml", "fixtures",
                    "b3_initial_assignment.xml")
    parsed_ia <- sbml_run_bridge(ia)
    expect_equal(length(parsed_ia$initialAssignments), 3L)
    expect_equal(parsed_ia$initialAssignments[[1]]$symbol, "p1")
    expect_equal(parsed_ia$initialAssignments[[1]]$formula, "2")
    expect_length(parsed_ia$rules, 0L)
})
