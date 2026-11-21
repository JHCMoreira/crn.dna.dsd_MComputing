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
    expect_true(is.function(crn$ki[[1]]))
    expect_true(isTRUE(attr(crn$ki[[1]], "sbml_full_propensity")))
    expect_equal(crn$ki[[1]](0, c(S1 = 1, S2 = 0)), 0.96)
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

test_that("boundaryCondition/constant become fixed_species and stay forced", {
    skip_if_not(sbml_available())

    path <- file.path("..", "..", "examples", "sbml", "fixtures",
                      "b2_contorno.xml")
    crn <- suppressWarnings(read_sbml_to_dnar(path))
    expect_equal(crn$fixed_species, "S1")
    expect_true(crn$species_info$S1$boundaryCondition)
    expect_false(crn$species_info$S1$constant)

    forced <- sbml_forced_concentrations(crn)
    expect_equal(names(forced), "S1")
    expect_equal(forced$S1(0), 1)
    expect_equal(forced$S1(10), 1)

    out <- react4(
        species = crn$species,
        ci = crn$ci,
        reactions = crn$reactions,
        ki = crn$ki,
        t = seq(0, 1, length.out = 11),
        forced_concentrations = forced
    )
    expect_true(all(abs(out$S1 - 1) < 1e-10))
    expect_gt(tail(out$S2, 1), 0.5)

    multi <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="multi"><listOfCompartments>',
        '<compartment id="c1" spatialDimensions="3" size="1" constant="true"/>',
        '<compartment id="c2" spatialDimensions="3" size="1" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c1" initialConcentration="1"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies></model></sbml>'
    ), multi)
    expect_error(read_sbml_to_dnar(multi), "não suportado: múltiplos compartimentos")
})

test_that("0D compartment keeps amount without dividing by volume", {
    skip_if_not(sbml_available())

    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="zero_d"><listOfCompartments>',
        '<compartment id="c" spatialDimensions="0" size="2" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialAmount="4"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="true" constant="false"/>',
        '<species id="S2" compartment="c" initialAmount="0"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies><listOfParameters>',
        '<parameter id="k" value="1" constant="true"/>',
        '</listOfParameters><listOfReactions>',
        '<reaction id="r1" reversible="false"><listOfReactants>',
        '<speciesReference species="S1" stoichiometry="1" constant="true"/>',
        '</listOfReactants><listOfProducts>',
        '<speciesReference species="S2" stoichiometry="1" constant="true"/>',
        '</listOfProducts><kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">',
        '<apply><times/><ci>k</ci><ci>S1</ci></apply></math></kineticLaw>',
        '</reaction></listOfReactions></model></sbml>'
    ), path)
    crn <- suppressWarnings(read_sbml_to_dnar(path))
    expect_equal(unname(crn$ci[["S1"]]), 4)
    expect_equal(crn$fixed_species, "S1")
    expect_equal(crn$species_info$S1$spatialDimensions, 0)
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
    expect_equal(crn$outputs$compartment, 2)
    expect_equal(crn$outputs$k, 1)
    expect_false("S1" %in% names(crn$outputs))
    expect_equal(crn$species_info$S1$volume, 2)
    expect_equal(crn$species_info$S2$volume, 2)
    # Law stays amount/time (V*k at S=1 => 2); dx_scale is CF/V = 1/2.
    expect_equal(crn$ki, 2)
    expect_equal(unname(crn$conversion_factors[c("S1", "S2")]), c(0.5, 0.5))
    out <- react4(
        species = crn$species,
        ci = crn$ci,
        reactions = crn$reactions,
        ki = crn$ki,
        t = c(0, 1),
        dx_scale = crn$conversion_factors
    )
    expect_equal(out$S1[2], exp(-1), tolerance = 1e-4)

    stoich <- file.path("..", "..", "examples", "sbml", "fixtures",
                        "b5_species_reference.xml")
    crn_stoich <- suppressWarnings(read_sbml_to_dnar(stoich))
    expect_equal(crn_stoich$symbols$S1_stoich, 1)
    expect_false("S1_stoich" %in% names(crn_stoich$outputs))
    expect_equal(crn_stoich$ki, 0.5)

    env <- sbml_math_env(list(x = 1.5, flag = TRUE))
    expect_equal(eval(parse(text = "ceil(x)"), envir = env), 2)
    expect_equal(eval(parse(text = "piecewise(10, flag, 0)"), envir = env), 10)
    expect_equal(eval(parse(text = "piecewise(42)"), envir = env), 42)
    expect_equal(eval(parse(text = "piecewise(1, FALSE, 2, FALSE, 7)"), envir = env), 7)
})

test_that("conversionFactor becomes dx_scale without changing mass-action order", {
    skip_if_not(sbml_available())

    path <- file.path("..", "..", "examples", "sbml", "fixtures",
                      "b6_conversion_factor.xml")
    crn <- suppressWarnings(read_sbml_to_dnar(path))
    expect_equal(unname(crn$ci[["S1"]]), 2)
    # V = 1, so dx_scale is just the conversion factor.
    expect_equal(crn$conversion_factors[["S1"]], 2)
    expect_equal(crn$conversion_factors[["S2"]], 0.5)
    expect_equal(crn$reactions, "S1 -> S2")
    expect_equal(crn$outputs$s1conv, 2)
    expect_equal(crn$outputs$modelconv, 0.5)

    out <- react4(
        species = crn$species,
        ci = crn$ci,
        reactions = crn$reactions,
        ki = crn$ki,
        t = c(0, 1),
        dx_scale = crn$conversion_factors
    )
    # Mass-action S1->S2 with ki=1 and CF_S1=2 => dS1/dt = -2*S1.
    expect_equal(out$S1[2], 2 * exp(-2), tolerance = 1e-4)
})

test_that("hasOnlySubstanceUnits keeps amount state in ci", {
    skip_if_not(sbml_available())

    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="hosu"><listOfCompartments>',
        '<compartment id="c" spatialDimensions="3" size="2" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialAmount="8"',
        ' hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>',
        '<species id="S2" compartment="c" initialAmount="4"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies></model></sbml>'
    ), path)
    crn <- read_sbml_to_dnar(path)
    expect_equal(unname(crn$ci[["S1"]]), 8)
    expect_equal(unname(crn$ci[["S2"]]), 2)
    expect_true(crn$species_info$S1$hasOnlySubstanceUnits)
    expect_false(crn$species_info$S2$hasOnlySubstanceUnits)
    expect_equal(crn$conversion_factors[["S1"]], 1)
    expect_equal(crn$conversion_factors[["S2"]], 0.5)
})

test_that("hosu + initialConcentration multiplies by volume", {
    skip_if_not(sbml_available())

    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="hosu_conc"><listOfCompartments>',
        '<compartment id="c" spatialDimensions="3" size="2" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialConcentration="3"',
        ' hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies></model></sbml>'
    ), path)
    crn <- read_sbml_to_dnar(path)
    expect_equal(unname(crn$ci[["S1"]]), 6)
    expect_true(crn$species_info$S1$hasOnlySubstanceUnits)
})

test_that("initialAmount uses compartment size from initialAssignment", {
    skip_if_not(sbml_available())

    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="ia_vol"><listOfCompartments>',
        '<compartment id="c" spatialDimensions="3" size="1" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialAmount="8"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies><listOfInitialAssignments>',
        '<initialAssignment symbol="c"><math xmlns="http://www.w3.org/1998/Math/MathML">',
        '<cn>4</cn></math></initialAssignment>',
        '</listOfInitialAssignments></model></sbml>'
    ), path)
    crn <- read_sbml_to_dnar(path)
    expect_equal(crn$symbols$c, 4)
    expect_equal(unname(crn$ci[["S1"]]), 2)
})

test_that("IA formulas see concentration from initialAmount (table example)", {
    # Tabela: c size=2, S1 initialAmount=8 → conc=4; IA p := S1 → p=4 (não 8).
    skip_if_not(sbml_available())

    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="ia_species_conc"><listOfCompartments>',
        '<compartment id="c" spatialDimensions="3" size="2" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialAmount="8"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies><listOfParameters>',
        '<parameter id="p" constant="true"/>',
        '</listOfParameters><listOfInitialAssignments>',
        '<initialAssignment symbol="p"><math xmlns="http://www.w3.org/1998/Math/MathML">',
        '<ci>S1</ci></math></initialAssignment>',
        '</listOfInitialAssignments></model></sbml>'
    ), path)
    crn <- read_sbml_to_dnar(path)
    expect_equal(unname(crn$ci[["S1"]]), 4)
    expect_equal(crn$symbols$p, 4)
    expect_equal(crn$symbols$S1, 4)
})

test_that("IA-only compartment size still converts initialAmount", {
    skip_if_not(sbml_available())

    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="ia_only_vol"><listOfCompartments>',
        '<compartment id="c" spatialDimensions="3" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialAmount="8"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies><listOfInitialAssignments>',
        '<initialAssignment symbol="c"><math xmlns="http://www.w3.org/1998/Math/MathML">',
        '<cn>4</cn></math></initialAssignment>',
        '</listOfInitialAssignments></model></sbml>'
    ), path)
    crn <- read_sbml_to_dnar(path)
    expect_equal(crn$symbols$c, 4)
    expect_equal(unname(crn$ci[["S1"]]), 2)
})

test_that("non-finite kinetic law rates error instead of becoming 0", {
    skip_if_not(sbml_available())

    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="nan_rate"><listOfCompartments>',
        '<compartment id="c" spatialDimensions="3" size="1" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialConcentration="1"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '<species id="S2" compartment="c" initialConcentration="0"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies><listOfParameters>',
        '<parameter id="k" value="NaN" constant="true"/>',
        '</listOfParameters><listOfReactions>',
        '<reaction id="r1" reversible="false"><listOfReactants>',
        '<speciesReference species="S1" stoichiometry="1" constant="true"/>',
        '</listOfReactants><listOfProducts>',
        '<speciesReference species="S2" stoichiometry="1" constant="true"/>',
        '</listOfProducts><kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">',
        '<apply><plus/><ci>k</ci><ci>S1</ci></apply></math></kineticLaw>',
        '</reaction></listOfReactions></model></sbml>'
    ), path)
    crn <- suppressWarnings(read_sbml_to_dnar(path))
    expect_true(is.function(crn$ki[[1]]))
    expect_error(crn$ki[[1]](0, c(S1 = 1, S2 = 0)), "Non-finite")
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

    crn_ia <- read_sbml_to_dnar(ia)
    expect_equal(crn_ia$symbols$p1, 2)
    expect_equal(crn_ia$symbols$p2, 6)
    expect_equal(unname(crn_ia$ci[["S1"]]), 6)
    expect_equal(crn_ia$ki, 1)
})

test_that("rules are exposed as metadata without being applied", {
    skip_if_not(sbml_available())

    case_00296 <- file.path(
        "..", "..", "external", "sbml-test-suite", "cases", "semantic",
        "00296", "00296-sbml-l3v2.xml"
    )
    skip_if_not(file.exists(case_00296))

    crn <- suppressWarnings(read_sbml_to_dnar(case_00296))
    expect_gt(length(crn$rules), 0L)
    expect_equal(crn$rules[[1]]$type, "assignment")
    expect_true(nzchar(crn$rules[[1]]$variable))
    expect_true(nzchar(crn$rules[[1]]$formula))
    # Not applied: species initials come from the SBML species attributes.
    expect_true(all(is.finite(crn$ci)))
})

test_that("initialAssignment Kahn order detects cycles", {
    expect_error(
        sbml_resolve_initial_assignments(
            list(
                list(symbol = "a", formula = "b"),
                list(symbol = "b", formula = "a")
            ),
            list()
        ),
        "cíclico"
    )
    resolved <- sbml_resolve_initial_assignments(
        list(
            list(symbol = "k2", formula = "k1 * 3"),
            list(symbol = "k1", formula = "2"),
            list(symbol = "S1", formula = "k2")
        ),
        list()
    )
    expect_equal(resolved$k1, 2)
    expect_equal(resolved$k2, 6)
    expect_equal(resolved$S1, 6)
})

test_that("speciesReference stoichiometry from initialAssignment is applied", {
    skip_if_not(sbml_available())

    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="ia_stoich"><listOfCompartments>',
        '<compartment id="c" spatialDimensions="3" size="1" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialAmount="1"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '<species id="S2" compartment="c" initialAmount="0"',
        ' hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies><listOfParameters>',
        '<parameter id="p1" value="1" constant="true"/>',
        '<parameter id="k" value="0.3" constant="true"/>',
        '</listOfParameters><listOfInitialAssignments>',
        '<initialAssignment symbol="ref1"><math xmlns="http://www.w3.org/1998/Math/MathML">',
        '<apply><times/><cn>2</cn><ci>p1</ci></apply></math></initialAssignment>',
        '</listOfInitialAssignments><listOfReactions>',
        '<reaction id="r1" reversible="false"><listOfReactants>',
        '<speciesReference id="ref1" species="S1" constant="true"/>',
        '</listOfReactants><listOfProducts>',
        '<speciesReference species="S2" stoichiometry="1" constant="true"/>',
        '</listOfProducts><kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">',
        '<apply><times/><ci>k</ci><ci>S1</ci></apply></math></kineticLaw>',
        '</reaction></listOfReactions></model></sbml>'
    ), path)
    crn <- suppressWarnings(read_sbml_to_dnar(path))
    expect_equal(crn$reactions, "2S1 -> S2")
    expect_equal(crn$symbols$ref1, 2)
})
