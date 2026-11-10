library(testthat)

testbench_path <- file.path("examples", "sbml", "testbench.R")
if (!file.exists(testbench_path)) {
    testbench_path <- file.path("..", "..", "examples", "sbml", "testbench.R")
}
source(testbench_path)

test_that("settings and metadata from the semantic suite are parsed", {
    case_dir <- file.path(tempdir(), "00001")
    dir.create(case_dir, showWarnings = FALSE)
    settings <- file.path(case_dir, "00001-settings.txt")
    metadata <- file.path(case_dir, "00001-model.m")
    writeLines(c(
        "start: 0", "duration: 5", "steps: 50",
        "variables: S1, S2", "absolute: 1e-7", "relative: 1e-4",
        "amount: S1", "concentration: S2"
    ), settings)
    writeLines(c(
        "(*", "category: Test", "synopsis: Simple reaction",
        "componentTags: Compartment, Species, Reaction",
        "testTags: Amount, NonUnityCompartment", "testType: TimeCourse",
        "levels: 2.4, 3.2", "generatedBy: Analytic", "packagesPresent:", "*)"
    ), metadata)

    parsed_settings <- parse_sbml_settings(settings)
    parsed_metadata <- parse_sbml_metadata(metadata)
    expect_equal(parsed_settings$steps, 50L)
    expect_equal(parsed_settings$variables, c("S1", "S2"))
    expect_equal(parsed_settings$absolute, 1e-7)
    expect_equal(parsed_metadata$component_tags, c("Compartment", "Species", "Reaction"))
    expect_equal(parsed_metadata$levels, c("2.4", "3.2"))
})

test_that("comparison uses the official combined pointwise tolerance", {
    reference <- data.frame(time = 0, S1 = 100)
    within <- data.frame(time = 0, S1 = 100.0105)
    outside <- data.frame(time = 0, S1 = 100.012)

    expect_true(compare_sbml_results(
        within, reference, "S1", 1e-3, 1e-4
    )$passed)
    failed <- compare_sbml_results(outside, reference, "S1", 1e-3, 1e-4)
    expect_false(failed$passed)
    expect_equal(failed$failing_points, 1L)
})

test_that("non-finite reference values are compared explicitly", {
    reference <- data.frame(time = 0:2, S1 = c(NaN, Inf, -Inf))
    same <- data.frame(time = 0:2, S1 = c(NaN, Inf, -Inf))
    wrong <- data.frame(time = 0:2, S1 = c(0, Inf, Inf))
    expect_true(compare_sbml_results(same, reference, "S1")$passed)
    expect_equal(compare_sbml_results(wrong, reference, "S1")$failing_points, 2L)
})

test_that("nested cases are discovered and the preferred SBML variant is selected", {
    root <- file.path(tempdir(), "semantic-suite")
    case_dir <- file.path(root, "00042")
    dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("start: 0", file.path(case_dir, "00042-settings.txt"))
    writeLines("<sbml/>", file.path(case_dir, "00042-sbml-l2v4.xml"))
    writeLines("<sbml/>", file.path(case_dir, "00042-sbml-l3v2.xml"))

    cases <- discover_sbml_cases(root)
    model <- find_sbml_case_files("00042", cases$case_dir, "l3v2")
    expect_equal(cases$case_id, "00042")
    expect_equal(model$variant, "l3v2")
    expect_equal(model$level, 3L)
    expect_equal(model$version, 2L)
})

test_that("known unsupported SBML constructs are explicit", {
    metadata <- empty_sbml_metadata()
    metadata$test_type <- "TimeCourse"
    metadata$component_tags <- c("Species", "RateRule")
    expect_match(sbml_known_unsupported(metadata), "RateRule")
    expect_equal(sbml_failure_cause(metadata), "regras")
})

test_that("unsupported output variables retain the complete case record", {
    root <- tempfile("missing-output-")
    dir.create(root)
    writeLines(c("start: 0", "duration: 1", "steps: 1",
                 "variables: parameter1", "absolute: 1e-6", "relative: 1e-4"),
               file.path(root, "00001-settings.txt"))
    writeLines(c("testType: TimeCourse", "componentTags: Parameter"),
               file.path(root, "00001-model.m"))
    writeLines("<sbml/>", file.path(root, "00001-sbml-l3v2.xml"))
    writeLines("time,parameter1\n0,1\n1,1", file.path(root, "00001-results.csv"))
    isolated <- new.env(parent = environment(run_sbml_case))
    isolated$sbml_available <- function() TRUE
    isolated$validate_sbml <- function(path) list(valid = TRUE, errors = list())
    isolated$read_sbml_to_dnar <- function(path) list(species = "S1")
    runner <- run_sbml_case
    environment(runner) <- isolated
    result <- runner("00001", root, case_dir = root)
    expect_equal(result$case_id, "00001")
    expect_equal(result$status, "nao suportado")
    expect_match(result$reason, "parameter1")
    expect_equal(nrow(sbml_result_row(result)), 1L)
})

test_that("different time grids and missing numerical values cannot pass", {
    reference <- data.frame(time = 0:1, S1 = c(0, 1))
    expect_error(compare_sbml_results(data.frame(time = 1:2, S1 = c(0, 1)), reference), "Tempos")
    actual <- data.frame(time = 0:1, S1 = c(NA_real_, 1))
    expect_false(compare_sbml_results(actual, reference)$passed)
})

test_that("output amounts and concentrations use compartment size", {
    skip_if_not(sbml_available())
    path <- tempfile(fileext = ".xml")
    writeLines(c(
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">',
        '<model id="units"><listOfCompartments>',
        '<compartment id="c" size="2" constant="true"/>',
        '</listOfCompartments><listOfSpecies>',
        '<species id="S1" compartment="c" initialConcentration="3" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>',
        '<species id="S2" compartment="c" initialAmount="8" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>',
        '</listOfSpecies></model></sbml>'
    ), path)
    settings <- list(variables = c("S1", "S2"), amount = "S1", concentration = "S2")
    output <- sbml_project_output(data.frame(time = 0, S1 = 3, S2 = 8), settings, path)
    expect_equal(output$S1, 6)
    expect_equal(output$S2, 4)
})

test_that("infrastructure failures with unknown SBML version remain in totals", {
    item <- list(case_id = "00001", sbml_level = NA_integer_,
                 sbml_version = NA_integer_, sbml_variant = "",
                 status = "falhou", cause = "infraestrutura_do_runner",
                 validation = NULL, comparison = NULL,
                 settings = list(absolute = NA_real_, relative = NA_real_),
                 metadata = empty_sbml_metadata(), warnings = character(),
                 reason = "missing model", xml_file = "")
    suite <- summarize_sbml_cases(list(item))
    expect_equal(sum(suite$status_counts$count), 1L)
    expect_equal(sum(suite$compatibility$count), 1L)
    expect_equal(suite$compatibility$sbml_level, "unknown")
})
