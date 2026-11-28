# Gráficos de um caso da SBML Test Suite: referência x CRN x circuito DNA.
#
# Uso, na raiz do repositório:
#   Rscript examples/sbml/dna_plots.R                         # casos do relatório
#   Rscript examples/sbml/dna_plots.R --cases=00062,00007 --out=results/plots
#
# Também pode ser carregado com source() para usar plot_case_comparison().

library(ggplot2)

# Paleta de referência (dataviz): cores categóricas em ordem fixa, tintas e grade.
DNA_PLOT_COLORS <- c("#257de0", "#f47230", "#00b894", "#f2a600",
                     "#ee7fac", "#008a17", "#6a4ce2", "#ec4647")
DNA_PLOT_INK <- "#191919"
DNA_PLOT_INK2 <- "#666666"
DNA_PLOT_GRID <- "#e7e7e7"

#' Plot reference, CRN and DNA series of one case, with the DNA error below.
#'
#' Top panel: SBML reference as points, CRN as solid lines and DNA circuit as
#' dashed lines, one color per species. Bottom panel: DNA minus reference,
#' which shows the deviation that is invisible at the scale of the top panel.
#'
#' @param result Output of `run_sbml_case_dna()` with `reference`, `crn`, `dna`.
#' @param file PNG path; when `NULL` the plot is only returned.
#' @return The ggplot object.
plot_case_comparison <- function(result, file = NULL) {
    if (is.null(result$dna)) {
        stop("caso ", result$case_id, " sem série do DNA (", result$status_dna, ")")
    }
    species <- setdiff(intersect(names(result$dna), names(result$reference)), "time")
    if (length(species) > length(DNA_PLOT_COLORS)) {
        stop("mais de ", length(DNA_PLOT_COLORS), " variáveis: use small multiples")
    }
    panel_top <- "concentração"
    panel_err <- "erro do DNA (DNA − referência)"
    long <- function(df, source, panel) {
        data.frame(time = rep(df$time, length(species)),
                   species = rep(species, each = nrow(df)),
                   value = unlist(df[species], use.names = FALSE),
                   source = source, panel = panel, stringsAsFactors = FALSE)
    }
    err <- result$dna
    err[species] <- result$dna[species] - result$reference[species]
    data <- rbind(long(result$reference, "referência SBML", panel_top),
                  long(result$crn, "CRN (DNAr)", panel_top),
                  long(result$dna, "circuito DNA", panel_top),
                  long(err, "circuito DNA", panel_err))
    data$species <- factor(data$species, levels = species)
    data$panel <- factor(data$panel, levels = c(panel_top, panel_err))
    data$source <- factor(data$source,
                          levels = c("referência SBML", "CRN (DNAr)", "circuito DNA"))
    colors <- setNames(DNA_PLOT_COLORS[seq_along(species)], species)

    m <- result$metrics
    subtitle <- sprintf(
        "DNA: %s | NRMSE %.1e | erro relativo máx. %.1e | circuito com %d espécies e %d reações | qmax = %.2g, cmax = %.2g",
        result$status_dna, m$nrmse_dna, m$max_rel_err_dna,
        result$n_dna_species, result$n_dna_reactions,
        result$params$qmax, result$params$cmax)

    p <- ggplot(data, aes(time, value, colour = species)) +
        geom_hline(data = data.frame(panel = factor(panel_err, levels = levels(data$panel)),
                                     y = 0),
                   aes(yintercept = y), colour = "#c2c2c2", linewidth = 0.4) +
        geom_line(data = subset(data, source == "CRN (DNAr)"),
                  aes(linetype = source), linewidth = 0.8) +
        geom_line(data = subset(data, source == "circuito DNA"),
                  aes(linetype = source), linewidth = 0.8) +
        geom_point(data = subset(data, source == "referência SBML"),
                   aes(shape = source), size = 2.2, stroke = 0.6,
                   fill = "white") +
        facet_grid(panel ~ ., scales = "free_y", switch = "y") +
        scale_colour_manual(values = colors, name = "variável") +
        scale_linetype_manual(values = c("CRN (DNAr)" = "solid", "circuito DNA" = "22"),
                              name = NULL) +
        scale_shape_manual(values = c("referência SBML" = 21), name = NULL) +
        labs(x = "tempo", y = NULL,
             title = sprintf("Caso %s: referência SBML x CRN x circuito DNA", result$case_id),
             subtitle = subtitle) +
        theme_minimal(base_size = 11) +
        theme(text = element_text(colour = DNA_PLOT_INK),
              axis.text = element_text(colour = DNA_PLOT_INK2),
              axis.title = element_text(colour = DNA_PLOT_INK2),
              plot.subtitle = element_text(colour = DNA_PLOT_INK2, size = 8.5),
              strip.text = element_text(colour = DNA_PLOT_INK2),
              strip.placement = "outside",
              panel.grid.major = element_line(colour = DNA_PLOT_GRID, linewidth = 0.3),
              panel.grid.minor = element_blank(),
              panel.spacing = unit(1, "lines"),
              legend.position = "right",
              plot.background = element_rect(fill = "white", colour = NA)) +
        guides(colour = guide_legend(order = 1),
               shape = guide_legend(order = 2),
               linetype = guide_legend(order = 3, override.aes = list(colour = DNA_PLOT_INK2)))
    if (!is.null(file)) {
        dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
        ggsave(file, p, width = 9, height = 6.5, dpi = 150)
    }
    p
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    option <- function(name, default) {
        hit <- grep(paste0("^--", name, "="), args, value = TRUE)
        if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
    }
    suite <- option("suite", file.path("external", "sbml-test-suite", "cases", "semantic"))
    out_dir <- option("out", file.path("results", "plots"))
    # 00062: reversível bimolecular; 00007: espécie de contorno; 01232: formação 0 -> X
    cases <- strsplit(option("cases", "00062,00007,01232"), ",", fixed = TRUE)[[1]]

    source(file.path("examples", "sbml", "testbench_dna.R"))
    for (id in cases) {
        r <- run_sbml_case_dna(id, suite)
        if (is.null(r$dna)) {
            cat(id, ": sem gráfico (", r$status_dna, ":", r$reason, ")\n")
            next
        }
        plot_case_comparison(r, file.path(out_dir, paste0(id, ".png")))
        cat(id, ": ", file.path(out_dir, paste0(id, ".png")), "\n", sep = "")
    }
}
