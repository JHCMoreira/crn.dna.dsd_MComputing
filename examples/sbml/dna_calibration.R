# Calibração de qmax e cmax do circuito DNA (papel D).
#
# Para cada caso, lê o SBML uma vez e simula o circuito numa grade de fatores:
#   qmax = fator_q * max(ki)       cmax = fator_c * beta * material
# e mede o erro em relação à referência da SBML Test Suite.
#
# Uso, na raiz do repositório (depois de run_dna_suite.R):
#   Rscript examples/sbml/dna_calibration.R
#   Rscript examples/sbml/dna_calibration.R --cases=00001,00062 --factors=1e2,1e4,1e6
#   Rscript examples/sbml/dna_calibration.R --so-graficos  # refaz só os gráficos
#
# Casos: por padrão, os de results/dna_summary.csv em que a CRN passou
# (assim o erro medido é só o da aproximação do DNA).
#
# Saídas em results/:
#   calibracao.csv                uma linha por caso x fator_q x fator_c
#   calibracao_erro.png           NRMSE x fator (qmax e cmax juntos, só qmax, só cmax)
#   calibracao_aprovacao.png      % de casos aprovados x fator
#   calibracao_grade.png          % aprovados na grade fator_q x fator_c

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
suite <- option("suite", file.path("external", "sbml-test-suite", "cases", "semantic"))
out_dir <- option("out", "results")
factors <- as.numeric(strsplit(option("factors", "1e2,1e3,1e4,1e5,1e6"), ",")[[1]])
reference_factor <- as.numeric(option("fixed", "1e4"))
cases <- option("cases", "")
plots_only <- "--so-graficos" %in% args

source(file.path("examples", "sbml", "testbench_dna.R"))

if (plots_only) {
    case_ids <- character()
} else if (nzchar(cases)) {
    case_ids <- strsplit(cases, ",", fixed = TRUE)[[1]]
} else {
    summary_file <- file.path(out_dir, "dna_summary.csv")
    if (!file.exists(summary_file)) {
        stop("Rode antes examples/sbml/run_dna_suite.R (falta ", summary_file, ")")
    }
    s <- read.csv(summary_file, colClasses = c(case_id = "character"))
    case_ids <- s$case_id[s$status_crn %in% "passou"]
}

run_calibration <- function(case_ids, factors) {
    rows <- list()
    discovered <- discover_sbml_cases(suite, case_ids)
    for (i in seq_len(nrow(discovered))) {
        id <- discovered$case_id[[i]]
        loaded <- tryCatch(
            load_sbml_case_dna(id, suite, case_dir = discovered$case_dir[[i]]),
            error = function(e) NULL)
        if (is.null(loaded) || is.null(loaded$context)) {
            cat(sprintf("[%d/%d] %s: pulado (%s)\n", i, nrow(discovered), id,
                        if (is.null(loaded)) "erro ao ler" else loaded$result$reason))
            next
        }
        for (fq in factors) {
            for (fc in factors) {
                r <- simulate_case_dna(loaded, qmax_factor = fq, cmax_factor = fc,
                                       stiff_retry = FALSE)
                rows[[length(rows) + 1L]] <- data.frame(
                    case_id = id, fator_qmax = fq, fator_cmax = fc,
                    status_dna = r$status_dna,
                    rigido = grepl("sistema rigido", r$reason),
                    nrmse = if (is.null(r$metrics)) NA_real_ else r$metrics$nrmse_dna,
                    erro_rel_max = if (is.null(r$metrics)) NA_real_ else r$metrics$max_rel_err_dna,
                    segundos = r$sim_seconds,
                    stringsAsFactors = FALSE)
            }
        }
        cat(sprintf("[%d/%d] %s\n", i, nrow(discovered), id))
    }
    do.call(rbind, rows)
}

if (plots_only) {
    cal <- read.csv(file.path(out_dir, "calibracao.csv"), colClasses = c(case_id = "character"))
    factors <- sort(unique(cal$fator_qmax))
} else {
    cat("Casos:", length(case_ids), "| fatores:", paste(factors, collapse = ", "), "\n")
    started <- Sys.time()
    cal <- run_calibration(case_ids, factors)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(cal, file.path(out_dir, "calibracao.csv"), row.names = FALSE)
    cat("Tempo:", format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
}

# ---- Resumos ---------------------------------------------------------------

# Três varreduras tiradas da grade: os dois juntos, só qmax e só cmax.
sweeps <- rbind(
    transform(subset(cal, fator_qmax == fator_cmax),
              varredura = "qmax e cmax juntos", fator = fator_qmax),
    transform(subset(cal, fator_cmax == reference_factor),
              varredura = "só qmax", fator = fator_qmax),
    transform(subset(cal, fator_qmax == reference_factor),
              varredura = "só cmax", fator = fator_cmax)
)
sweep_levels <- unique(sweeps$varredura)
sweeps$varredura <- factor(sweeps$varredura, levels = sweep_levels)

by_sweep <- do.call(rbind, lapply(split(sweeps, list(sweeps$varredura, sweeps$fator), drop = TRUE),
                                  function(d) data.frame(
    varredura = d$varredura[[1]], fator = d$fator[[1]], casos = nrow(d),
    aprovados = mean(d$status_dna == "passou"),
    rigidos = mean(d$rigido),
    nrmse_mediana = median(d$nrmse, na.rm = TRUE),
    nrmse_p10 = unname(quantile(d$nrmse, 0.1, na.rm = TRUE)),
    nrmse_p90 = unname(quantile(d$nrmse, 0.9, na.rm = TRUE)),
    segundos_mediana = median(d$segundos, na.rm = TRUE),
    segundos_max = max(d$segundos, na.rm = TRUE))))
by_sweep <- by_sweep[order(by_sweep$varredura, by_sweep$fator), ]
write.csv(by_sweep, file.path(out_dir, "calibracao_resumo.csv"), row.names = FALSE)
print(by_sweep, row.names = FALSE, digits = 3)

grid <- do.call(rbind, lapply(split(cal, list(cal$fator_qmax, cal$fator_cmax), drop = TRUE),
                              function(d) data.frame(
    fator_qmax = d$fator_qmax[[1]], fator_cmax = d$fator_cmax[[1]],
    aprovados = mean(d$status_dna == "passou"), rigidos = mean(d$rigido))))

# ---- Gráficos --------------------------------------------------------------
# Paleta de referência (dataviz): 3 primeiras cores categóricas, tinta e grade.

library(ggplot2)
series_colors <- setNames(c("#257de0", "#f47230", "#00b894"), sweep_levels)
# Legendas com potências desenhadas pelo gráfico (plotmath).
fixed_power <- round(log10(reference_factor))
series_labels <- c(expression("qmax e cmax juntos"),
                   bquote("só qmax (cmax fixo em" ~ 10^.(fixed_power) * ")"),
                   bquote("só cmax (qmax fixo em" ~ 10^.(fixed_power) * ")"))
ink <- "#191919"; ink2 <- "#666666"; grid_line <- "#e7e7e7"
base_theme <- theme_minimal(base_size = 12) +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
          axis.title = element_text(colour = ink2), plot.subtitle = element_text(colour = ink2),
          panel.grid.major = element_line(colour = grid_line, linewidth = 0.3),
          panel.grid.minor = element_blank(), legend.position = "top",
          legend.title = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
pow10 <- function(x) parse(text = sprintf("10^%d", round(log10(as.numeric(x)))))
log_x <- scale_x_log10(breaks = factors, labels = pow10)

ok <- subset(by_sweep, is.finite(nrmse_mediana))
p_err <- ggplot(ok, aes(fator, nrmse_mediana, colour = varredura, fill = varredura)) +
    geom_ribbon(aes(ymin = nrmse_p10, ymax = nrmse_p90), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    scale_y_log10(labels = pow10) + log_x +
    scale_colour_manual(values = series_colors, labels = series_labels) +
    scale_fill_manual(values = series_colors, labels = series_labels) +
    labs(x = "fator sobre max(ki) (qmax) e sobre o material (cmax)",
         y = "NRMSE do circuito DNA (mediana; faixa 10%-90%)",
         title = "Erro do circuito DNA em função de qmax e cmax",
         subtitle = sprintf("%d casos da SBML Test Suite em que a CRN passa; simulações rígidas ficam de fora",
                            length(unique(cal$case_id)))) +
    base_theme
ggsave(file.path(out_dir, "calibracao_erro.png"), p_err, width = 8, height = 5, dpi = 150)

p_pass <- ggplot(by_sweep, aes(fator, 100 * aprovados, colour = varredura)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    log_x + scale_y_continuous(limits = c(0, 100)) +
    scale_colour_manual(values = series_colors, labels = series_labels) +
    labs(x = "fator", y = "% de casos aprovados no DNA",
         title = "Casos aprovados no circuito DNA (mesma tolerância da CRN)",
         subtitle = "Queda com fatores altos nos dois parâmetros: sistema rígido, o resolvedor não termina") +
    base_theme
ggsave(file.path(out_dir, "calibracao_aprovacao.png"), p_pass, width = 8, height = 5, dpi = 150)

pct <- function(x) ifelse(x > 0 & x < 0.01, "<1%", sprintf("%.0f%%", 100 * x))
grid$rotulo <- paste0(pct(grid$aprovados),
                      ifelse(grid$rigidos > 0, paste0("\n(", pct(grid$rigidos), " rígidos)"), ""))
p_grid <- ggplot(grid, aes(factor(fator_qmax), factor(fator_cmax), fill = 100 * aprovados)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = rotulo, colour = aprovados > 0.6), size = 3.2, lineheight = 0.9,
              show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "#ffffff", `FALSE` = ink)) +
    scale_fill_gradient(low = "#eaf2fc", high = "#0d4e9c", limits = c(0, 100), name = "% aprovados") +
    scale_x_discrete(labels = pow10) + scale_y_discrete(labels = pow10) +
    labs(x = "fator de qmax", y = "fator de cmax",
         title = "% de casos aprovados no DNA para cada par (qmax, cmax)") +
    base_theme + theme(legend.position = "right", legend.title = element_text(colour = ink2),
                       panel.grid.major = element_blank())
ggsave(file.path(out_dir, "calibracao_grade.png"), p_grid, width = 7, height = 5.5, dpi = 150)

cat("Gráficos salvos em", out_dir, "\n")
