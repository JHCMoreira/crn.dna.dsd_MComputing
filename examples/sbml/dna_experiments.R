# Experimentos do papel D: o circuito DNA (4-domain, Soloveichik 2010)
# reproduz a CRN? Qual regra de reversão de escala vale?
# Uso, na raiz do projeto: Rscript examples/sbml/dna_experiments.R
source("R/parser.R"); source("R/util_functions.R")
source("R/crn_reactor.R"); source("R/4domain_reactor.R")
dir.create("results", showWarnings = FALSE)

# react3() tem um bug (usa "..." sem declarar, R/crn_reactor.R:1001);
# por isso a CRN de referência usa react4(), o mesmo motor do react_4domain().
erro_rel <- function(a, b) max(abs(a - b)) / max(abs(b))

# ------------------------------------------------------------------
# Experimento 1: A + B -> C (parâmetros do teste existente)
# ------------------------------------------------------------------
sp <- c("A", "B", "C"); ci <- c(1e3, 1e3, 0)
rx <- "A + B -> C";      k  <- 1e-7
t  <- seq(0, 72000, 100)

# 1a) CRN pura e circuito DNA sem escala
crn <- as.data.frame(react4(sp, ci, rx, k, t))
dna <- react_4domain(sp, ci, rx, k, qmax = 1e-3, cmax = 1e5,
                     alpha = 1, beta = 1, t = t)
b <- as.data.frame(dna$behavior)
cat("Colunas do behavior:", paste(names(b), collapse = ", "), "\n")
erro1_A <- erro_rel(b$A, crn$A)
erro1_C <- erro_rel(b$C, crn$C)
cat(sprintf("1a) sem escala:  erro rel. máx  A = %.3g | C = %.3g\n", erro1_A, erro1_C))

# 1b) Com escala: alpha = 2 (2x mais lento), beta = 10 (10x mais concentrado)
alpha <- 2; beta <- 10
dna2 <- react_4domain(sp, ci, rx, k, qmax = 1e-3, cmax = 1e6,
                      alpha = alpha, beta = beta, t = t * alpha)
b2 <- as.data.frame(dna2$behavior)
erro2_A <- erro_rel(b2$A / beta, crn$A)
erro2_C <- erro_rel(b2$C / beta, crn$C)
cat(sprintf("1b) alpha=2, beta=10, revertido: erro A = %.3g | C = %.3g\n", erro2_A, erro2_C))

png("results/exp1_ApB.png", width = 900, height = 500)
plot(crn$time, crn$C, type = "l", lwd = 2, ylim = c(0, max(ci)),
     xlab = "tempo (s)", ylab = "concentração",
     main = "A + B -> C: CRN x circuito DNA")
lines(crn$time, crn$A, lwd = 2, col = "grey40")
lines(b$time,  b$C,  lty = 2, lwd = 2, col = "red")
lines(b$time,  b$A,  lty = 2, lwd = 2, col = "orange")
lines(t, b2$C / beta, lty = 3, lwd = 2, col = "blue")
legend("right",
       c("CRN: C", "CRN: A", "DNA: C", "DNA: A",
         "DNA (alpha=2, beta=10) revertido: C"),
       lty = c(1, 1, 2, 2, 3), col = c("black", "grey40", "red", "orange", "blue"),
       lwd = 2)
dev.off()

# ------------------------------------------------------------------
# Experimento 2: lambda_1 importa?
# Duas reações bimoleculares com sigmas diferentes -> cria buffer.
# Varremos qmax como múltiplo de sigma_max: perto de 1 o lambda_1 é
# grande, mas o circuito também fica ruim (reações auxiliares lentas).
# ------------------------------------------------------------------
sp2 <- c("A", "B", "C"); ci2 <- c(1e3, 1e3, 0)
rx2 <- c("A + B -> C", "C + B -> A"); k2 <- c(1e-7, 5e-8)
sigma_max <- 1e-7
crn2 <- as.data.frame(react4(sp2, ci2, rx2, k2, t))
cat("2) Erro relativo máximo por espécie (A B C)\n")
cat(sprintf("%8s %8s | %-20s | %-20s\n", "q/sigma", "lambda_1",
            "sem dividir", "dividindo por lambda"))
exp2 <- list()
for (f in c(1.5, 2, 5, 11, 101, 1e3, 1e4)) {
  q   <- sigma_max * f
  buf <- get_buff_modules(check_crn(sp2, ci2, rx2, k2, t), k2, q, 1e5)
  lambda_1 <- if (is.null(buf)) 1 else buf$lambda_1
  d <- as.data.frame(react_4domain(sp2, ci2, rx2, k2, qmax = q, cmax = 1e5,
                                   alpha = 1, beta = 1, t = t)$behavior)
  sem <- sapply(sp2, function(s) erro_rel(d[[s]], crn2[[s]]))
  com <- sapply(sp2, function(s) erro_rel(d[[s]] / lambda_1, crn2[[s]]))
  cat(sprintf("%8g %8.4f | %.3f %.3f %.3f  | %.3f %.3f %.3f\n",
              f, lambda_1, sem[1], sem[2], sem[3], com[1], com[2], com[3]))
  exp2[[length(exp2) + 1]] <- data.frame(
    fator_q_sigma = f, lambda_1 = lambda_1, especie = sp2,
    erro_sem_dividir = sem, erro_dividindo = com)
}
write.csv(do.call(rbind, exp2), "results/exp2_lambda.csv", row.names = FALSE)

# ------------------------------------------------------------------
# Experimento 3: casos reais da SBML Test Suite pelo circuito DNA
# Fluxo à mão: SBML -> read_sbml_to_dnar() -> react3() (CRN) e
# react_4domain() (DNA) -> comparação com o -results.csv da suíte.
# Precisa da suíte em external/sbml-test-suite
# (Rscript scripts/download_sbml_test_suite.R, do papel A).
# ------------------------------------------------------------------
source("examples/sbml/testbench.R")   # read_sbml_to_dnar(), parse_sbml_settings()

suite_dir <- file.path("external", "sbml-test-suite", "cases", "semantic")
casos_exp3 <- c("00048", "00060", "00062")   # acrescente outros casos aqui

passa_tol <- function(sim, ref, especies, abs_tol, rel_tol) {
  all(vapply(especies, function(s)
    all(abs(sim[[s]] - ref[[s]]) <= abs_tol + rel_tol * abs(ref[[s]])),
    logical(1)))
}

# atol/rtol do deSolve: o padrão (1e-6) é do tamanho dos valores de vários
# casos da suíte (~1e-6) e gera ~1% de erro numérico (ex.: 00062). Usamos os
# mesmos valores da rodada do papel A.
rodar_caso_dna <- function(caso, fator_qmax = 1e4, fator_cmax = 1e4,
                           atol = 1e-12, rtol = 1e-10) {
  pasta <- file.path(suite_dir, caso)
  xml <- file.path(pasta, paste0(caso, "-sbml-l3v2.xml"))
  if (!file.exists(xml)) {                # sem L3V2: usa a versão mais nova
    xml <- tail(sort(list.files(pasta, "-sbml-l.*\\.xml$", full.names = TRUE)), 1)
  }
  avisos <- character()
  m <- withCallingHandlers(read_sbml_to_dnar(xml), warning = function(w) {
    avisos <<- c(avisos, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  settings_file <- file.path(pasta, paste0(caso, "-settings.txt"))
  s <- parse_sbml_settings(settings_file)
  linhas <- readLines(settings_file, warn = FALSE)
  valor <- function(chave) as.numeric(sub(".*:", "", grep(paste0("^", chave, ":"), linhas, value = TRUE)))
  tempos <- seq(s$start, s$start + s$duration, length.out = s$steps + 1)
  qmax <- max(m$ki) * fator_qmax
  cmax <- max(max(m$ci), 1e-12) * fator_cmax

  crn <- as.data.frame(react3(m$species, unname(m$ci), m$reactions, unname(m$ki), tempos,
                              atol = atol, rtol = rtol))
  d <- suppressWarnings(react_4domain(m$species, unname(m$ci), m$reactions, unname(m$ki),
                                      qmax = qmax, cmax = cmax, alpha = 1, beta = 1, t = tempos,
                                      atol = atol, rtol = rtol))
  dna <- as.data.frame(d$behavior)[, c("time", m$species)]
  ref <- read.csv(file.path(pasta, paste0(caso, "-results.csv")), check.names = FALSE)
  esp <- intersect(m$species, names(ref))
  erro_rel_dna <- max(vapply(esp, function(x) erro_rel(dna[[x]], ref[[x]]), numeric(1)))

  png(file.path("results", paste0("caso", caso, "_dna.png")), width = 900, height = 500)
  cores <- seq_along(esp)
  matplot(ref$time, ref[esp], pch = 16, col = cores, xlab = "tempo", ylab = "valor",
          main = paste("Caso", caso, ": pontos = referência, linha = CRN, tracejado = DNA"))
  matlines(crn$time, crn[esp], lty = 1, lwd = 2, col = cores)
  matlines(dna$time, dna[esp], lty = 2, lwd = 2, col = cores)
  legend("right", esp, col = cores, lwd = 2)
  invisible(dev.off())

  data.frame(
    caso = caso, reacoes = paste(m$reactions, collapse = "; "),
    avisos_leitura = length(avisos), qmax = qmax, cmax = cmax,
    reacoes_dna = length(d$reactions),
    erro_rel_max_dna = erro_rel_dna,
    crn_passa_tol_caso = passa_tol(crn, ref, esp, valor("absolute"), valor("relative")),
    dna_passa_tol_caso = passa_tol(dna, ref, esp, valor("absolute"), valor("relative")),
    dna_passa_tol_grupo = passa_tol(dna, ref, esp, 1e-6, 1e-4),
    dna_passa_tol_dna   = passa_tol(dna, ref, esp, 1e-6, 1e-2)
  )
}

exp3 <- do.call(rbind, lapply(casos_exp3, function(caso)
  tryCatch(rodar_caso_dna(caso), error = function(e) {
    cat("Caso", caso, "deu erro:", conditionMessage(e), "\n"); NULL })))
cat("3) Casos da SBML Test Suite pelo circuito DNA\n")
print(exp3[, setdiff(names(exp3), "reacoes")], row.names = FALSE)
write.csv(exp3, "results/exp3_casos_suite.csv", row.names = FALSE)
