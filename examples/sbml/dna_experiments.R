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
