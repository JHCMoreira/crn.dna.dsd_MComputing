# Classificacao de leis cineticas SBML (papel C).
#
# Pergunta respondida para cada reacao: "a lei e uma constante vezes a
# concentracao de cada reagente elevada a sua estequiometria?"
#   - sim -> "mass_action", com a constante k extraida;
#   - sim, mas a lei tem a forma f - g (reversivel) e cada metade e acao de
#     massas -> "reversible_mass_action", com k_forward e k_backward;
#   - nao -> "general", com uma razao especifica (vai para a tabela de suporte
#     do relatorio).
#
# A analise e feita sobre a arvore sintatica da formula (parse), nao por regex,
# entao ordem dos fatores, parenteses, `A^2` ou `A*A` dao o mesmo resultado.

# ---- nucleo: achata um produto de fatores -----------------------------------

# Devolve list(ok, const, powers, reason). `powers` e um vetor nomeado com o
# expoente de cada especie; `const` e o produto numerico dos fatores constantes.
sbml_flatten_product <- function(expr, species, env, power = 1) {
    fail <- function(reason) {
        list(ok = FALSE, const = NA_real_, powers = numeric(), reason = reason)
    }
    merge <- function(a, b) {
        if (!a$ok) return(a)
        if (!b$ok) return(b)
        powers <- a$powers
        for (s in names(b$powers)) {
            powers[s] <- (if (s %in% names(powers)) powers[[s]] else 0) +
                b$powers[[s]]
        }
        list(ok = TRUE, const = a$const * b$const, powers = powers, reason = NULL)
    }
    vars <- all.vars(expr)
    if (!any(vars %in% species)) {
        value <- tryCatch(eval(expr, envir = env), error = function(e) NULL)
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
            return(fail("fator constante nao finito (Inf/NaN) ou parametro desconhecido"))
        }
        return(list(ok = TRUE, const = as.numeric(value)^power,
                    powers = numeric(), reason = NULL))
    }
    if (is.symbol(expr)) {
        name <- as.character(expr)
        return(list(ok = TRUE, const = 1,
                    powers = setNames(power, name), reason = NULL))
    }
    op <- as.character(expr[[1]])
    if (op == "(") {
        return(sbml_flatten_product(expr[[2]], species, env, power))
    }
    if (op == "*") {
        return(merge(
            sbml_flatten_product(expr[[2]], species, env, power),
            sbml_flatten_product(expr[[3]], species, env, power)
        ))
    }
    if (op == "/") {
        return(merge(
            sbml_flatten_product(expr[[2]], species, env, power),
            sbml_flatten_product(expr[[3]], species, env, -power)
        ))
    }
    if (op == "^") {
        if (any(all.vars(expr[[3]]) %in% species)) {
            return(fail("expoente depende de uma especie"))
        }
        n <- tryCatch(eval(expr[[3]], envir = env), error = function(e) NULL)
        if (!is.numeric(n) || length(n) != 1L || !is.finite(n)) {
            return(fail("expoente nao avaliavel"))
        }
        return(sbml_flatten_product(expr[[2]], species, env, power * n))
    }
    if (op == "-" && length(expr) == 2L) {
        inner <- sbml_flatten_product(expr[[2]], species, env, power)
        if (inner$ok) inner$const <- -inner$const
        return(inner)
    }
    if (op %in% c("+", "-")) {
        return(fail("soma ou diferenca de termos que dependem de especies (ex.: Michaelis-Menten, reversivel f - g)"))
    }
    fail(paste0("funcao '", op, "' aplicada a uma especie"))
}

# Compara os expoentes encontrados com a estequiometria dos reagentes.
sbml_match_stoichiometry <- function(powers, stoich, species_names) {
    expected <- numeric()
    for (s in names(stoich)) expected[s] <- stoich[[s]]
    names_all <- union(names(expected), names(powers))
    problems <- character()
    for (s in names_all) {
        e <- if (s %in% names(expected)) expected[[s]] else 0
        g <- if (s %in% names(powers)) powers[[s]] else 0
        if (abs(e - g) > 1e-9) {
            problems <- c(problems, if (e > 0 && g == 0) {
                sprintf("falta o reagente %s (esperado expoente %g, a lei nao o usa)", s, e)
            } else if (e == 0) {
                sprintf("a especie %s aparece na lei mas nao e reagente", s)
            } else {
                sprintf("expoente de %s e %g, mas a estequiometria e %g", s, g, e)
            })
        }
    }
    problems
}

# Estequiometria nomeada (soma se a especie aparece em mais de uma referencia).
sbml_stoich_vector <- function(references) {
    out <- numeric()
    for (ref in references) {
        s <- ref$species
        v <- if (is.null(ref$stoich)) 1 else as.numeric(ref$stoich)
        out[s] <- (if (s %in% names(out)) out[[s]] else 0) + v
    }
    out
}

#' Classify one SBML kinetic law as mass action or general.
#'
#' @param formula Kinetic-law formula (libSBML L3 string).
#' @param reactants,products Lists of references (`species`, `stoich`).
#' @param constants Named list of numeric values for every non-species symbol
#'   (global parameters, compartment sizes and local parameters already merged).
#' @param species Character vector with all species ids.
#' @return A list with `class` ("mass_action", "reversible_mass_action" or
#'   "general"), `k` (numeric, NA when general), `k_forward`/`k_backward`
#'   (reversible only) and `reason` (empty string when supported).
#' @export
sbml_classify_kinetic_law <- function(formula, reactants, products = list(),
                                      constants = list(), species = character()) {
    general <- function(reason) {
        list(class = "general", k = NA_real_, k_forward = NA_real_,
             k_backward = NA_real_, reason = reason)
    }
    if (is.null(formula) || !nzchar(formula)) {
        return(general("reacao sem lei cinetica"))
    }
    expr <- tryCatch(parse(text = formula)[[1]], error = function(e) NULL)
    if (is.null(expr)) {
        return(general("formula nao pode ser interpretada"))
    }
    env <- sbml_math_env(constants)

    as_mass_action <- function(e, refs) {
        flat <- sbml_flatten_product(e, species, env)
        if (!flat$ok) return(list(ok = FALSE, reason = flat$reason))
        problems <- sbml_match_stoichiometry(flat$powers, sbml_stoich_vector(refs), species)
        if (length(problems)) {
            return(list(ok = FALSE, reason = paste(problems, collapse = "; ")))
        }
        if (!is.finite(flat$const) || flat$const == 0) {
            return(list(ok = FALSE, reason = "constante de velocidade invalida (zero ou nao finita)"))
        }
        # `k` fica com sinal: em R, `-k2 * B` e `(-k2) * B`, entao o sinal
        # do termo pode estar dentro da constante.
        list(ok = TRUE, k = flat$const)
    }

    direct <- as_mass_action(expr, reactants)
    if (direct$ok && direct$k < 0) {
        direct <- list(ok = FALSE, reason = "constante de velocidade negativa")
    }
    if (direct$ok) {
        return(list(class = "mass_action", k = direct$k, k_forward = NA_real_,
                    k_backward = NA_real_, reason = ""))
    }

    # f - g: tenta como par direto/inverso (so faz sentido com produtos).
    terms <- sbml_additive_terms(expr)
    if (length(terms) == 2L && length(products) > 0L) {
        # sinal efetivo do termo = sinal do termo * sinal da constante
        eff <- lapply(terms, function(tm) {
            flat <- sbml_flatten_product(tm$expr, species, env)
            sg <- if (flat$ok && is.finite(flat$const) && flat$const < 0) -1 else 1
            list(sign = tm$sign * sg, expr = tm$expr)
        })
        signs <- vapply(eff, `[[`, numeric(1), "sign")
        if (identical(sort(signs), c(-1, 1))) {
            pos <- eff[[which(signs > 0)]]$expr
            neg <- eff[[which(signs < 0)]]$expr
            fwd <- as_mass_action(pos, reactants)
            bwd <- as_mass_action(neg, products)
            if (fwd$ok && bwd$ok) {
                return(list(class = "reversible_mass_action", k = NA_real_,
                            k_forward = abs(fwd$k), k_backward = abs(bwd$k),
                            reason = ""))
            }
            return(general(paste0(
                "forma f - g, mas ",
                if (!fwd$ok) paste0("termo direto: ", fwd$reason) else "",
                if (!fwd$ok && !bwd$ok) "; " else "",
                if (!bwd$ok) paste0("termo inverso: ", bwd$reason) else ""
            )))
        }
    }
    general(direct$reason)
}

# Termos de uma soma/diferenca de nivel superior, com sinal.
sbml_additive_terms <- function(expr, sign = 1) {
    if (is.call(expr)) {
        op <- as.character(expr[[1]])
        if (op == "(") return(sbml_additive_terms(expr[[2]], sign))
        if (op == "+" && length(expr) == 3L) {
            return(c(sbml_additive_terms(expr[[2]], sign),
                     sbml_additive_terms(expr[[3]], sign)))
        }
        if (op == "-" && length(expr) == 3L) {
            return(c(sbml_additive_terms(expr[[2]], sign),
                     sbml_additive_terms(expr[[3]], -sign)))
        }
        if (op == "-" && length(expr) == 2L) {
            return(sbml_additive_terms(expr[[2]], -sign))
        }
    }
    list(list(sign = sign, expr = expr))
}

#' Classify every kinetic law of a model returned by `read_sbml_to_dnar()`.
#'
#' @param parsed Bridge output (`sbml_run_bridge()`), with `reactions`.
#' @param symbols Resolved symbol table (`model$symbols`).
#' @param species Character vector of species ids.
#' @return List with one classification per reaction (see
#'   [sbml_classify_kinetic_law()]) plus `id` and `reversible`.
#' @keywords internal
sbml_classify_model <- function(parsed, symbols, species) {
    lapply(parsed$reactions, function(reaction) {
        values <- sbml_merge_values(symbols, reaction$parameters)
        constants <- values[setdiff(names(values), species)]
        info <- sbml_classify_kinetic_law(
            reaction$formula, reaction$reactants, reaction$products,
            constants, species
        )
        info$id <- reaction$id
        info$reversible <- isTRUE(reaction$reversible)
        info$n_reactants <- sum(vapply(reaction$reactants, function(r) {
            if (is.null(r$stoich)) 1 else as.numeric(r$stoich)
        }, numeric(1)))
        info
    })
}

#' Support table: which reactions can go to the CRN and to DNA.
#'
#' `crn_ok` is TRUE for every reaction (the CRN engine accepts a full-propensity
#' function). `dna_ok` is TRUE only when the law is mass action (or a reversible
#' pair of mass-action laws) and has at most two reactants, which is what the
#' 4-domain translation accepts.
#'
#' @param model Object returned by [read_sbml_to_dnar()].
#' @return A data.frame, one row per reaction.
#' @export
sbml_kinetic_report <- function(model) {
    info <- model$kinetic_info
    if (is.null(info) || length(info) == 0L) {
        return(data.frame())
    }
    data.frame(
        id = vapply(info, `[[`, character(1), "id"),
        reaction = model$reactions,
        kinetic_class = vapply(info, `[[`, character(1), "class"),
        k = vapply(info, function(i) {
            if (i$class == "mass_action") i$k else NA_real_
        }, numeric(1)),
        k_forward = vapply(info, `[[`, numeric(1), "k_forward"),
        k_backward = vapply(info, `[[`, numeric(1), "k_backward"),
        crn_ok = TRUE,
        dna_ok = vapply(info, function(i) {
            i$class != "general" && i$n_reactants <= 2
        }, logical(1)),
        reason = vapply(info, function(i) {
            if (i$class != "general" && i$n_reactants > 2) {
                "mais de 2 reagentes: a traducao de 4 dominios nao suporta"
            } else {
                i$reason
            }
        }, character(1)),
        stringsAsFactors = FALSE
    )
}

#' Split reversible mass-action reactions into forward and backward reactions.
#'
#' The default model keeps a reversible reaction as one net-flux reaction
#' (`A -> B`). The DNA translation needs elementary steps, so this returns a
#' copy of `model` where each `reversible_mass_action` reaction `A -> B` becomes
#' `A -> B` (ki = k_forward) followed by `B -> A` (ki = k_backward).
#' Other reactions are kept untouched. `ki` of the result is numeric when every
#' reaction ends up mass action, otherwise a list.
#'
#' @param model Object returned by [read_sbml_to_dnar()].
#' @return A model of the same class with `reactions`, `ki`, `kinetic_info`
#'   and `kinetic_laws` expanded; `split_from` maps new reactions to originals.
#' @export
sbml_split_reversible <- function(model) {
    info <- model$kinetic_info
    if (is.null(info)) {
        stop("modelo sem kinetic_info: gere-o com read_sbml_to_dnar()")
    }
    ki_in <- as.list(model$ki)
    reactions <- character()
    ki <- list()
    new_info <- list()
    laws <- list()
    split_from <- integer()
    for (i in seq_along(info)) {
        inf <- info[[i]]
        if (inf$class == "reversible_mass_action") {
            sides <- strsplit(model$reactions[[i]], "->", fixed = TRUE)[[1]]
            left <- trimws(sides[[1]])
            right <- trimws(sides[[2]])
            reactions <- c(reactions, paste(left, "->", right),
                           paste(right, "->", left))
            ki <- c(ki, list(inf$k_forward, inf$k_backward))
            fwd <- inf; fwd$class <- "mass_action"; fwd$k <- inf$k_forward
            bwd <- inf; bwd$class <- "mass_action"; bwd$k <- inf$k_backward
            fwd$id <- paste0(inf$id, "_fwd"); bwd$id <- paste0(inf$id, "_rev")
            fwd$n_reactants <- bwd$n_reactants <- NA_real_
            new_info <- c(new_info, list(fwd, bwd))
            laws <- c(laws, list(model$kinetic_laws[[i]], model$kinetic_laws[[i]]))
            split_from <- c(split_from, i, i)
        } else {
            reactions <- c(reactions, model$reactions[[i]])
            ki <- c(ki, list(ki_in[[i]]))
            new_info <- c(new_info, list(inf))
            laws <- c(laws, list(model$kinetic_laws[[i]]))
            split_from <- c(split_from, i)
        }
    }
    # n_reactants do passo inverso = soma dos coeficientes do lado direito
    for (j in seq_along(new_info)) {
        if (is.na(new_info[[j]]$n_reactants)) {
            lhs <- trimws(strsplit(reactions[[j]], "->", fixed = TRUE)[[1]][[1]])
            new_info[[j]]$n_reactants <- sbml_count_reactants(lhs)
        }
    }
    if (all(vapply(ki, function(x) is.numeric(x) && length(x) == 1L, logical(1)))) {
        ki <- unlist(ki, use.names = FALSE)
    }
    model$reactions <- reactions
    model$ki <- ki
    model$kinetic_info <- new_info
    model$kinetic_laws <- laws
    model$split_from <- split_from
    model
}

# Numero de moleculas do lado esquerdo "2 A + B" -> 3 ("0" -> 0).
sbml_count_reactants <- function(side) {
    side <- trimws(side)
    if (side %in% c("", "0")) return(0)
    parts <- trimws(strsplit(side, "+", fixed = TRUE)[[1]])
    sum(vapply(parts, function(p) {
        m <- regmatches(p, regexpr("^[0-9.]+", p))
        if (length(m) && nzchar(m)) as.numeric(m) else 1
    }, numeric(1)))
}

# Funcao de propensao completa para uma lei que nao e acao de massas.
sbml_general_rate_function <- function(formula, constants) {
    expr <- parse(text = formula)[[1]]
    fun <- function(t, y) {
        env <- sbml_math_env(c(constants, as.list(y)))
        env$time <- t
        val <- suppressWarnings(eval(expr, envir = env))
        if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
            return(0)
        }
        as.numeric(val)
    }
    attr(fun, "sbml_full_propensity") <- TRUE
    fun
}
