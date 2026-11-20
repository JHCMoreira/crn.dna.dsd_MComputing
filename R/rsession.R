# cat("===== R VERSION =====\n")
# sessionInfo()

# cat("\n===== INSTALLED R PACKAGES =====\n")
# pkgs <- c(
#   "JuliaCall",
#   "diffeqr",
#   "pak",
#   "devtools",
#   "remotes"
# )

# for (p in pkgs) {
#   cat("\n", p, ": ")
#   if (requireNamespace(p, quietly = TRUE)) {
#     cat(as.character(packageVersion(p)), "\n")
#   } else {
#     cat("NOT INSTALLED\n")
#   }
# }

if (!isNamespace(topenv())) {
    library(JuliaCall)

    JuliaCall:::julia_version()
}

