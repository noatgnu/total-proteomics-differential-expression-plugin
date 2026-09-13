#!/usr/bin/env Rscript
# Resolves the R dependency tree and system requirements for the current OS into deps-partial.json.

suppressMessages({
  if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
})

pkgs <- character(0)
if (file.exists("r-packages.txt")) {
  pkgs <- readLines("r-packages.txt")
  pkgs <- trimws(pkgs)
  pkgs <- pkgs[nzchar(pkgs)]
}

packages_out <- list()
sysreqs_out <- list()

if (length(pkgs) > 0) {
  refs <- paste0("any::", pkgs)

  deps <- tryCatch(pak::pkg_deps(refs, dependencies = NA), error = function(e) {
    message("pak::pkg_deps failed: ", conditionMessage(e))
    NULL
  })
  if (!is.null(deps) && nrow(deps) > 0) {
    packages_out <- lapply(seq_len(nrow(deps)), function(i) {
      list(name = deps$package[i], version = as.character(deps$version[i]), repository = as.character(deps$type[i]))
    })
  }

  sysreqs <- tryCatch(pak::pkg_sysreqs(refs), error = function(e) {
    message("pak::pkg_sysreqs failed: ", conditionMessage(e))
    NULL
  })
  message("Raw sysreqs result (for debugging if the parsed list below looks wrong):")
  print(sysreqs)

  if (!is.null(sysreqs) && is.data.frame(sysreqs$packages) && nrow(sysreqs$packages) > 0) {
    sp <- sysreqs$packages
    sysreqs_out <- tryCatch(
      lapply(seq_len(nrow(sp)), function(i) {
        pkg_names <- sp$packages[[i]]
        list(name = as.character(sp$sysreq[i]), packages = as.character(unlist(pkg_names)))
      }),
      error = function(e) {
        message("Structured sysreqs parse failed (", conditionMessage(e), "); falling back to pre_install commands")
        list(list(name = "unparsed", packages = as.character(sysreqs$pre_install)))
      }
    )
  }
}

out <- list(
  generatedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  rVersion = paste(R.version$major, R.version$minor, sep = "."),
  packages = packages_out,
  systemRequirements = sysreqs_out
)

jsonlite::write_json(out, "deps-partial.json", auto_unbox = TRUE, pretty = TRUE, null = "null")
