#!/usr/bin/env Rscript
#
# Generate pkgdown article vignettes from the repo-root markdown docs, so the
# compatibility and architecture references render inside the gsemr docs site
# without duplicating their content. The repo-root .md files stay the single
# source of truth; the generated vignettes/*.Rmd are gitignored.
#
# Run from bindings/r/:  Rscript tools/build-vignettes.R

pkg_root  <- normalizePath(file.path(dirname(sub("--file=", "",
               grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
repo_root <- normalizePath(file.path(pkg_root, "..", ".."))
vig_dir   <- file.path(pkg_root, "vignettes")
dir.create(vig_dir, showWarnings = FALSE)

# (source .md at repo root, output vignette slug, article title)
docs <- list(
  list(src = "API_COMPAT.md",  out = "compatibility", title = "Compatibility with R GenomicSEM"),
  list(src = "ARCHITECTURE.md", out = "architecture",  title = "Architecture & Differences from R GenomicSEM")
)

for (d in docs) {
  src_path <- file.path(repo_root, d$src)
  if (!file.exists(src_path)) {
    message("skip: ", d$src, " not found at ", src_path)
    next
  }
  body <- readLines(src_path, warn = FALSE)

  # Drop the document's own leading H1 — the vignette title supplies it, and a
  # duplicate H1 reads oddly in the rendered article.
  first_h1 <- which(grepl("^# ", body))[1]
  if (!is.na(first_h1)) body <- body[-seq_len(first_h1)]

  # Rewrite cross-doc links so they resolve within the rendered site:
  #   ./ARCHITECTURE.md(#anchor) -> architecture.html(#anchor)
  #   ./API_COMPAT.md(#anchor)   -> compatibility.html(#anchor)
  # Other ./*.md links (README, CONTRIBUTING) point back at the GitHub repo.
  repo_url <- "https://github.com/PoHsuanLai/gsem/blob/master"
  body <- gsub("\\./ARCHITECTURE\\.md", "architecture.html", body)
  body <- gsub("architecture\\.html#", "architecture.html#", body)
  body <- gsub("\\./API_COMPAT\\.md", "compatibility.html", body)
  body <- gsub("\\./(README|CONTRIBUTING|CHANGELOG|TODO)\\.md",
               paste0(repo_url, "/\\1.md"), body)
  # Images live under bench/ at the repo root; point them at the raw GitHub URL.
  raw_url <- "https://raw.githubusercontent.com/PoHsuanLai/gsem/master"
  body <- gsub("\\]\\((bench/[^)]+)\\)", paste0("](", raw_url, "/\\1)"), body)

  header <- c(
    "---",
    paste0("title: \"", d$title, "\""),
    "output: rmarkdown::html_vignette",
    "vignette: >",
    paste0("  %\\VignetteIndexEntry{", d$title, "}"),
    "  %\\VignetteEngine{knitr::rmarkdown}",
    "  %\\VignetteEncoding{UTF-8}",
    "---",
    ""
  )

  # "Edit this page" footer pointing at the TRUE source (the repo-root .md),
  # not the generated vignette (which is gitignored).
  edit_url <- paste0(repo_url, "/", d$src)
  footer <- c(
    "",
    "----",
    "",
    paste0("*This page is generated from [`", d$src, "`](", edit_url,
           ") in the repository — [edit it there](", edit_url, ").*")
  )

  out_path <- file.path(vig_dir, paste0(d$out, ".Rmd"))
  writeLines(c(header, body, footer), out_path)
  message("wrote ", out_path)
}
