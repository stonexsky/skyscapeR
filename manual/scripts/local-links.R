# Preserve the book title's home link for offline use and GitHub project sites.
for (page in list.files("docs", pattern = "[.]html$", full.names = TRUE)) {
  lines <- readLines(page, encoding = "UTF-8", warn = FALSE)
  revised <- gsub('href="(\\./|/)"', 'href="index.html"', lines)
  revised <- gsub('href="/([a-zA-Z0-9_-]+[.]html)([#"])', 'href="\\1\\2', revised)
  revised <- gsub('href="([a-zA-Z0-9_-]+)[.]qmd([#"])', 'href="\\1.html\\2', revised)
  if (!identical(revised, lines)) writeLines(revised, page, useBytes = TRUE)
}

# Resolve published single-page study permalinks when hosting the chapter edition.
home <- "docs/index.html"
marker <- "<!-- study-permalinks -->"
lines <- readLines(home, encoding = "UTF-8", warn = FALSE)
if (!any(grepl(marker, lines, fixed = TRUE))) {
  script <- paste(c(marker, "<script>",
    "const studyPages = {",
    "  '#a02-ephemeris-accuracy': 'a02-ephemeris.html',",
    "  '#h01-field-validation': 'study-register.html',",
    "  '#h02-skyline-comparison': 'study-register.html'",
    "};",
    "function openStudy() {",
    "  const page = studyPages[location.hash];",
    "  if (page) location.replace(page + location.hash);",
    "}",
    "openStudy();",
    "window.addEventListener('hashchange', openStudy);",
    "</script>"), collapse = "\n")
  writeLines(sub("</head>", paste0(script, "\n</head>"), lines, fixed = TRUE), home, useBytes = TRUE)
}
