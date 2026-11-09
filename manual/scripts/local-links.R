# Preserve the book title's home link for offline use and GitHub project sites.
for (page in list.files("docs", pattern = "[.]html$", full.names = TRUE)) {
  lines <- readLines(page, encoding = "UTF-8", warn = FALSE)
  revised <- gsub('href="(\\./|/)"', 'href="index.html"', lines)
  revised <- gsub('href="/([a-zA-Z0-9_-]+[.]html)([#"])', 'href="\\1\\2', revised)
  revised <- gsub('href="([a-zA-Z0-9_-]+)[.]qmd([#"])', 'href="\\1.html\\2', revised)
  if (!identical(revised, lines)) writeLines(revised, page, useBytes = TRUE)
}
