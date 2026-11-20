# Resolve DOI macros at build time so the offline manual has working links.
body_html <- function(source) capture.output(tools::Rd2HTML(source,
  standalone=FALSE, no_links=TRUE, texmath="text", stages=c("build","install","render")))
write_help <- function(source, output) writeLines(c("```{=html}",body_html(source),"```"),output)
write_help("source/star.phases.Rd","_help.qmd")
write_help("source/createHor.Rd","_horizon-help.qmd")
for (n in tools::file_path_sans_ext(list.files("source",pattern="[.]Rd$"))) {
  slug <- tolower(gsub("[.]","-",n))
  body <- body_html(file.path("source",paste0(n,".Rd")))
  base <- "https://stonexsky.github.io/skyscapeR/manual/"
  if (file.exists(file.path("source",paste0(slug,"-help.qmd"))))
    writeLines(c("```{=html}",gsub(base,"",body,fixed=TRUE),"```"),file.path("source",paste0(slug,"-help.qmd")))
  path <- file.path("help",paste0(n,".html"))
  header <- sub("(</nav>).*","\\1",readLines(path,n=1))
  body <- gsub("../help/figures/","../figures/",body,fixed=TRUE)
  body <- sub("<h3>Scientific checks</h3>",'<h3 id="scientific-checks">Scientific checks</h3>',body,fixed=TRUE)
  writeLines(c(header,gsub(base,"../",body,fixed=TRUE),"</body></html>"),path)
}
