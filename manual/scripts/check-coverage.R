# Run from the Quarto project root with the documented package installed.
# Checks documentation coverage. It does not run the vignette or its studies.
ns <- loadNamespace("skyscapeR")
version <- as.character(getNamespaceVersion(ns))
archive <- file.path("archives", paste0("skyscapeR_", version, ".tar.gz"))
stopifnot(file.exists(archive))
w <- tempfile()
dir.create(w)
utils::untar(archive, files=c("skyscapeR/NAMESPACE", "skyscapeR/vignettes/skyscapeR.Rmd"), exdir=w)
nsp <- readLines(file.path(w,"skyscapeR","NAMESPACE"))
exports <- sub("^export\\((.*)\\)$","\\1",grep("^export",nsp,value=TRUE))
methods <- sub("^S3method\\(([^,]+),(.+)\\)$","\\1.\\2",grep("^S3method",nsp,value=TRUE))
read_text <- function(file) paste(readLines(file,warn=FALSE),collapse="\n")
entry_files <- Sys.glob("fn-*.qmd")
entries <- vapply(entry_files, read_text, character(1))
titles <- vapply(entry_files, function(path) readLines(path, n=1L), character(1))
names(entries) <- sub("^# `", "", sub("[(][)].*$", "", titles))
api <- c(exports, methods)
stopifnot(setequal(names(entries),api), length(entries)==length(api))
formal_count <- 0L
for (name in api) {
  args <- names(formals(get(name,ns)))
  formal_count <- formal_count+length(args)
  entry <- entries[[name]]
  rows <- regmatches(entry,gregexpr('<th scope="row"><code>[^<]+</code></th>',entry))[[1]]
  found <- sub('</code></th>','',sub('<th scope="row"><code>','',rows,fixed=TRUE),fixed=TRUE)
  if (!identical(found,args)) stop("Argument coverage changed for ",name)
}
vignette <- readLines(file.path(w,"skyscapeR","vignettes","skyscapeR.Rmd"))
guide <- read_text("vignette-guide.qmd")
titles <- character()
sections <- list()
current <- "Setup"
in_chunk <- FALSE
chunk <- character()
calls <- function(e) {
  if (!is.call(e)) return(character())
  head <- if(is.symbol(e[[1]])) as.character(e[[1]]) else character()
  unique(c(head,unlist(lapply(as.list(e)[-1],calls))))
}
chunk_count <- 0L
for (line in vignette) {
  if (in_chunk) {
    if (!grepl("^```",line)) {chunk <- c(chunk,line); next}
    parsed <- parse(text=chunk)
    explicit <- intersect(unique(unlist(lapply(parsed,calls))),api)
    sections[[current]] <- union(sections[[current]],explicit)
    in_chunk <- FALSE
    chunk_count <- chunk_count+1L
  } else if (grepl("^```\\{r",line)) {
    in_chunk <- TRUE
    chunk <- character()
  } else if (grepl("^#{2,3} ",line)) {
    current <- sub(" \\{#.*$","",sub("^#{2,3} ","",line))
    titles <- c(titles,current)
  }
}
guide_parts <- strsplit(guide,"\n### ",fixed=TRUE)[[1]][-1]
stopifnot(!in_chunk, identical(sub(" \\{#.*","",sub("\n.*","",guide_parts)),titles))
for(i in seq_along(titles)) {
  for(n in sections[[titles[i]]]) {
    link <- paste0("fn-",gsub(".","-",tolower(n),fixed=TRUE),".html)")
    if (!grepl(link,guide_parts[i],fixed=TRUE)) stop("Vignette route missing for ",n," in ",titles[i])
  }
}
unlink(w,recursive=TRUE)
cat(length(api),"entries;",formal_count,"formal arguments;",length(titles),
    "vignette topics;",chunk_count,"R chunks parsed.\n")
