#' Import user-defined constellation lines
#' @description Read custom line connections for use in a sky sketch.
#' @param file Path to a Stellarium legacy constellationship.fab text file.
#' @param set Label for this collection, default "custom".
#' @return A data frame with the same columns as \code{\link{constellations}}, accepted by
#' the \code{display$constellations} setting of \code{\link{sketchSky}}. Subset rows to select individual figures.
#' @details Each non-comment line contains a unique figure identifier, the number
#' of segments, then two Hipparcos numbers per segment. Endpoints describe
#' separate pairs, not one continuous polygon. Blank lines and comments beginning
#' with a hash are ignored. Non-negative segment counts and positive integer HIP identifiers are validated.
#' Zero-segment records are retained with an empty HIP field and skipped when drawn.
#' Repeated endpoints represent a single-star figure, drawn as a small ring.
#' This supports Stellarium's legacy constellationship.fab format, not the JSON
#' sky-culture format introduced in Stellarium 25.1. Only lines are imported,
#' not artwork, labels, boundaries or star positions. Endpoint stars must be in
#' skyscapeR's stellar catalogue; unavailable stars produce warnings when drawn.
#' The import does not modify the bundled catalogue or register global names.
#' Save the returned object with saveRDS(), or load the file in each script.
#' Respect the source sky culture's attribution and redistribution terms.
#' @references Stellarium User Guide, Adding Sky Cultures (legacy editions) and
#' Sky Cultures (current editions): \url{https://stellarium.org/guide/}.
#' Legacy example: \url{https://github.com/Stellarium/stellarium/blob/v0.22.2/skycultures/western/constellationship.fab}.
#' @seealso \code{\link{sketchSky}}, \code{\link{constellations}}
#' @examples
#' file <- tempfile(fileext=".fab")
#' writeLines("MyFigure 2 27989 25336 25336 25930", file)
#' figures <- importConstellations(file)
#' unlink(file)
#' figures[, c("name1", "HIP")]
#' \dontrun{
#' sky <- sketchSky(malta_horizons[[1]], time="2026/12/21 20:00",
#'   timezone="UTC",
#'   display = list(constellations = figures))
#' plot(malta_horizons[[1]], sky=sky)
#' }
#' @export
importConstellations <- function(file, set="custom") {
  if (!is.character(file) || length(file)!=1L || is.na(file) || !file.exists(file))
    stop("Supply an existing constellation line file.",call.=FALSE)
  if (!is.character(set) || length(set)!=1L || is.na(set) || !nzchar(trimws(set)))
    stop("Supply one non-empty collection label.",call.=FALSE)
  text <- readLines(file,warn=FALSE,encoding="UTF-8")
  text <- sub("^\ufeff","",text)
  text <- trimws(sub("#.*$","",text))
  rows <- which(nzchar(text))
  if (!length(rows)) stop("No constellation lines found.",call.=FALSE)
  tokens <- strsplit(text[rows],"\\s+")
  for (i in seq_along(tokens)) {
    x <- tokens[[i]]
    valid <- length(x)>=2L && all(grepl("^[0-9]+$",x[-1L]))
    nums <- suppressWarnings(as.numeric(x[-1L]))
    valid <- valid && all(is.finite(nums)) && nums[1L]>=0 && all(nums[-1L]>0) &&
      all(nums<=.Machine$integer.max) && length(x)==2+2*nums[1L]
    if (!valid) stop("Invalid constellation record on line ",rows[i],
      ": expected ID, non-negative segment count, and two positive HIP numbers per segment.",call.=FALSE)
  }
  ids <- vapply(tokens,`[`,character(1),1L)
  if (anyDuplicated(ids)) stop("Constellation identifiers must be unique within a file.",call.=FALSE)
  data.frame(set=set,name1=ids,name2=ids,name3=ids,name4=ids,
    HIP=vapply(tokens,function(x) paste(as.integer(x[-c(1L,2L)]),collapse=" "),character(1)),
    stringsAsFactors=FALSE)
}
