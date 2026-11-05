#' Reconnect saved horizons to moved panorama PNGs
#' @description Update the external PNG paths in one horizon or a list after moving a project.
#' @param hor A skyscapeR.horizon object or a list of horizon objects.
#' @param path A directory containing the moved PNGs, or a character vector of PNG filenames in any order.
#' @param recursive Search subdirectories when path names a directory. Defaults to FALSE.
#' @details A saved horizon records the checksum of its displayed PNG. Copy that PNG with the RDS, including the levelled copy for a calibrated photograph. Matching checks file contents, so a renamed file can be recovered and different panoramas with the same filename can be distinguished. For identical copies, the original basename is preferred, then the first matching candidate. A skyline-only horizon passes through unchanged. If any required image cannot be verified, the function stops before returning a result. Assign the returned object and use saveRDS() to keep the updated paths. The original PeakFinder archive used for native-layout export is separate.
#' @return The horizon or list, retaining its names, order and metadata with updated absolute PNG paths.
#' @examples
#' \dontrun{
#' sites <- readRDS("sites.rds")
#' sites <- relinkHor(sites, "panoramas", recursive = TRUE)
#' saveRDS(sites, "sites.rds")
#' # A single renamed image
#' hh <- relinkHor(hh, "photos/west-ridge.png")
#' }
#' @export
relinkHor <- function(hor, path, recursive=FALSE) {
  single <- inherits(hor, "skyscapeR.horizon")
  sites <- if (single) list(hor) else hor
  if (!is.list(sites) || !all(vapply(sites, inherits, logical(1), "skyscapeR.horizon")))
    stop("Supply a horizon or a list of horizons.")
  if (!is.logical(recursive) || length(recursive)!=1L || is.na(recursive))
    stop("recursive must be TRUE or FALSE.")
  if (!is.character(path) || !length(path) || anyNA(path) || any(!nzchar(path)))
    stop("Supply a directory or one or more PNG paths.")
  if (length(path)==1L && dir.exists(path)) {
    path <- list.files(path, pattern="[.]png$", full.names=TRUE,
                       recursive=recursive, ignore.case=TRUE)
    path <- path[!dir.exists(path)]
  } else if (any(!file.exists(path) | dir.exists(path))) {
    stop("A supplied image path is missing or is a directory.")
  }
  files <- unique(normalizePath(path, mustWork=TRUE))
  hashes <- rep(NA_character_, length(files))
  for (i in seq_along(sites)) {
    im <- sites[[i]]$metadata$image
    if (is.null(im)) next
    if (!is.character(im$md5) || length(im$md5)!=1L || is.na(im$md5))
      stop("Horizon ",i," has no saved PNG checksum.")
    found <- FALSE
    for (j in order(basename(files)!=basename(im$path))) {
      if (is.na(hashes[j])) hashes[j] <- unname(tools::md5sum(files[j]))
      if (identical(hashes[j], unname(im$md5))) {
        sites[[i]]$metadata$image$path <- files[j]
        found <- TRUE
        break
      }
    }
    if (!found) stop("Cannot find a checksum match for horizon ",i," (",im$path,").")
  }
  if (single) sites[[1L]] else sites
}
