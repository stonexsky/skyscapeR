#' Catalogue of star names and designations
#' @description Relate stellar proper names and other celestial object names to the catalogue identifiers needed to calculate
#' stellar positions.
#' @details Use this lookup to resolve a named target before comparing its position across epochs.
#' Each entry associates a common name with the identifier used by the fixed-star routines.
#' Positional and photometric calculations use the stellar catalogue; this table describes
#' the naming convention and does not establish that a historical observer recognised the
#' same figure.
#'
#' The stellar proper names come from the IAU Catalog of Star Names, compiled by
#' the International Astronomical Union Working Group on Star Names. The table also
#' includes lookup entries for supernovae or their remnants, deep-sky objects,
#' the Galactic Pole and the Galactic Centre.
#' @references International Astronomical Union, Working Group on Star Names. IAU Catalog of Star Names. \url{https://exopla.net/star-names/modern-iau-star-names/}.
#' @seealso \code{\link{constellations}}
#' @examples data(star.names)
#' star.names$identifier[star.names$Western == "Aldebaran"]
#' @docType data
#' @usage data(star.names)
#' @format A data frame with 304 rows and two character columns.
#' \describe{
#'   \item{identifier}{Catalogue identifier used to look up the object in the ephemeris.}
#'   \item{Western}{Object name used in the package.}
#' }
#' @keywords datasets
"star.names"
