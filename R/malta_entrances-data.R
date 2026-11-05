#' Entrance bearings of Maltese temples
#' @description Compare the directions and angular widths of entrance corridors at the prehistoric
#' Maltese temples.
#' @details The left and right bearings delimit the field of view from a structure. Preserve their
#' order around north when deriving a centre and width, and specify how that sector is
#' represented as a density. Latitude, longitude and phase support site-specific
#' transformation and grouping. The companion horizon profiles provide the horizon
#' altitudes used to translate bearings into possible declinations.
#' @seealso \code{\link{malta_horizons}}
#' @docType data
#' @usage data(malta_entrances)
#' @format A tibble (also a data frame) with 32 rows and six columns.
#' \describe{
#'   \item{Name}{Structure name.}
#'   \item{Latitude}{Site latitude in degrees north.}
#'   \item{Longitude}{Site longitude in degrees east.}
#'   \item{Az L}{Azimuth of the left edge of the entrance, in degrees.}
#'   \item{Az R}{Azimuth of the right edge of the entrance, in degrees.}
#'   \item{Phase}{Archaeological phase assigned to the structure.}
#' }
#' @keywords datasets
#' @references Lomsdalen, T. (2021). \emph{Viewscapes and Cosmology in the Prehistoric Temples of Malta}. PhD Dissertation - University of Malta.
#' @references Silva, F. and Lomsdalen, T. (2025). ‘No easy way from the earth to the stars’: a new statistical approach to the orientation of the Maltese Temples. \emph{Archaeological and Anthropological Sciences}  \doi{10.1007/s12520-025-02208-4}
#' @examples
#' data(malta_entrances)
#' mid_az <- (malta_entrances$`Az L` + malta_entrances$`Az R`)/2
#' plotAzimuth(mid_az)
"malta_entrances"
