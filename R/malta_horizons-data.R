#' DEM-derived horizon profiles for Maltese temples
#' @description Account for the terrain visible from each Maltese temple when interpreting the
#' orientation of its entrance.
#' @details Tore Lomsdalen generated the profiles for his PhD research using Andrew Smith’s Horizon software and a 5 m DEM made available through Heritage Malta. The profiles describe how horizon altitude and its uncertainty vary with azimuth at each
#' site. The list is unnamed; inspect each profile's metadata before pairing it with an
#' entrance observation, especially after subsetting or reordering. Use the same matched
#' profile when transforming a bearing distribution and when simulating its
#' random-orientation null.
#' @seealso \code{\link{malta_entrances}}
#' @docType data
#' @usage data(malta_horizons)
#' @format An unnamed list of 32 skyscapeR.horizon objects. Access a profile by position, for example malta_horizons[[1]].
#' \describe{
#'   \item{metadata}{Site name and georeference (latitude, longitude, elevation).}
#'   \item{data}{Horizon samples with columns az, alt and alt.unc. All three angles are in degrees.}
#' }
#' @keywords datasets
#' @references Lomsdalen, T. (2021). \emph{Viewscapes and Cosmology in the Prehistoric Temples of Malta}. PhD Dissertation - University of Malta.
#' @references Silva, F. and Lomsdalen, T. (2025). ‘No easy way from the earth to the stars’: a new statistical approach to the orientation of the Maltese Temples. \emph{Archaeological and Anthropological Sciences} \doi{10.1007/s12520-025-02208-4}
#' @examples
#' data(malta_horizons)
#' plot(malta_horizons[[1]])
#' plot(malta_horizons[[10]])
"malta_horizons"