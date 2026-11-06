utils::globalVariables(c("constellations"))
#' Line segments for drawing Western and Chinese constellations
#' @description Identify the line segments used to represent Western and Chinese constellations in a
#' reconstruction of the sky.
#' @details The segment list supplies a drawing convention for overlaying figures on a sky sketch or planisphere.
#' Within a row, read successive pairs of Hipparcos identifiers as separate lines.
#' Identifiers are endpoints, not a polygon to join sequentially or a boundary on the
#' celestial sphere. Epoch and observer location enter when the stars are plotted.
#'
#' All 88 ordered lists of Hipparcos line endpoints match the Western constellation
#' line file distributed with Stellarium 0.22.2. This comparison identifies the
#' drawing convention used here.
#' Credit: Stellarium's team. The upstream Western sky culture specifies
#' CC BY-SA 4.0 International Public License + Free Art License.
#' The line endpoints are unchanged and repackaged as an R data frame.
#' See the installed COPYRIGHTS file for sources and licence links.
#'
#' The Chinese records preserve all 318 entries in Stellarium 0.22.2's Chinese
#' line file: 312 have segments (including 57 single-star figures), and six have
#' no segments. Empty HIP fields are skipped; repeated endpoints mark a single
#' star with a small ring. The 28 mansion rows are duplicated as a separately
#' selectable set, ordered from Horn through Chariot. These are asterism figures,
#' not mansion sectors. Selecting both sets does not draw duplicate figures.
#'
#' English names and Chinese pinyin are retained from upstream. Pinyin without
#' tones is ambiguous: use CHN identifiers or subset by set and English name.
#' Credit: Karrie Berglund, Sun Shuwei and Stellarium contributors; Chinese data
#' are CC BY-SA 4.0. The tradition follows later imperial catalogues as described
#' by upstream, including Yi Shitong's Chinese and Western Contrast Star Chart
#' and Catalogue 1950.0. It is not a reconstruction for every historical epoch.
#'
#' All Chinese endpoints resolve in the bundled stellar catalogue. The
#' 2026-11-04 supplement uses SIMBAD ICRS/J2000 astrometry, with missing V
#' magnitudes from Hipparcos via VizieR. Five unavailable radial velocities
#' are stored as zero placeholders, not measured zeros; HIP 108797 uses
#' Hipparcos 2007 proper motion and parallax from SIMBAD's measurement tables.
#' Binary/system identifications and long-epoch extrapolation retain the
#' limitations of those source measurements. The installed ephemeris file
#' star-supplement-2026-11-04.csv records per-field references and assumptions.
#' @references Stellarium 0.22.2. Western constellation lines. \url{https://github.com/Stellarium/stellarium/blob/v0.22.2/skycultures/western/constellationship.fab}.
#' Chinese source: \url{https://github.com/Stellarium/stellarium/tree/v0.22.2/skycultures/chinese}.
#' @seealso \code{\link{sketchSky}}, \code{\link{importConstellations}},
#' \code{\link{star.names}}
#' @examples data(constellations)
#' constellations[constellations$name2 == "Ori", c("name1", "HIP")]
#' subset(constellations, set == "chinese_mansions", select = c(name1, name2))
#' system.file("ephemeris", "star-supplement-2026-11-04.csv", package = "skyscapeR")
#' @docType data
#' @usage data(constellations)
#' @format A data frame with 434 rows and six columns.
#' \describe{
#'   \item{set}{Sky-culture identifier. Use "modern" (88 rows), "chinese" (318 rows), or "chinese_mansions" (28 rows).}
#'   \item{name1}{Full constellation name.}
#'   \item{name2}{Western abbreviation or stable Chinese identifier (CHN001 to CHN318).}
#'   \item{name3}{Chinese pinyin; missing for Western records.}
#'   \item{name4}{Repeated CHN identifier for Chinese records; missing for Western records.}
#'   \item{HIP}{Space-separated Hipparcos identifiers. Each successive pair gives the two endpoints of one line segment.}
#' }
#' @keywords datasets
"constellations"
