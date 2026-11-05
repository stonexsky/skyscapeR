#' Named elevation palettes for DEM terrain panoramas
#' @description Choose land and water colours for a site panorama. The presets
#' use illustrative elevation tints for different landscapes.
#' @details Pass a preset name to the terrain.palette argument of \code{\link{createHor}}.
#' Use \code{demPalettes("standalone")} to inspect the colours in that preset.
#' Temperate uses green and ochre, arid uses sand and red-brown, alpine runs
#' from forest greens to pale rock and white, and coastal pairs sand and grass
#' with a deep-blue to turquoise water ramp. Relief gives land one neutral tint
#' for inspecting slope shading. Each preset has its own sea_cols vector.
#' Land stops are ordered low to high and water stops deep to shallow.
#' \code{\link{createHor}} treats elevations at or below sea_level as water for colouring.
#' An inland depression below that threshold also receives the water tint.
#' sea_level=NULL uses the land ramp for all elevations.
#' @param name Optional exact preset name. Omit to retrieve the full catalogue.
#' @return With a name, a list containing land_cols and sea_cols. Without a
#' name, a named list of those pairs.
#' @seealso \code{\link{createHor}}
#' @examples
#' names(demPalettes())
#' demPalettes("coastal")$sea_cols
#' @export
demPalettes <- function(name=NULL) {
  presets <- list(
    standalone=list(land_cols=c("#2e5d3a","#6b8f5a","#b6b08a","#9a7b5c",
                                "#6d625b","#cfcfcf","#ffffff"),
                    sea_cols=c("#082567","#0b3d91","#2b5aa5","#4fa3d1")),
    temperate=list(land_cols=c("#395F36","#557C42","#739550","#8EA663","#ADAD78"),
                   sea_cols=c("#153D54","#40758A","#84B9BB")),
    arid=list(land_cols=c("#A66A47","#CE9469","#E1B989","#AF8569","#DBC9AC"),
              sea_cols=c("#15495A","#327C86","#A0D7CC")),
    alpine=list(land_cols=c("#324E3E","#6B855D","#AAA48D","#898781","#D0D1D0","#FAFAF7"),
                sea_cols=c("#183952","#396D86","#91CBD2")),
    coastal=list(land_cols=c("#CFBF95","#95A578","#637D59","#A6A494","#DCD3B8"),
                 sea_cols=c("#0B2948","#15536C","#31899A","#91D4BD")),
    relief=list(land_cols="#C6C2B6",sea_cols=c("#425E6B","#8CAEB3"))
  )
  if (is.null(name)) return(presets)
  if (!is.character(name) || length(name)!=1L || is.na(name) || !name %in% names(presets))
    stop("Choose a terrain palette from ",paste(names(presets),collapse=", "),".")
  presets[[name]]
}
