

#' Find celestial targets within declination and time ranges
#' @description Screen solar, lunar and stellar candidates whose declinations overlap a measured
#' orientation range during a proposed period of use.
#' @details The search tests overlap in declination. Solar and lunar limiting declinations are
#' sampled across the years of the interval. Site coordinates enable zenith targets
#' and a horizon enables the spatial equinox. The named equinoxes and solstices use
#' geocentric apparent solar longitudes of 0, 90, 180 and 270 degrees. Dates are UT
#' civil days in the requested calendar. Different dates for the first and last year
#' appear in endpoint order separated by a slash. Intermediate years can differ.
#' The sunrise/set date ranges also use the two endpoint years.
#' For each star passing the catalogue magnitude cut, the search computes its
#' declination at noon UT on 1 January Gregorian of every astronomical year,
#' including both endpoints. It tests overlap before rounding and linearly
#' interpolates between annual samples to estimate entry and exit. Separate
#' passages are reported separately, with boundaries rounded to the nearest year.
#' Shorter variations within a sampling year are unresolved. Magnitude describes
#' catalogue brightness and does not model atmospheric visibility.
#'
#' Silva and Lomsdalen (2025, table 5) give an example of stellar target screening for the Maltese temples.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @references Silva, F. and Lomsdalen, T. (2025). 'No easy way from the earth to the stars': a new statistical approach to the orientation of the Maltese temples. Archaeological and Anthropological Sciences 17, 96. \doi{10.1007/s12520-025-02208-4}.
#' @seealso \code{\link{star.phases}}, \code{\link{riseset}}, \code{\link{sketchSky}}, \code{\link{spatial.equinox}}, \code{\link{solar.date}}, \code{\link{EFM}}, \code{\link{planet.extremes}}
#' @return A list of candidate solar, lunar and stellar matches. The tables give target names and declination limits in degrees. Stellar matches also carry catalogue identifiers and visual magnitudes. Set dates to TRUE to request stellar date ranges. A category with no matches may be absent.
#' @param decrange Range of declination to consider.
#' @param timerange Earlier argument name for \code{year}, retained for existing scripts.
#' @param year Temporal range to consider. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param max.mag Maximum magnitude of stars to consider. Defaults to 2.5
#' @param loc FALSE for geocentric declinations, a vector of latitude, longitude and elevation, or a skyscapeR.horizon to also evaluate the local horizon.
#' @param dates Whether to display time period when identified celestial targets
#' occur within declination range. Default is FALSE.
#' @param calendar Calendar for civil dates. Use \code{"G"} for Gregorian or \code{"J"} for Julian.
#' @export
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   \dontrun{
#'   findTargets(c(-25,-17.5), c(-2500,-1750))
#'   
#'   # if a location is given then the zenith and anti-zenith sun will also be looked at:
#'   findTargets(c(3,12), c(-2500,-1750), loc=c(8.6, 7.3, 200))
#'   
#'   # if a horizon profile is given then the spatial equinox will also be looked at:
#'   hor <- malta_horizons[[1]]
#'   findTargets(c(-7,2), c(-2500,-1750), loc=hor)
#'   }
#' }
#' @section Candidate selection and visibility:
#' A match is a declination overlap, not evidence that a celestial object was
#' visible, deliberately observed or associated with the archaeological feature.
#' The magnitude filter uses fixed catalogue values; it includes integrated
#' magnitudes for non-stellar entries such as M45. It does not predict naked-eye
#' visibility through extinction or twilight, or resolve the members of a cluster.
#'
#' Choose the dating interval and measurement uncertainty from the archaeological
#' and survey evidence. Then examine candidate stars with \code{\link{star.phases}},
#' \code{\link{riseset}} and \code{\link{sketchSky}}. For a statistically selected
#' local declination band, candidate screening is exploratory; the significance
#' of that band does not identify its cause. Physical plausibility and independent
#' archaeological evidence are separate parts of the interpretation.
findTargets <- function(decrange, year, loc=FALSE, max.mag=2.5, dates=FALSE,
                        calendar=skyscapeR.env$calendar, timerange) {
  if (!missing(timerange)) {
    if (!missing(year)) stop("Supply year only once.",call.=FALSE)
    year <- timerange
  }
  if (is.character(year)) year <- astroYear(year)
  targets <- c()

  years <- seq(min(year), max(year))
  seasons <- vapply(unique(range(years)), .season_dates, character(4), calendar=calendar)
  labels <- apply(seasons,1,function(x) paste(unique(x),collapse=" / "))

  ## solar
  targets$solar <- data.frame(name='Test', min.dec=12, max.dec=12, date1='1 Jan', date2='1 Feb', stringsAsFactors=F)
  targets$solar[1,] <- c('june solstice', range(jS(seq(min(year), max(year),1), loc, verbose=FALSE)), labels[2], NA)
  targets$solar[2,] <- c('december solstice', range(dS(seq(min(year), max(year),1), loc, verbose=FALSE)), labels[4], NA)
  targets$solar[3,] <- c('astronomical equinox', 0, 0, labels[1], labels[3])
  if (class(loc)[1] == 'skyscapeR.horizon') {
    aux <- solar.date(sort(as.numeric(spatial.equinox(loc)$declination)), mean(year), calendar, verbose=F)
    targets$solar[4,] <- c('spatial equinox', as.numeric(aux[1,]), paste(aux[2,1],'/',aux[2,2]), paste(aux[3,2],'/',aux[3,1]))
  }
  if(class(loc)[1] != 'logical') {
    if (!is.null(zenith(loc))) {
      aux <- solar.date(sort(as.numeric(zenith(loc))), mean(year), calendar, verbose=F)
      targets$solar[5,] <- c('zenith sun', as.numeric(aux[1,]), NA, aux[2,], aux[3,])
      aux <- solar.date(sort(as.numeric(antizenith(loc))), mean(year), calendar, verbose=F)
      targets$solar[6,] <- c('antizenith sun', as.numeric(aux[1,]), NA, aux[2,], aux[3,])
    }
  }

  aux1 <- .try_target(solar.date(seq(min(decrange),max(decrange),0.01), min(year), calendar, verbose=F))
  aux2 <- .try_target(solar.date(seq(min(decrange),max(decrange),0.01), max(year), calendar, verbose=F))

  ttt <- c()
  if (class(aux1)[1] != 'try-error') { ttt <- c(ttt, aux1[-1,,drop=FALSE]) }
  if (class(aux2)[1] != 'try-error') { ttt <- c(ttt, aux2[-1,,drop=FALSE]) }
  ttt <- ttt[!is.na(ttt)]
  if (length(ttt)) {
    # A leap-year reference preserves 29 February and all month ends.
    parts <- strsplit(trimws(ttt)," +")
    day <- as.integer(vapply(parts,`[`,character(1),1L))
    month <- match(vapply(parts,`[`,character(1),2L),month.abb)
    template <- as.Date("2000-01-01") + 0:365
    ordinal <- sort(unique(match(sprintf("%02d-%02d",month,day),format(template,"%m-%d"))))
    groups <- split(ordinal,cumsum(c(TRUE,diff(ordinal)>2)))
    if (length(groups)>1L && 366L %in% ordinal && 1L %in% ordinal) {
      groups[[1]] <- c(groups[[length(groups)]],groups[[1]])
      groups <- head(groups,-1L)
    }
    windows <- vapply(groups,function(g) {
      endpoints <- as.POSIXlt(template[c(g[1],tail(g,1L))],tz="UTC")
      paste(paste(endpoints$mday,month.abb[endpoints$mon+1L]),collapse=" - ")
    },character(1))
    targets$solar[7,] <- c("sunrise/set",range(decrange),windows[1],
      if(length(windows)>1L) paste(windows[-1L],collapse=" / ") else NA_character_)
  }

  targets$solar[,2] <- as.numeric(targets$solar[,2])
  targets$solar[,3] <- as.numeric(targets$solar[,3])

  targets$solar$min.dec <- round(targets$solar$min.dec,2)
  targets$solar$max.dec <- round(targets$solar$max.dec,2)

  targets$solar <- targets$solar[!is.na(targets$solar[,1]),,drop=FALSE]
  for (i in seq_len(NROW(targets$solar))) {
    if ((min(decrange) > max(targets$solar[i,c(2,3)], na.rm=T))  | (max(decrange) < min(targets$solar[i,c(2,3)], na.rm=T))) {
      targets$solar[i,] <- c(NA,NA,NA,NA,NA)
    }
  }
  targets$solar <- targets$solar[!is.na(targets$solar[,1]),,drop=FALSE]
  rownames(targets$solar) <- NULL
  if (NROW(targets$solar) ==0) { targets$solar <- c() }



  ## lunar
  targets$lunar <- data.frame(name='Test', min.dec=12, max.dec=12, stringsAsFactors=F)

  time <- seq(min(year), max(year),1)
  aux <- sMjLX(time, loc)


  targets$lunar[1,] <- c('southern major lunar extreme', range(sMjLX(seq(min(year), max(year),1), loc, verbose=FALSE)))
  targets$lunar[2,] <- c('southern minor lunar extreme', range(smnLX(seq(min(year), max(year),1), loc, verbose=FALSE)))
  targets$lunar[3,] <- c('northern minor lunar extreme', range(nmnLX(seq(min(year), max(year),1), loc, verbose=FALSE)))
  targets$lunar[4,] <- c('northern major lunar extreme', range(nMjLX(seq(min(year), max(year),1), loc, verbose=FALSE)))
  ## TODO add EFM distribution

  targets$lunar[,2] <- as.numeric(targets$lunar[,2])
  targets$lunar[,3] <- as.numeric(targets$lunar[,3])

  targets$lunar$min.dec <- round(targets$lunar$min.dec,2)
  targets$lunar$max.dec <- round(targets$lunar$max.dec,2)

  for (i in seq_len(NROW(targets$lunar))) {
    if ((targets$lunar$min.dec[i] < min(decrange) & targets$lunar$max.dec[i] < min(decrange)) | (targets$lunar$min.dec[i] > max(decrange) & targets$lunar$max.dec[i] > max(decrange))) {
      targets$lunar[i,] <- c(NA,NA,NA)
    }
  }
  targets$lunar <- targets$lunar[!is.na(targets$lunar[,1]),,drop=FALSE]
  rownames(targets$lunar) <- NULL
  if (NROW(targets$lunar) ==0) { targets$lunar <- c() }



  ## stellar
  ss <- .read_star_catalogue()
  ss <- ss[which(ss$magV<=max.mag),,drop=FALSE] # magnitude filter

  targets$stellar <- data.frame(name='Test', identifier='Test', magV=1.2, min.dec=12, max.dec=12, daterange=NA, stringsAsFactors=F)[FALSE,,drop=FALSE]
  for (i in seq_len(NROW(ss))) {
    d <- star(ss$identifier[i],years)$coord$Dec
    bounds <- range(d)
    if (bounds[2]<min(decrange) || bounds[1]>max(decrange)) next
    interval <- if (dates) .stellar_dates(years,d,range(decrange)) else NA_character_
    targets$stellar[nrow(targets$stellar)+1L,] <- list(ss$name[i],ss$identifier[i],
      ss$magV[i],round(bounds[1],2),round(bounds[2],2),interval)
  }
  for (i in seq_len(NROW(targets$stellar))) {
    aux <- star.names$Western[which(star.names$identifier == targets$stellar[i,2])]
    if (length(aux)>0) { targets$stellar[i,1] <- aux } else { targets$stellar[i,1] <- '' }
  }

  rownames(targets$stellar) <- NULL
  if (NROW(targets$stellar) ==0) { targets$stellar <- c() }

  ## for "security" reasons
  if (dates==FALSE) targets$stellar$daterange <- NULL

  return(targets)
}
