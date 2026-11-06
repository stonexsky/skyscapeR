#' Prepare a reusable sky sketch for a horizon plot
#' @description Calculate the sky above one observing site at one instant, then reuse it
#' as the background of a horizon plot. Natural style illustrates daylight, twilight,
#' moonlight and atmospheric extinction. Publication styles show the selected objects
#' without brightness-dependent visibility filtering.
#' @param hor A \code{skyscapeR.horizon} with finite latitude, longitude and observer
#' elevation in metres in \code{metadata$georef}. The numerical profile supplies the
#' foreground when plotted. A panorama image is optional.
#' @param jd One finite UT Julian day. Supply exactly one of \code{jd} and \code{time}.
#' @param time A civil timestamp accepted by \code{\link{time2jd}}, or a POSIXt instant.
#' @param timezone,calendar Civil timezone and calendar for character \code{time}, as
#' in \code{\link{time2jd}}. Defaults are the session settings. A POSIXt value already
#' identifies an instant; these arguments do not reinterpret it. Use UTC or a fixed
#' offset for ancient civil dates.
#' @param display A named list of logical switches: \code{sun}, \code{moon},
#' \code{planets}, \code{stars}, \code{constellation.lines}, \code{milky.way},
#' \code{galactic.equator}, \code{ecliptic}, \code{equator}, \code{poles} and
#' \code{zenith}. Sun, Moon and stars default to TRUE; other switches default to
#' FALSE, except that the Milky Way follows \code{stars} unless explicitly set.
#' The additional \code{constellations} entry accepts names or an imported data
#' frame, rather than a logical value. See Constellations below. A named logical
#' vector is sufficient when supplying only switches.
#' @param style One of \code{"natural"}, \code{"publication-colour"} or
#' \code{"publication-mono"}. Publication styles have a white background and dark
#' guides. Colour uses muted stellar and planetary colours and a warm solar disc;
#' mono uses black and white. Neither publication style predicts visibility.
#' @param max.mag Finite catalogue magnitude cutoff, default 6. Applied before
#' extinction to star symbols and the integrated magnitudes of supported clusters.
#' It does not limit constellation endpoints, the Sun, Moon or planets. Raising
#' this value cannot supply stars absent from the catalogue or overcome the
#' symbol-size cutoff described below.
#' @param k Non-negative visual extinction coefficient in magnitudes per airmass
#' for natural style. Defaults to \code{\link{skyscapeR.vars}()$k}, initially 0.2.
#' This clear-sky baseline is an assumption, not an estimate of the site's weather.
#' Zero removes extinction and the modelled lunar scattering. Publication styles
#' do not apply extinction to symbols.
#' @details The function computes an object without drawing. Use
#' \code{plot(hor, sky = sky, xrange = ..., yrange = ...)} to choose the sector.
#' Replotting the object changes framing without recalculating the instant. Use
#' \code{\link{skyStep}} to recalculate at another time. Astronomical calculations
#' require Swiss Ephemeris files covering the epoch; see
#' \code{\link{skyscapeR-package}} for installation.
#' @section Positions, atmosphere and symbol sizes:
#' Swiss Ephemeris supplies apparent topocentric stellar and planetary positions,
#' solar and lunar angular diameters, and lunar phase. Atmospheric refraction is
#' separate from astronomical apparent-place corrections. Refraction, pressure and
#' temperature come from \code{\link{skyscapeR.vars}} when the sketch is calculated;
#' the object records these settings. Azimuth increases clockwise from true north.
#'
#' Star symbols use fixed catalogue visual magnitudes and spectral-class colours.
#' There is no variable-star light curve or time-dependent cluster magnitude.
#' Selected combined-system entries avoid duplicate symbols for Dubhe and Zeta
#' Herculis; other catalogue components can remain distinct. Non-stellar catalogue
#' entries are omitted from the point-star layer but remain available to
#' \code{\link{findTargets}} and \code{\link{star.phases}}.
#'
#' Stellar and planetary point sizes use the package display mapping
#' \code{cex = 4.45 - 1.57 * log(m + 10)}, where \code{m} is the catalogue or
#' planetary magnitude in publication styles and the extinction-adjusted magnitude
#' in natural style. Magnitudes greater than 7 have zero size. Symbols are graphical
#' markers, not angular stellar diameters or calibrated retinal images. The five
#' planets are Mercury, Venus, Mars, Jupiter and Saturn; their magnitudes are
#' calculated for the instant. The Sun and Moon use angular disc diameters rather
#' than this point-size mapping. The illuminated lunar limb faces the Sun even
#' when the Sun is hidden. Discs are drawn in azimuth-altitude coordinates, so wide
#' fields and high altitudes introduce projection distortion.
#' @section Constellations:
#' Use \code{display = list(constellations = "Orion")} to select and enable line
#' figures. \code{"modern"} selects all 88 Western figures, \code{"chinese"} the
#' full Chinese set, and \code{"chinese_mansions"} the 28 lunar-mansion figures.
#' Names and abbreviations are case-sensitive; for example \code{c("Ori", "Tau")}.
#' Alternatively, supply a data frame from \code{\link{importConstellations}}.
#' See \code{\link{constellations}} for sources and identifiers.
#'
#' Supplying a selection enables lines unless \code{constellation.lines = FALSE}
#' is also supplied. Setting \code{constellation.lines = TRUE} without a selection
#' uses modern figures. Endpoint coordinates are calculated independently of the
#' star-symbol magnitude cut. Lines therefore remain when their stars are faint,
#' extincted, or switched off. A single-star figure is drawn as a ring. Normal
#' plot clipping and the horizon foreground still apply. These are reference
#' figures; their visibility is not an observing prediction.
#' @section Stellar extinction and twilight visibility:
#' Natural style adds \eqn{k X(h)} to magnitude, using Rozenberg's approximation
#' \eqn{X(h) = 1 / (\sin h + 0.025 \exp(-11 \sin h))}, with angles converted
#' from degrees. This is the extinction-airmass expression discussed by Krisciunas
#' and Schaefer (1991, equation 14). Stellar point symbols use apparent altitude for
#' this optical path calculation even if the display uses geometric altitude.
#' Altitude is clamped at zero for this calculation. Planet symbols use their
#' plotted altitude. Near-horizon extinction remains sensitive to the real atmosphere.
#'
#' For effective sky brightness \eqn{B} in nanolamberts, the display limiting
#' magnitude is \eqn{6 - 1.45 \log_{10}(B/60)}. Symbols fade smoothly over the last
#' 0.75 magnitude below this limit, with zero opacity at or beyond the limit.
#' The smooth opacity is \eqn{f^2(3-2f)}, where \eqn{f} is the clipped fraction
#' \eqn{(m_{lim}-m)/0.75}. This is a package display rule, not a published visual
#' threshold law. Symbols are suppressed when the altitude used for visibility testing is
#' negative.
#'
#' During twilight, stars and cluster centres use a shared adaptation term based
#' on solar zenith brightness plus local scattered moonlight. Relative to the dark
#' baseline of 60 nanolamberts, this brightness is raised to the power 0.40.
#' Its weight increases linearly from zero at solar altitude -12 degrees to one
#' at -9 degrees, stays one to -5 degrees, and falls to zero at -1 degree.
#' Local and adapted brightness are blended logarithmically. This empirical
#' wide-field display adjustment was informed by comparisons with Stellarium;
#' it is not Stellarium's implementation or a universal calibration. Planets use
#' the local brightness calculation without this stellar adaptation adjustment.
#' @section Solar and lunar sky brightness:
#' Solar zenith brightness is interpolated logarithmically between representative
#' values at solar altitudes -18, -12, -6, -0.833, 15 and 90 degrees: respectively
#' 60, 240, 100000, 3000000, 300000000 and 1000000000 nanolamberts. The scale is
#' informed by Schaefer (1993); the interpolation and angular gradient are package
#' approximations. Solar altitude for twilight is geometric. Sky colours and the
#' Belt of Venus are illustrative approximations informed by Lee (1994, 2015),
#' not calibrated photometry.
#'
#' Lunar scattering follows Krisciunas and Schaefer (1991, equations 20-23), using
#' lunar phase angle, angular separation, viewing altitude and \code{k}. It uses
#' their scattering-airmass expression, distinct from the stellar extinction
#' expression. Separation is bounded below at 10 degrees to avoid extrapolation
#' near the disc. Only a Moon above the astronomical horizontal contributes.
#' A ridge can hide the disc without removing its modelled sky illumination.
#' Terrain does not cast shadows into the atmospheric scattering calculation.
#'
#' The Sun and Moon switches control both their discs and their respective sky
#' illumination. With both off, natural style uses its dark baseline while
#' extinction remains active. Use \code{display = list(moon = FALSE)} for a
#' moonless comparison. These are deliberate display scenarios rather than
#' changes to the calculated astronomical instant.
#' @section Unresolved clusters:
#' In natural style with stars enabled, 17 catalogue-matched open and globular
#' clusters can contribute a soft patch. Angular diameters come from OpenNGC;
#' missing minor diameters use the major diameter. An equal-area circle has radius
#' \code{sqrt(major.arcmin * minor.arcmin) / 120} degrees. Spherical angular
#' separation defines its footprint, including across north and near the zenith.
#' The integrated catalogue magnitude supplies an upper visibility bound through
#' the stellar display rule. No integrated M45 point symbol is drawn.
#'
#' The patch fades as visible point symbols account for projected light within
#' its footprint. This is a flux-weighted estimate using all projected stars,
#' including possible foreground overlaps, not a cluster-membership calculation.
#' The central opacity is \code{0.28 * integrated * remaining^3 / sqrt(1 + radius^2)},
#' with a smooth radial taper to zero at the edge. These constants are rendering
#' choices. Publication styles omit diffuse patches and retain individual stars.
#'
#' Integrated magnitude alone does not determine extended-object visibility.
#' Surface brightness, angular size and contrast matter (Crumey 2014). This
#' rendering does not implement Crumey's threshold model, resolve actual cluster
#' structure, or establish a first-detection date. Only clusters in the supplied
#' size table receive a patch; it is not a complete deep-sky renderer.
#' @section Milky Way and reference guides:
#' The schematic Milky Way outline is adapted from NASA's Deep Star Maps 2020,
#' credited to NASA/Goddard Scientific Visualization Studio and ESA/Gaia/DPAC.
#' Its fixed shape is rotated to the observing epoch; stellar evolution and the
#' Galaxy's physical evolution are not reconstructed. The main dark rift is
#' retained, but the band is not a calibrated surface-brightness map.
#'
#' In natural style, solar fading decreases from full strength at -18 degrees
#' to zero at -12 degrees. Local sky brightness reduces contrast, and both fill
#' and outline are attenuated with altitude by transmission
#' \eqn{10^{-0.4 k X(h)}}. Band altitude is clamped at zero for airmass. Publication
#' styles show the schematic band without atmospheric dimming. The Galactic
#' equator is a separate dash-dot guide. The ecliptic is dashed and the celestial
#' equator dotted; these reference curves use geometric altitude. Pole markers
#' are geometric and the zenith is labelled at the centre of the plot width.
#' @section Interpretation and limitations:
#' Twilight thresholds and unresolved-cluster glow are display approximations,
#' not a universally calibrated naked-eye detection model. Agreement in selected
#' comparisons with Stellarium or \code{\link{star.phases}} does not establish
#' agreement across other objects, sites, epochs or atmospheric conditions.
#' Neither the presence nor absence of a cluster patch establishes its first
#' or last detection date, and its opacity is not measured surface brightness.
#' Report the visibility assumptions separately when interpreting a proposed
#' observation.
#'
#' A sketch is a visualisation, not a test of naked-eye detectability. There is no
#' observer-specific acuity, adaptation history, light-pollution map, cloud model,
#' wavelength-dependent extinction, or variable-star brightness model. Very low
#' altitudes, twilight thresholds and narrow terrain clearances deserve particular
#' care. Screen size, symbol size and projection also affect what appears visible.
#'
#' \code{\link{star.phases}} uses a different visibility criterion: a fixed
#' limiting magnitude, atmospheric extinction and a star-Sun altitude-separation
#' rule on a synthetic annual cycle. Shared positions and extinction do not make
#' its event dates identical to this display's threshold. A rise time from
#' \code{\link{riseset}} is a centre crossing of the skyline, not first visibility.
#' Inspect a short sequence after crossing with \code{\link{skyStep}} while the
#' sky is still dark enough. An absent sketch symbol alone cannot rule out an
#' observationally plausible heliacal sighting. Cluster predictions based on
#' integrated magnitude require additional caution.
#'
#' For reproducibility, retain the sky object, horizon and external image, package
#' and ephemeris versions, civil calendar/timezone, and atmosphere settings.
#' Recalculation with \code{skyStep} reads the session refraction, pressure and
#' temperature again; keep these fixed when comparing time steps.
#' @return An object of classes \code{skyscapeR.sketch} and \code{skyscapeR.sky}.
#' It stores \code{jd}, \code{loc}, \code{style}, \code{display}, \code{settings}
#' (calendar, magnitude cutoff and extinction), and \code{atmosphere}.
#' \code{stars} and \code{coords} contain azimuth, altitude and symbol size;
#' \code{star.opacity} records display opacity. \code{bodies}, \code{diameters},
#' \code{planet.size} and \code{lunar} describe selected Solar System bodies.
#' \code{clusters} holds supported cluster footprints and opacity when applicable;
#' \code{guides}, \code{light}, \code{background} and twilight fields support drawing.
#' Components can be empty when a layer is disabled. No image file is written.
#' @references Krisciunas, K. and Schaefer, B. E. (1991). A model of the brightness
#' of moonlight. PASP 103, 1033-1039. \doi{10.1086/132921}.
#' @references Schaefer, B. E. (1986). Atmospheric Extinction Effects on Stellar
#' Alignments. Journal for the History of Astronomy 17, S32-S42.
#' \doi{10.1177/002182868601701003}.
#' @references Schaefer, B. E. (1993). Astronomy and the Limits of Vision.
#' Vistas in Astronomy 36, 311-361. \doi{10.1016/0083-6656(93)90113-X}.
#' @references Lee, R. L., Jr. (1994). Twilight and daytime colors of the clear sky.
#' Applied Optics 33, 4629-4638. \doi{10.1364/AO.33.004629}.
#' @references Lee, R. L., Jr. (2015). Measuring and modeling twilight's Belt of
#' Venus. Applied Optics 54, B194-B203. \doi{10.1364/AO.54.00B194}.
#' @references Crumey, A. (2014). Human contrast threshold and astronomical
#' visibility. MNRAS 442, 2600-2619. \doi{10.1093/mnras/stu992}.
#' @references NASA/Goddard Space Flight Center Scientific Visualization Studio
#' (2020). Deep Star Maps 2020. \url{https://svs.gsfc.nasa.gov/4851/}.
#' @references Verga, M. and contributors. OpenNGC. Source of cluster angular
#' sizes. \url{https://github.com/mattiaverga/OpenNGC}. See installed COPYRIGHTS
#' for data attribution and licences.
#' @references Astrodienst. Swiss Ephemeris programmer's manual, fixed stars,
#' horizontal coordinates and planetary phenomena.
#' \url{https://www.astro.com/swisseph/swephprg.htm}.
#' @seealso \code{\link{star.phases}}, \code{\link{riseset}},
#' \code{\link{findTargets}}, \code{\link{skyStep}},
#' \code{\link{stellariumSky}}, \code{\link{plot.skyscapeR.horizon}},
#' \code{\link{importConstellations}}, \code{\link{skyscapeR.vars}}
#' @examples
#' \dontrun{
#' hor <- malta_horizons[[1]]
#' sector <- c(0, 180)
#' sky <- sketchSky(hor, time = "2026/12/21 20:00", timezone = "UTC")
#' plot(hor, sky = sky, xrange = sector)
#' sky <- skyStep(hor, sky, minutes = 10, xrange = sector)
#'
#' # Compare publication styles at exactly the same instant.
#' colour <- sketchSky(hor, jd = sky$jd, style = "publication-colour",
#'                     display = list(constellations = "modern"))
#' mono <- sketchSky(hor, jd = sky$jd, style = "publication-mono",
#'                   display = list(constellations = "modern"))
#' plot(hor, sky = colour, xrange = sector)
#' plot(hor, sky = mono, xrange = sector)
#' }
#' @export
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
sketchSky <- function(hor, jd=NULL, time=NULL, timezone, calendar,
                      display=list(), style=c("natural","publication-colour","publication-mono"), max.mag=6, k=skyscapeR.vars()$k) {
  style <- match.arg(style)
  if (!inherits(hor,"skyscapeR.horizon")) stop("Supply a georeferenced site horizon.")
  loc <- hor$metadata$georef
  if (!is.numeric(loc) ||
      length(loc)!=3L || any(!is.finite(loc))) stop("Supply a georeferenced site horizon.")
  if (missing(calendar)) calendar <- skyscapeR.env$calendar
  jd <- .sky_instant(jd,time,timezone,calendar)
  flags <- c(sun=TRUE,moon=TRUE,planets=FALSE,stars=TRUE,constellation.lines=FALSE,
             milky.way=FALSE,galactic.equator=FALSE,ecliptic=FALSE,equator=FALSE,poles=FALSE,zenith=FALSE)
  if (anyDuplicated(names(display))) stop("Use unique names in display.")
  display <- as.list(display)
  selection <- display[["constellations"]]
  if (is.null(selection)) {
    selection <- "modern"
  } else {
    display[["constellations"]] <- NULL
    if (is.null(display[["constellation.lines"]])) display$constellation.lines <- TRUE
  }
  if (length(display)) {
    if (is.null(names(display)) || anyDuplicated(names(display)) ||
        !all(names(display) %in% names(flags)) ||
        !all(vapply(display,function(x) is.logical(x) && length(x)==1L && !is.na(x),logical(1))))
      stop("Use uniquely named logical display switches. See ?sketchSky.")
    flags[names(display)] <- unlist(display,use.names=FALSE)
  }
  if (!"milky.way" %in% names(display)) flags["milky.way"] <- flags["stars"]
  if (!is.numeric(max.mag) || length(max.mag)!=1L || !is.finite(max.mag)) stop("Supply one magnitude limit.")
  if (!is.numeric(k) || length(k)!=1L || !is.finite(k) || k<0) stop("Supply a non-negative extinction coefficient k.")
  atmosphere <- skyscapeR.vars()[c("refraction","atm","temp")]
  position <- function(id) unlist(body.position(jd,obj=id,loc=loc,verbose=FALSE)$horizontal,use.names=FALSE)
  solar <- if (flags["sun"] || flags["moon"]) position("sun")
  visibility <- 1
  twilight <- background <- NULL
  if (style=="natural") twilight <- "night"
  if (style=="natural" && flags["sun"]) {
    altitude <- if (atmosphere$refraction) {
      body.position(jd,obj="sun",loc=loc,refraction=FALSE,verbose=FALSE)$horizontal$alt
    } else solar[2]
    visibility <- .milky_way_visibility(altitude)
    twilight <- .sketch_twilight(altitude)
  }
  light <- NULL
  bodies <- list()
  diameters <- planet.size <- numeric()
  if (flags["sun"]) {
    bodies$sun <- solar
    diameters["sun"] <- .sketch_pheno(jd,"sun",loc)[4]
  }
  lunar <- NULL
  if (flags["moon"]) {
    bodies$moon <- position("moon")
    pheno <- .sketch_pheno(jd,"moon",loc)
    diameters["moon"] <- pheno[4]
    lunar <- list(fraction=pheno[2], angle=.lunar_light_angle(bodies$moon,solar))
  }
  if (style=="natural") {
    light <- list(sun=if (flags["sun"]) c(solar[1],altitude),moon=bodies$moon,
                  phase=if (flags["moon"]) pheno[1],k=k)
    background <- .sky_colour(light,0,90,opaque=TRUE)
  }
  planet.opacity <- numeric()
  if (flags["planets"]) {
    ids <- c("mercury","venus","mars","jupiter","saturn")
    bodies[ids] <- lapply(ids,position)
    for (id in ids) {
      pos <- bodies[[id]]
      mag <- .sketch_pheno(jd,id,loc)[5]
      if (!is.null(light)) {
        mag <- mag + .stellar_extinction(max(0,pos[2]),k)
        planet.opacity[id] <- .point_visibility(light,pos[1],pos[2],mag)
      }
      planet.size[id] <- mag2size(mag)
    }
  }

  star.colours <- NULL
  stars <- character()
  edges <- NULL
  clusters <- NULL
  if (flags["stars"]) {
    catalogue <- .read_star_catalogue()
    if (style == "natural") {
      sizes <- utils::read.csv(system.file("extdata","cluster-sizes.csv",package="skyscapeR"))
      clusters <- merge(sizes,catalogue[,c("identifier","magV")],by="identifier")
      clusters <- clusters[is.finite(clusters$magV) & clusters$magV <= max.mag,,drop=FALSE]
    }
    catalogue <- .plot_star_catalogue(catalogue)
    keep <- which(catalogue$magV<=max.mag)
    stars <- as.character(catalogue$identifier[keep])
    if (length(keep) && !is.null(catalogue$stellar)) star.colours <- .sketch_star_colours(catalogue$stellar[keep],catalogue$magV[keep])
  }
  if (flags["constellation.lines"]) edges <- do.call(rbind,constLines(selection))
  ids <- unique(as.character(c(stars,edges,clusters$identifier)))
  ids <- ids[!is.na(ids)]
  opacity <- stats::setNames(numeric(length(ids)),ids)
  coords <- matrix(NA_real_,length(ids),3,dimnames=list(ids,c("az","alt","size")))
  if (length(ids)) {
    .file_position(jd,0L,2L)
    swephR::swe_set_topo(loc[2],loc[1],loc[3])
    adaptation <- 0
    if (!is.null(light$sun)) {
      adaptation <- stats::approx(c(-12,-9,-5,-1),c(0,1,1,0),xout=light$sun[2],rule=2)$y
    }
    for (i in seq_along(ids)) {
      key <- paste0(",",ids[i])
      target <- swephR::swe_fixstar2_ut(key,jd,32768L+2048L+2L)
      if (target$return<0L) next
      if (bitwAnd(target$return,7L)!=2L) .ephemeris_failure(jd,key,target$serr)
      pos <- swephR::swe_azalt(jd,1L,loc[c(2,1,3)],atmosphere$atm,atmosphere$temp,target$xx[1:3])$xaz
      magnitude <- swephR::swe_fixstar2_mag(key)$mag
      opacity[i] <- 1
      if (style=="natural") {
        # Optical path length depends on apparent altitude even on a geometric plot.
        magnitude <- magnitude + .stellar_extinction(max(0,pos[3]),k)
        opacity[i] <- .point_visibility(light,(pos[1]+180)%%360,pos[3],magnitude,adaptation)
      }
      coords[i,] <- c((pos[1]+180)%%360,pos[if (atmosphere$refraction) 3L else 2L],mag2size(magnitude))
    }
  }
  if (NROW(clusters)) {
    centres <- coords[clusters$identifier,,drop=FALSE]
    clusters$az <- centres[,1]
    clusters$alt <- centres[,2]
    # A circle with the catalogue area avoids assuming a position angle.
    clusters$radius <- sqrt(clusters$major.arcmin * clusters$minor.arcmin) / 120
    extinction <- .stellar_extinction(pmax(0,clusters$alt),k)
    integrated <- .point_visibility(light,clusters$az,clusters$alt,
                                    clusters$magV + extinction,adaptation)
    clusters$alpha <- numeric(nrow(clusters))
    points <- coords[stars,,drop=FALSE]
    for (i in which(is.finite(integrated) & integrated > 0)) {
      distance <- .cluster_distance(clusters$az[i],clusters$alt[i],points[,1],points[,2])
      inside <- which(distance <= clusters$radius[i] & is.finite(distance))
      # Fade the patch as point symbols account for the light in its footprint.
      # This is a projected-light estimate, not a membership assignment.
      fraction <- sum(opacity[stars[inside]] *
                        10^(-.4 * (catalogue$magV[keep[inside]] - clusters$magV[i])),
                      na.rm=TRUE)
      remaining <- max(0, 1 - fraction)
      # The integrated magnitude sets an upper visibility bound. Dilute the
      # display contrast of light spread over a wider area.
      clusters$alpha[i] <- .28 * integrated[i] * remaining^3 /
                           sqrt(1 + clusters$radius[i]^2)
    }
  }
  guides <- list()
  for (id in intersect(c("ecliptic","equator"),names(flags)[flags])) {
    frame <- if (id=="ecliptic") 0L else 1L
    guides[[id]] <- t(vapply(seq(0,360,by=.5),function(lon) {
      pos <- swephR::swe_azalt(jd,frame,loc[c(2,1,3)],atmosphere$atm,atmosphere$temp,c(lon,0,1))$xaz
      c((pos[1]+180)%%360,pos[2])
    },numeric(2)))
  }
  outline <- flags["milky.way"] && visibility > 0
  if (outline || flags["galactic.equator"]) {
    paths <- .galactic_paths(jd, outline, flags["galactic.equator"])
    for (id in names(paths)) {
      apparent <- id != "galactic.equator" && atmosphere$refraction
      guides[[id]] <- t(apply(paths[[id]], 1L, function(eq) {
        pos <- swephR::swe_azalt(jd,1L,loc[c(2,1,3)],atmosphere$atm,atmosphere$temp,c(eq,1))$xaz
        c((pos[1]+180)%%360,pos[if (apparent) 3L else 2L])
      }))
    }
  }
  .with_astronomy_provenance(structure(list(clusters=clusters,light=light,planet.opacity=planet.opacity,star.opacity=opacity,milky.way.visibility=visibility,guides=guides,diameters=diameters,planet.size=planet.size,star.colours=star.colours,
                 background=background, settings=list(calendar=calendar,max.mag=max.mag,k=k),
                 jd=jd,loc=loc,display=c(as.list(flags),list(constellations=selection)),style=style,day=style=="natural" && flags["sun"] && solar[2]>0,twilight=twilight,atmosphere=atmosphere,
                 lunar=lunar,stars=coords[match(stars,ids),,drop=FALSE],edges=edges,coords=coords,bodies=bodies),
            class=c("skyscapeR.sketch","skyscapeR.sky")))
}

.sketch_pheno <- function(jd,id,loc) {
  body <- as.integer(checkbody(id))
  .file_position(jd,body,2L)
  swephR::swe_set_topo(loc[2],loc[1],loc[3])
  pheno <- swephR::swe_pheno_ut(jd,body,32768L+2L)
  if (pheno$return<0L || bitwAnd(pheno$return,7L)!=2L)
    .ephemeris_failure(jd,id,pheno$serr)
  pheno$attr
}

.sketch_star_colours <- function(type,mag) {
  palette <- c(O="#ADCAFF",B="#BECFFF",A="#EFF4FF",F="#FFF8EC",
               G="#FFE9BA",K="#FFC38E",M="#FFAC95")
  tint <- unname(palette[substr(trimws(type),1,1)])
  tint[is.na(tint)] <- "white"
  weight <- pmax(0,pmin(1,(6-mag)/5))
  rgb <- 255-(255-grDevices::col2rgb(tint))*rep(weight,each=3)
  grDevices::rgb(rgb[1,],rgb[2,],rgb[3,],maxColorValue=255)
}

# Clear-sky luminance in nanolamberts. The solar angular gradient is illustrative.
.sky_luminance <- function(light,az,alt,adaptation=0) {
  rad <- pi/180
  alt <- pmax(0,pmin(90,alt))
  separation <- function(p) {
    cosine <- sin(alt*rad)*sin(p[2]*rad)+cos(alt*rad)*cos(p[2]*rad)*cos((az-p[1])*rad)
    acos(pmin(1,pmax(-1,cosine)))/rad
  }
  solar <- 0*(az+alt)
  if (!is.null(light$sun)) {
    zenith <- 10^stats::approx(c(-18,-12,-6,-.833,15,90),
      log10(c(60,240,1e5,3e6,3e8,1e9)),xout=light$sun[2],rule=2)$y
    solar <- (zenith-60)*exp(1.4*(cos(separation(light$sun)*rad)-sin(light$sun[2]*rad)))

  }
  lunar <- 0
  if (!is.null(light$moon) && light$moon[2]>0 && light$k>0) {
    # Krisciunas & Schaefer (1991), equations 20-23. Bound the near-disc extrapolation.
    rho <- pmax(10,separation(light$moon))
    scatter <- 10^5.36*(1.06+cos(rho*rad)^2)+10^(6.15-rho/40)
    intensity <- 10^(-.4*(3.84+.026*light$phase+4e-9*light$phase^4))
    transmission <- function(h) 10^(-.4*light$k/sqrt(1-.96*cos(h*rad)^2))
    lunar <- scatter*intensity*transmission(light$moon[2])*(1-transmission(alt))
  }
  brightness <- 60+pmax(0,solar)+lunar
  if (adaptation > 0 && !is.null(light$sun)) {
    local <- log(brightness/60)
    adapted <- 0.40*log((zenith+lunar)/60)
    brightness <- 60*exp((1-adaptation)*local + adaptation*adapted)
  }
  brightness
}

.point_visibility <- function(light,az,alt,mag,adaptation=0) {
  brightness <- .sky_luminance(light,az,alt,adaptation)
  # Smooth display threshold anchored at magnitude six under a dark sky.
  limit <- 6-1.45*log10(brightness/60)
  fade <- pmax(0,pmin(1,(limit-mag)/.75))
  ifelse(alt<0,0,fade^2*(3-2*fade))
}

.sky_colour <- function(light,az,alt,opaque=FALSE) {
  level <- log10(.sky_luminance(light,az,alt)/60)
  palette <- grDevices::colorRamp(c("#546789","#708BAE","#A8CCE7","lightblue"))
  channel <- palette(pmin(1,level/7))
  if (!is.null(light$sun) && light$sun[2] > -12 && light$sun[2] < 10) {
    h <- pmin(90,pmax(0,alt))
    solar.alt <- light$sun[2]
    # Smooth screen-colour approximation. Directional terms vanish at the zenith.
    twilight <- stats::approx(c(-12,-5,-2,0,10),c(0,.65,1,.9,0),solar.alt)$y
    direction <- cos((az-light$sun[1])*pi/180)*cos(h*pi/180)
    warm <- twilight*pmax(0,direction)^2*exp(-h/15)
    belt <- exp(-((solar.alt+2)/2)^2)*pmax(0,-direction)^2*
      exp(-((h-(3-1.5*solar.alt))/4)^2)
    lab <- grDevices::convertColor(channel/255,"sRGB","Lab")
    lab[,2] <- lab[,2]+24*warm+17*belt
    lab[,3] <- lab[,3]+58*warm+12*belt-7*twilight*sin(h*pi/180)
    channel <- 255*grDevices::convertColor(lab,"Lab","sRGB")
  }
  alpha <- 1-10^(-level*.28)
  if (opaque) {
    base <- matrix(rep(c(16,23,41),each=length(alpha)),ncol=3)
    channel <- base*(1-alpha)+channel*alpha
    alpha <- 1
  }
  grDevices::rgb(channel[,1],channel[,2],channel[,3],255*alpha,maxColorValue=255)
}

.sky_veil <- function(light,box) {
  az <- seq(box[1],box[2],length.out=241)
  alt <- seq(box[4],box[3],length.out=121)
  grid <- expand.grid(alt=alt,az=az)
  rgba <- .sky_colour(light,grid$az,grid$alt)
  graphics::rasterImage(grDevices::as.raster(matrix(rgba,length(alt))),box[1],box[3],box[2],box[4],interpolate=TRUE)
}

.milky_way_visibility <- function(altitude) {
  pmin(1,pmax(0,(-12-altitude)/6))
}

.sketch_twilight <- function(altitude) {
  stages <- c("night","astronomical","nautical","civil","daylight")
  stages[1L+findInterval(altitude,c(-18,-12,-6,-50/60))]
}

.lunar_light_angle <- function(moon, sun) {
  # Project the solar direction onto the Moon's local right/up tangent axes.
  rad <- pi/180
  da <- (sun[1]-moon[1])*rad
  ma <- moon[2]*rad
  sa <- sun[2]*rad
  atan2(sin(sa)*cos(ma)-cos(sa)*sin(ma)*cos(da), cos(sa)*sin(da))
}

.lunar_lit_polygon <- function(fraction, angle, radius=.25) {
  # Limb and terminator share the two horns. Their separation sets the lit area.
  t <- seq(-pi/2,pi/2,length.out=97)
  u <- radius*c(cos(t),(1-2*fraction)*rev(cos(t)))
  v <- radius*c(sin(t),rev(sin(t)))
  list(x=u*cos(angle)-v*sin(angle), y=u*sin(angle)+v*cos(angle))
}

.draw_sketch_sky <- function(sky) {
  box <- graphics::par("usr")
  paper <- sky$style!="natural"
  colour.paper <- sky$style=="publication-colour"
  bg <- if (paper) "white" else if (sky$day) "lightblue" else "midnightblue"
  if (!is.null(sky$background)) {
    bg <- sky$background
  } else if (!paper && !is.null(sky$twilight)) {
    colours <- c(daylight="lightblue",civil="#7897BD",nautical="#40577F",
                 astronomical="#222F51",night="#101729")
    bg <- colours[[sky$twilight]]
  }
  ink <- if (paper) "black" else "white"
  stars.ink <- ink
  if ((!paper || colour.paper) && length(sky$star.colours)) {
    stars.ink <- sky$star.colours
    if (colour.paper) {
      channels <- grDevices::col2rgb(stars.ink)*0.45
      stars.ink <- grDevices::rgb(channels[1,],channels[2,],channels[3,],maxColorValue=255)
    }
  }
  line.ink <- ink
  if (!paper && length(sky$star.opacity)) {
    if (NROW(sky$stars)) {
      rgb <- grDevices::col2rgb(rep_len(stars.ink,NROW(sky$stars)))
      alpha <- sky$star.opacity[rownames(sky$stars)]
      stars.ink <- grDevices::rgb(rgb[1,],rgb[2,],rgb[3,],255*alpha,maxColorValue=255)
    }
  }
  if (!is.null(sky$light)) bg <- "#101729"
  graphics::rect(box[1],box[3],box[2],box[4],col=bg,border=NA)
  visibility <- sky$milky.way.visibility
  if (is.null(visibility)) visibility <- 1  # sketches saved before twilight fading
  band.k <- if (is.null(sky$light)) 0 else sky$light$k
  .draw_galactic_band(sky$guides$south,sky$guides$north,visibility=visibility,k=band.k)
  if (!is.null(sky$light)) .sky_veil(sky$light,box)
  .draw_cluster_glows(sky$clusters,box)
  for (k in seq(floor(min(box[1:2])/360)-1,ceiling(max(box[1:2])/360)+1)) {
    shift <- 360*k
    for (id in names(sky$guides)) {
      line <- sky$guides[[id]]
      az <- line[,1]
      # Draw each continuous section on either side of the north seam.
      jumps <- which(abs(diff(az)) > 180)
      groups <- cumsum(c(0L,abs(diff(az)) > 180))
      outline <- id %in% c("north","south")
      if (outline && !is.null(sky$light)) {
        dx <- (diff(az)+180)%%360-180
        starts <- seq_along(dx)
        brightness <- .sky_luminance(sky$light,az[starts]+dx/2,(line[starts,2]+line[starts+1,2])/2)
        height <- (line[starts,2]+line[starts+1,2])/2
        alpha <- .55*visibility*(brightness/60)^(-.28)*.band_transmission(height,band.k)
        graphics::segments(az[starts]+shift,line[starts,2],az[starts]+dx+shift,line[starts+1,2],
                           col=grDevices::rgb(1,1,1,alpha),lwd=.7)
        next
      }
      col <- if (outline) grDevices::adjustcolor(ink,alpha.f=.55*visibility) else ink
      lty <- switch(id, ecliptic=2, equator=3, galactic.equator=4, 1)
      for (rows in split(seq_along(az),groups))
        graphics::lines(az[rows]+shift,line[rows,2],col=col,lty=lty,lwd=.7)
      if (length(jumps)) {
        dx <- (az[jumps+1L]-az[jumps]+180)%%360-180
        graphics::segments(az[jumps]+shift,line[jumps,2],
                             az[jumps]+shift+dx,line[jumps+1L,2],col=col,lty=lty,lwd=.7)
      }
    }
    if (isTRUE(sky$display[["poles"]])) {
      graphics::points(c(0,180)+shift,c(1,-1)*sky$loc[1],pch=3,col=ink)
      graphics::text(c(0,180)+shift,c(1,-1)*sky$loc[1],c("NCP","SCP"),pos=4,cex=.7,col=ink)
    }
    if (NROW(sky$edges)) {
      a <- sky$coords[match(sky$edges[,1],rownames(sky$coords)),,drop=FALSE]
      b <- sky$coords[match(sky$edges[,2],rownames(sky$coords)),,drop=FALSE]
      dx <- (b[,1]-a[,1]+180)%%360-180
      graphics::segments(a[,1]+shift,a[,2],a[,1]+shift+dx,b[,2],col=line.ink,lwd=.5)
      single <- !is.na(sky$edges[,1]) & !is.na(sky$edges[,2]) & sky$edges[,1]==sky$edges[,2]
      if (any(single)) graphics::points(a[single,1]+shift,a[single,2],
        pch=1,cex=.8,col=rep(line.ink,length.out=nrow(a))[single],lwd=.5)
    }
    s <- sky$stars
    graphics::points(s[,1]+shift,s[,2],cex=s[,3],pch=19,col=stars.ink)
    for (id in names(sky$bodies)) {
      p <- sky$bodies[[id]]
      if (id=="moon" && !is.null(sky$lunar)) {
        theta <- seq(0,2*pi,length.out=193)
        radius <- if (is.null(sky$diameters)) .25 else sky$diameters[["moon"]]/2
        x <- p[1]+shift
        y <- p[2]
        graphics::polygon(x+radius*cos(theta),y+radius*sin(theta),
                          col=if (paper) "black" else if (!is.null(sky$light)) .sky_colour(sky$light,p[1],p[2],opaque=TRUE) else bg,border=NA)
        lit <- .lunar_lit_polygon(sky$lunar$fraction,sky$lunar$angle,radius)
        graphics::polygon(x+lit$x,y+lit$y,col=if (paper) "white" else "ivory",border=NA)
        if (paper) graphics::lines(x+radius*cos(theta),y+radius*sin(theta),col="black",lwd=.5)
      } else if (id %in% c("sun","moon")) {
        theta <- seq(0,2*pi,length.out=61)
        radius <- if(is.null(sky$diameters)) .26 else sky$diameters[[id]]/2
        fill <- "white"
        if ((!paper || colour.paper) && id=="sun") {
          ramp <- grDevices::colorRamp(c("#FFAC60","#FFE2A9","#FFFAE8"))
          fill <- grDevices::rgb(ramp(pmax(0,pmin(1,p[2]/20))),maxColorValue=255)
        }
        graphics::polygon(p[1]+shift+radius*cos(theta),p[2]+radius*sin(theta),
                          col=fill,border=if(paper) "black" else NA,lwd=.5)
      } else {
        colours <- c(mercury="#DAD7D3",venus="#FFF9E8",mars="#FFAF88",jupiter="#F5DDC1",saturn="#EDD9A0")
        if (colour.paper) colours <- c(mercury="#525456",venus="#746337",mars="#A04429",jupiter="#836144",saturn="#89702D")
        size <- if (is.null(sky$planet.size)) 1 else sky$planet.size[[id]]
        colour <- if (paper && !colour.paper) ink else colours[[id]]
        if (length(sky$planet.opacity)) colour <- grDevices::adjustcolor(colour,alpha.f=sky$planet.opacity[[id]])
        graphics::points(p[1]+shift,p[2],pch=19,cex=size,col=colour)
      }
    }
  }
  if (isTRUE(sky$display[["zenith"]])) {
    graphics::points(mean(box[1:2]),90,pch=3,col=ink)
    graphics::text(mean(box[1:2]),90,"Zenith",pos=1,col=ink,cex=.7)
  }
}

# Haversine distance in degrees, including fields that cross the north seam.
.cluster_distance <- function(az,alt,x,y) {
  rad <- pi/180
  a <- sin((y-alt)*rad/2)^2 + cos(y*rad)*cos(alt*rad)*sin((x-az)*rad/2)^2
  2*asin(sqrt(pmin(1,pmax(0,a))))/rad
}

.draw_cluster_glows <- function(clusters,box) {
  if (!NROW(clusters)) return(invisible(NULL))
  for (i in which(clusters$alpha > 0)) {
    r <- clusters$radius[i]
    alt <- clusters$alt[i]
    lo <- max(box[3],alt-r,-90)
    hi <- min(box[4],alt+r,90)
    if (hi <= lo) next
    span <- if (abs(alt)+r >= 90) 180 else
      asin(min(1,sin(r*pi/180)/cos(alt*pi/180)))*180/pi
    for (turn in seq(floor((box[1]-clusters$az[i]-span)/360),
                     ceiling((box[2]-clusters$az[i]+span)/360))) {
      centre <- clusters$az[i]+360*turn
      left <- max(box[1],centre-span)
      right <- min(box[2],centre+span)
      if (right <= left) next
      x <- seq(left,right,length.out=81)
      y <- seq(hi,lo,length.out=81)
      distance <- outer(y,x,function(y,x) .cluster_distance(centre,alt,x,y))
      weight <- pmax(0,1-(distance/r)^2)^3
      ink <- grDevices::rgb(1,1,1,as.vector(weight)*clusters$alpha[i])
      graphics::rasterImage(matrix(ink,nrow=length(y)),left,lo,right,hi,interpolate=TRUE)
    }
  }
}
