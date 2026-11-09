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
#' the cutoff includes additional stars only where the catalogue supplies them.
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
#' @section Reading a natural sketch:
#' The sky is calculated for the specified instant. Stars, clusters, planets and
#' the Moon contribute light on a common visual (V-band) magnitude scale.
#' Diffuse Milky Way light and atmospheric background are measured per unit
#' sky area. Their light is added before brightness is compressed for display.
#' To inspect the values underlying a view, use \code{\link{skyRadiance}}.
#'
#' Stars remain unresolved as the field is enlarged. Mercury, Venus, Mars,
#' Jupiter, Saturn and the Moon use apparent disc sizes and illuminated
#' fractions from Swiss Ephemeris. A close view can therefore show the phase
#' of Venus. Star colours follow catalogue spectral classes; the diffuse band
#' has a neutral tint. Colour and projection influence the screen appearance.
#' Publication styles use magnitude-scaled star and planet symbols on white. In publication colour, stellar hue identifies the leading catalogue spectral class O through M; other or missing classes use neutral grey.
#'
#' Azimuth increases clockwise from true north. Refraction, pressure and
#' temperature come from \code{\link{skyscapeR.vars}} when the sky is calculated.
#' The Sun and Moon switches control their contributions to sky illumination as
#' well as their discs. Use \code{display = list(moon = FALSE)} for a moonless
#' comparison. The terrain masks the drawn sky; a local ridge hiding the Moon
#' can still leave scattered moonlight in the sky.
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
#' @section Twilight and extinction:
#' At twilight the sky brightens differently towards and away from the Sun.
#' Zenith brightness uses the Paranal V-band observations of Patat et al.
#' (2006) at solar depressions of 5 to 15 degrees. Sacramento Peak measurements
#' supply earlier twilight levels and the directional variation (Koomen et al.
#' 1952). These clear-sky estimates depend on site conditions. Moonlight uses
#' the scattering model of Krisciunas and Schaefer (1991), with its dependence
#' on lunar phase, altitude and angular separation.
#'
#' Direct light fades towards the horizon as its path through the atmosphere
#' lengthens. The coefficient \code{k} sets the loss per airmass. The model
#' uses Rozenberg's airmass approximation (Krisciunas and Schaefer 1991);
#' uncertainty is greatest close to the skyline. The twilight colours are
#' illustrative (Lee 1994, 2015).
#' @section Milky Way and clusters:
#' The Milky Way is traced from the diffuse light in Deep Star Maps 2020.
#' Within the smoothed envelope, retained dark patches have at least 10 percent
#' local contrast; two bright tiers begin at 1.10 and 1.75 times the mean band
#' tone after the dark patches are excluded. Retained patches cover at least
#' 10 square degrees. Natural sketches anchor the mean surface brightness
#' of a four-degree field in Carina to 20.7 V magnitudes per square arcsecond
#' (van der Kruit 2024). Brightness elsewhere is approximate. Atmospheric
#' extinction can make the band difficult to follow down to its calculated
#' horizon crossing. The Galactic equator is an optional geometric guide.
#'
#' A supported cluster contributes its catalogue light across its angular
#' extent. Light already assigned to plotted stars within that footprint is
#' subtracted before the remaining glow is drawn. This accounts for projected
#' overlap; it does not establish membership. Integrated cluster magnitude
#' alone is insufficient to predict first detection (Crumey 2014).
#' @section Using a sketch with a stellar phase date:
#' In 69 of 70 individual comparisons, the estimated first morning of visibility was within a week of the predicted heliacal-rising date.
#'
#' Use \code{\link{riseset}} to locate the skyline crossing on each morning,
#' then \code{\link{skyStep}} to follow the star as it climbs through a
#' brightening sky. A plotted pixel does not establish naked-eye visibility.
#' A quantitative test needs a point-source threshold and an assumed observer
#' limit (Crumey 2014). See also \code{\link{star.phases}}.
#' Record site, horizon, epoch, calendar, timezone, extinction
#' and display switches when reporting a reconstruction.
#' @return An object of classes \code{skyscapeR.sketch} and \code{skyscapeR.sky}.
#' It stores \code{jd}, \code{loc}, \code{style}, \code{display}, \code{settings}
#' (calendar, magnitude cutoff and extinction), and \code{atmosphere}.
#' \code{stars} and \code{coords} contain azimuth, altitude and symbol size;
#' \code{star.opacity} retains an earlier symbol-opacity diagnostic; natural
#' sketches draw from \code{light.budget}. \code{bodies}, \code{diameters},
#' \code{planet.size} and \code{lunar} describe selected Solar System bodies.
#' \code{clusters} holds supported cluster footprints and opacity when applicable. The natural view uses the integrated Pleiades magnitude as the light budget within its projected footprint. When the plotted stars exceed that budget, their fluxes are scaled together to preserve their relative brightness. Any residual light uses a broad, tapered profile;
#' \code{guides}, \code{light}, \code{background} and twilight fields support drawing.
#' Components can be empty when a layer is disabled. No image file is written.
#' @references Patat, F., Ugolnikov, O. S. and Postylyakov, O. V. (2006).
#' UBVRI twilight sky brightness at ESO-Paranal. Astronomy and Astrophysics 455, 385-393.
#' \doi{10.1051/0004-6361:20064992}.
#' @references Koomen, M. J. et al. (1952). Measurements of the brightness
#' of the twilight sky. JOSA 42, 353-356. \doi{10.1364/JOSA.42.000353}.
#' @references van der Kruit, P. C. (2024). Pannekoek's Galaxy.
#' Journal of Astronomical History and Heritage 27, 135-199.
#' \doi{10.3724/SP.J.1440-2807.2024.01.10}.
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
#' @section Validation:
#' W02 verified the reconstructed Gacrux bearing and altitude at the same instant and viewpoint using direct equatorial coordinates and spherical geometry. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w02-field-to-sky}{W02 check in the technical manual}.
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
  shared <- style == "natural"
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
    light$calibrated.twilight <- TRUE
  }
  body.flux <- numeric()
  planet.discs <- list()
  if (shared && flags["moon"]) {
    light$moon.mag <- pheno[5]
    body.flux["moon"] <- .source_flux(pheno[5], bodies$moon[2], k)
  }
  if (shared) background <- .sky_colour(light,0,90,opaque=TRUE)
  planet.opacity <- numeric()
  if (flags["planets"]) {
    ids <- c("mercury","venus","mars","jupiter","saturn")
    bodies[ids] <- lapply(ids,position)
    for (id in ids) {
      pos <- bodies[[id]]
      ph <- .sketch_pheno(jd,id,loc)
      mag <- ph[5]
      if (shared) {
        if (is.null(solar)) solar <- position("sun")
        diameters[id] <- ph[4]
        planet.discs[[id]] <- list(fraction=ph[2],angle=.lunar_light_angle(pos,solar))
        body.flux[id] <- .source_flux(mag,pos[2],k)
      }
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
    if (length(keep) && !is.null(catalogue$stellar)) {
      star.colours <- if (style=="publication-colour") {
        .publication_star_colours(catalogue$stellar[keep],stars)
      } else .sketch_star_colours(catalogue$stellar[keep],catalogue$magV[keep])
    }
  }
  if (flags["constellation.lines"]) edges <- do.call(rbind,constLines(selection))
  ids <- unique(as.character(c(stars,edges,clusters$identifier)))
  ids <- ids[!is.na(ids)]
  opacity <- stats::setNames(numeric(length(ids)),ids)
  star.flux <- opacity
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
      if (shared) star.flux[i] <- .source_flux(magnitude,pos[3],k)
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
    if (shared) {
      clusters$total.flux <- .source_flux(clusters$magV,clusters$alt,k)
      clusters$residual.flux <- clusters$total.flux
    }
    clusters$alpha <- numeric(nrow(clusters))
    points <- coords[stars,,drop=FALSE]
    for (i in which(is.finite(integrated) & (shared | integrated > 0))) {
      distance <- .cluster_distance(clusters$az[i],clusters$alt[i],points[,1],points[,2])
      inside <- which(distance <= clusters$radius[i] & is.finite(distance))
      if (shared) {
        members <- stars[inside]
        drawn <- sum(star.flux[members],na.rm=TRUE)
        if (clusters$identifier[i]=="M 45" && drawn>clusters$total.flux[i]) {
          # Treat the catalogue aperture magnitude as the total light in this view.
          star.flux[members] <- star.flux[members] * clusters$total.flux[i]/drawn
          drawn <- clusters$total.flux[i]
        }
        clusters$residual.flux[i] <- max(0,clusters$total.flux[i]-drawn)
      }
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
  outline <- flags["milky.way"] && (shared || visibility > 0)
  milky.way.frame <- NULL
  if (outline || flags["galactic.equator"]) {
    frame <- .galactic_frame(jd)
    if (outline) milky.way.frame <- .galactic_local_frame(jd,loc,frame)
    if (flags["galactic.equator"]) {
      equator <- .galactic_paths(jd,FALSE,TRUE,frame)$galactic.equator
      guides$galactic.equator <- t(apply(equator,1L,function(eq) {
        pos <- swephR::swe_azalt(jd,1L,loc[c(2,1,3)],atmosphere$atm,atmosphere$temp,c(eq,1))$xaz
        c((pos[1]+180)%%360,pos[2])
      }))
    }
  }
  .with_astronomy_provenance(structure(list(light.budget=if(shared) list(stars=star.flux[stars],bodies=body.flux,planet.discs=planet.discs) else NULL,clusters=clusters,light=light,planet.opacity=planet.opacity,star.opacity=opacity,milky.way.visibility=visibility,milky.way.frame=milky.way.frame,guides=guides,diameters=diameters,planet.size=planet.size,star.colours=star.colours,
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

.publication_star_colours <- function(type,ids) {
  # Keep spectral hues legible on white. Magnitude is encoded by symbol size.
  palette <- c(O="#294C9A",B="#416DB6",A="#637C9C",F="#847B60",
               G="#91691C",K="#B15B23",M="#AC4039")
  colours <- unname(palette[substr(trimws(type),1,1)])
  colours[is.na(colours)] <- "#61656A"
  stats::setNames(colours,ids)
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

# Moonless zenith in nL. Patat et al. (2006), V band, for Sun -15 to -5.
# Earlier twilight anchors are Sacramento Peak zenith photometry (Koomen
# et al. 1952, Table I); interpolate logarithmically to the adjacent regimes.
.twilight_zenith <- function(h) {
  patat <- function(d) {
    x <- d-5
    34.08*exp(20.7233-.92104*(11.84+1.518*x-.057*x^2))
  }
  if(h >= -15 && h <= -5) return(patat(-h))
  unit <- 10.76391*pi*1e5 # candles/ft^2 to nL
  10^stats::approx(c(-18,-15,-5,-3,0,5,15,90),
    log10(c(60,patat(15),patat(5),unit,8*unit,31*unit,3e8,1e9)),
    xout=h,rule=2)$y
}

# Measured angular gradients at Sacramento Peak. Interpolate log ratios so
# that the field stays positive and all azimuths converge at the zenith.
.twilight_ratio <- local({
  angles <- c(0,10,30,50,70,90)
  log_ratio <- function(x) log(x/x[,6])
  toward <- log_ratio(rbind(c(.9,.7,.18,.06,.033,.022),c(20,15,4.7,1.9,1.3,1),
                            c(150,150,37,16.5,9.5,8)))
  away <- log_ratio(rbind(c(.078,.068,.043,.035,.025,.022),c(2.4,2,2.1,1.6,1.1,1),
                          c(38,37,22,14,9.5,8)))
  function(h,delta,alt) {
    u <- min(2,max(0,(h+6)/3))
    i <- min(2,1+floor(u))
    mix <- u-(i-1)
    profile <- function(mat) stats::approx(angles,(1-mix)*mat[i,]+mix*mat[i+1,],
                                           xout=alt,rule=2)$y
    w <- (1+cos(delta*pi/180))/2
    exp(w*profile(toward)+(1-w)*profile(away))
  }
})

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
    if (isTRUE(light$calibrated.twilight)) {
      h <- light$sun[2]
      zenith <- .twilight_zenith(h)
      # Match the existing daytime gradient above civil twilight.
      measured <- .twilight_ratio(h,az-light$sun[1],alt)
      weight <- pmin(1,pmax(0,h/5))
      daytime <- exp(1.4*(cos(separation(light$sun)*rad)-sin(h*rad)))
      solar <- (zenith-60)*measured^(1-weight)*daytime^weight
    } else {
    zenith <- 10^stats::approx(c(-18,-12,-6,-.833,15,90),
      log10(c(60,240,1e5,3e6,3e8,1e9)),xout=light$sun[2],rule=2)$y
    solar <- (zenith-60)*exp(1.4*(cos(separation(light$sun)*rad)-sin(light$sun[2]*rad)))
    }

  }
  lunar <- 0
  if (!is.null(light$moon) && light$moon[2]>0 && light$k>0) {
    # Krisciunas & Schaefer (1991), equations 20-23. Bound the near-disc extrapolation.
    rho <- pmax(10,separation(light$moon))
    scatter <- 10^5.36*(1.06+cos(rho*rad)^2)+10^(6.15-rho/40)
    intensity <- 10^(-.4*(3.84+.026*light$phase+4e-9*light$phase^4))
    if (!is.null(light$moon.mag)) intensity <- 10^(-.4*(3.84 + light$moon.mag + 12.73))
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

.lunar_plot_path <- function(centre,east,up,shift) {
  rad <- pi/180
  h <- centre[2]*rad
  distance <- sqrt(east^2+up^2)*rad
  bearing <- atan2(east,up)
  alt <- asin(pmax(-1,pmin(1,sin(h)*cos(distance)+cos(h)*sin(distance)*cos(bearing))))
  offset <- atan2(sin(bearing)*sin(distance)*cos(h),cos(distance)-sin(h)*sin(alt))/rad
  # Follow the boundary continuously across the azimuth seam.
  az <- centre[1]+shift+offset[1]+c(0,cumsum((diff(offset)+180)%%360-180))
  list(x=az,y=alt/rad)
}

.lunar_plot_fill <- function(path,col) {
  x <- path$x
  y <- path$y
  if (abs(tail(x,1)-x[1])>180) {
    # A boundary enclosing a pole closes along the edge of the sky map.
    x <- c(x,tail(x,1),x[1])
    y <- c(y,rep(sign(mean(y))*90,2))
  }
  graphics::polygon(x,y,col=col,border=NA)
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
  if (colour.paper && NROW(sky$stars)) {
    ids <- rownames(sky$stars)
    if (!is.null(names(sky$star.colours)) && all(ids %in% names(sky$star.colours))) {
      stars.ink <- sky$star.colours[ids]
    } else {
      # Older sketches stored the natural, magnitude-faded tints.
      catalogue <- .read_star_catalogue()
      stars.ink <- .publication_star_colours(catalogue$stellar[match(ids,catalogue$identifier)],ids)
    }
  } else if (!paper && length(sky$star.colours)) stars.ink <- sky$star.colours
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
  shared <- !is.null(sky$light.budget)
  if (shared) {
    .draw_shared_light(sky,box)
  } else if (!is.null(sky$milky.way.frame)) {
    atmosphere <- sky$atmosphere
    atmosphere$jd <- sky$jd
    atmosphere$geopos <- sky$loc[c(2,1,3)]
    .draw_galactic_layers(sky$milky.way.frame,visibility,band.k,atmosphere,paper=paper)
  } else if (!is.null(sky$guides$south)) {
    # Retain the original geometry of sketches saved with the earlier outline.
    .draw_galactic_band(sky$guides$south,sky$guides$north,visibility=visibility,k=band.k)
  }
  if (!shared && !is.null(sky$light)) .sky_veil(sky$light,box)
  if (!shared) .draw_cluster_glows(sky$clusters,box)
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
    if (!shared) graphics::points(s[,1]+shift,s[,2],cex=s[,3],pch=19,col=stars.ink)
    for (id in names(sky$bodies)) {
      if (shared && id != "sun") next
      p <- sky$bodies[[id]]
      if (id=="moon" && !is.null(sky$lunar)) {
        theta <- seq(0,2*pi,length.out=193)
        radius <- if (is.null(sky$diameters)) .25 else sky$diameters[["moon"]]/2
        limb <- .lunar_plot_path(p,radius*cos(theta),radius*sin(theta),shift)
        dark <- if (paper) "black" else if (!is.null(sky$light)) .sky_colour(sky$light,p[1],p[2],opaque=TRUE) else bg
        .lunar_plot_fill(limb,dark)
        lit <- .lunar_lit_polygon(sky$lunar$fraction,sky$lunar$angle,radius)
        .lunar_plot_fill(.lunar_plot_path(p,lit$x,lit$y,shift),if (paper) "white" else "ivory")
        if (paper) graphics::lines(limb$x,limb$y,col="black",lwd=.5)
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
