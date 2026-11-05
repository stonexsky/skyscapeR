#' @description Analyse the directions of archaeological structures in relation to the local horizon and
#' to proposed celestial or terrestrial targets. The package supports reduction of field
#' measurements, coordinate transformation, visibility modelling and simulation-based
#' significance tests.
#' @details For a worked example from field measurements through to a multi-site test, open \code{vignette("skyscapeR", package = "skyscapeR")}.
#' Site coordinates are latitude and longitude in decimal degrees (north and east positive),
#' followed by elevation in metres where required. True bearings are clockwise from north.
#'
#' @section Recommended astronomical data:
#' Install \pkg{swephRdata} if you need astronomical positions and do not already have covering Swiss Ephemeris files.
#' Solar, lunar, planetary and fixed-star positions, and events derived from them,
#' require Swiss Ephemeris files covering the requested body and epoch. This applies
#' to modern and prehistoric dates. Install once, then restart R:
#' \preformatted{install.packages("swephRdata", repos = "https://rstub.r-universe.dev", type = "source")}
#' Alternatively, set \env{SE_EPHE_PATH} to existing Swiss Ephemeris data directories
#' before loading skyscapeR. The bundled star catalogue supplies stellar identifiers;
#' it cannot replace the planetary data used for stellar positions. Missing coverage
#' stops the calculation with installation instructions. Users can work with measured
#' horizons, calendar conversions and orientation statistics without these files.
#'
#' @section How do I reduce a field observation to a true bearing?:
#' Convert degrees, minutes and seconds with \code{\link{ten}}. For a compass survey, estimate magnetic declination with \code{\link{mag.dec}} and apply the correction with \code{\link{reduct.compass}}. For a solar sight, \code{\link{sunAz}} gives the solar bearing and \code{\link{reduct.theodolite}} reduces the instrument readings. Display the corrected bearings with \code{\link{plotAzimuth}}.
#'
#' @section How do I obtain, inspect and export a site horizon?:
#' \strong{Build or import.} Use \code{\link{createHor}} for survey bearings or a DEM,
#' \code{\link{importPanorama}} for a sky-masked photograph and \code{\link{importPeakFinder}}
#' for a PeakFinder landscape. For HeyWhatsThat, request a panorama with \code{\link{createHWT}}
#' or retrieve one by ID with \code{\link{downloadHWT}}.
#'
#' \strong{Inspect or share.} Add reconstructed structures with \code{\link{horFeature}},
#' plot the horizon, or read its altitude at a bearing with \code{\link{hor2alt}}.
#' Use \code{\link{exportPanorama}} for an equirectangular PNG or \code{\link{exportStellarium}}
#' for a landscape archive. To move an image-backed horizon, keep the PNG with the saved
#' object and reconnect it with \code{\link{relinkHor}}.
#'
#' @section Which declinations does a structure face, given measurement uncertainty?:
#' For a single sightline, use \code{\link{az2dec}} with the site and horizon altitude. For uncertain bearings, construct distributions with \code{\link{az.pdf}} and transform them with \code{\link{coordtrans}}. Use \code{\link{dec.pdf}} when the measurements are already in declination. Inspect individual curves with \code{\link[=plot.skyscapeR.pdf]{plot(pdf)}}, sum them with \code{\link{spd}} and display the result with \code{\link[=plot.skyscapeR.spd]{plot(summed)}}. Use \code{\link{hpdi}} to identify the highest-density intervals and \code{\link{plotBars}} to display orientation intervals against candidate declinations.
#'
#' @section Could the observed pattern arise under a null hypothesis?:
#' Use \code{\link{randomTest}} for uniformly random azimuths and \code{\link{modelTest}} for a specified distribution of azimuths or declinations. Inspect the global test and local intervals with \code{\link[=print.skyscapeR.sigTest]{print(fit)}} and \code{\link[=plot.skyscapeR.sigTest]{plot(fit)}}. For independent structures with a prespecified probability of a hit, \code{\link{bernoulli.trial}} gives an exact-count or upper-tail probability. \code{\link{pval2stars}} formats a probability for display.
#'
#' @section Which celestial targets fit the sightline and epoch?:
#' Screen a declination interval and year range with \code{\link{findTargets}}. For solar limits use \code{\link{jS}} and \code{\link{dS}} at the June and December solstices and \code{\link{eq}} at the equinox. \code{\link{spatial.equinox}} uses the local horizon; \code{\link{zenith}} and \code{\link{antizenith}} give the overhead and opposite solar limits. For northern and southern lunar standstill limits use \code{\link{nMjLX}} and \code{\link{sMjLX}} (major) or \code{\link{nmnLX}} and \code{\link{smnLX}} (minor). \code{\link{obliquity}} supplies the epoch-dependent tilt and \code{\link{parallax.corr}} the approximate site correction. \code{\link{planet.extremes}} estimates planetary limits. Resolve stellar names with \code{\link{star.names}} and obtain coordinates with \code{\link{star}}. Prepare targets for plotting with \code{\link{skyObject}}.
#'
#' @section When and where could an observer see a celestial event?:
#' Use \code{\link{body.position}} for coordinates at an instant, \code{\link{orbit}} for a path and \code{\link{riseset}} for crossings of the local horizon. \code{\link{solar.date}} finds dates when the Sun reaches a specified declination. Use \code{\link{moonphase}} for lunar illumination, \code{\link{EFM}} for equinoctial full moons and \code{\link{CC}} for solstitial crescent crossovers. Model seasonal stellar visibility with \code{\link{star.phases}} and display the seasons with \code{\link[=plot.skyscapeR.starphases]{plot(phases)}}. Use \code{\link{sketchSky}} or \code{\link{stellariumSky}} for a horizon plot, step through time with \code{\link{skyStep}} or use \code{\link{planisphere}} for a wider view. \code{\link{constellations}} supplies the constellation line segments.
#'
#' @section How do I specify historical dates and observing conditions?:
#' Inspect or set the session defaults with \code{\link{skyscapeR.vars}}. Assemble a timestamp with \code{\link{timestring}}, convert it to a UT Julian date with \code{\link{time2jd}} and convert back with \code{\link{jd2time}}. Use \code{\link{BC.AD}} to label an astronomical year and \code{\link{long.date}} to format a month and day. Year zero is 1 BC. Supply the intended calendar and an explicit UTC or fixed offset for dates before 1900.
#'
#' @section Which datasets can I use for a worked analysis?:
#' Pair \code{\link{malta_entrances}} with \code{\link{malta_horizons}} to connect temple entrance bearings with DEM-derived site horizons. To follow a single sightline and compare groups of Maltese temples, open \code{vignette("skyscapeR", package = "skyscapeR")}.
#'
#' @section Console output:
#' Set \code{verbose = FALSE} on functions that offer it to silence routine updates.
#' For a quiet session, set \code{options(skyscapeR.verbose = FALSE)}.
#' Warnings and errors remain visible. Progress bars normally appear in interactive
#' sessions; DEM tracing also accepts an explicit progress switch. Use
#' \code{print(fit)} for significance results in the console or a report.
#'
#' @section Recording an analysis:
#' Preserve original measurements, horizon profiles, feature layers and external PNG
#' files with the analysis script. Record the package and ephemeris versions, calendar,
#' timezone, refraction, uncertainty model, grid settings and simulation seed.
#' Refine numerical sampling and compare representative cases using matched conventions.
#' When interpreting a reconstruction, distinguish predicted events from visibility under local observing conditions.
#' See \code{\link{EFM}}, \code{\link{CC}} and \code{\link{star.phases}}
#' for the astronomical selection and empirical visibility assumptions.
#'
#' @section Updating older scripts:
#' The README migration guide identifies analyses and figures to rerun. Saved results
#' retain their original numbers until recalculated. \code{\link{skyObject}} is the
#' primary target constructor; \code{sky.objects} remains a compatibility alias.
#' The README migration guide documents the removal of retired functions, including sigTest() and curvigram().
#' @seealso \code{\link{skyscapeR.vars}}, \code{\link{createHor}}, \code{\link{az.pdf}}, \code{\link{dec.pdf}}, \code{\link{coordtrans}}, \code{\link{randomTest}}, \code{\link{modelTest}}
#' @author Fabio Silva (\email{fsilva@stonexsky.com}). Contributors: Billy Finnimore, Victor Reijs and Rita Gautschy.
"_PACKAGE"

## usethis namespace: start
#' @importFrom methods is
## usethis namespace: end
NULL
