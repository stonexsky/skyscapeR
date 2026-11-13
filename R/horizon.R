
#' Construct a site horizon from measurements or an elevation model
#' @description Build a continuous horizon profile from measurements of horizon altitude at known
#' bearings around a site, or predict the horizon and visible terrain from a DEM.
#' @details \strong{Measured bearings.} To evaluate sightlines between measured bearings, the constructor repeats the azimuths
#' one turn in each direction, removes duplicate azimuth entries and linearly interpolates
#' altitude and altitude uncertainty onto a 0.01-degree grid from -90 to 450 degrees.
#' Scalar uncertainty is repeated. Optional smoothing applies a running mean over the
#' specified number of grid samples to both series, with zero fill at the outer ends.
#' Latitude, longitude and elevation identify the observer; omitted elevation defaults to
#' zero. Interpolation cannot recover narrow features between measurements.
#'
#' \strong{Observer and DEM coverage.} For a DEM, supply a point with a CRS, projected x and y with loc.crs, or WGS84 latitude and longitude.
#' The optional third coordinate is ground elevation in metres in the DEM's vertical
#' datum. Otherwise ground elevation is bilinearly sampled at the site. Observer
#' elevation in the returned georeference is ground elevation plus eye.height.
#' A terrain model represents bare earth. A surface model includes the objects
#' represented in its elevations. Use a single continuous band in metres with a
#' defined horizontal CRS. Missing samples and terrain beyond the raster edge are skipped. A warning and the incomplete.bearings metadata identify incomplete coverage. Each bearing needs at least one usable sample.
#'
#' \strong{Sampling narrow ridges.} Terrain is sampled along true bearings from north, clockwise. A WGS84 local
#' azimuthal equidistant projection sets the distances and bearings. Sample coordinates
#' are transformed to the source grid for bilinear elevation extraction in batches.
#' Each horizon angle is the maximum sampled elevation angle along its ray. The
#' requested radius is included wherever covered. A missing observer ground elevation without an override or an observer outside the
#' DEM cause an error. Rays and distance samples are discrete. A finer constructor
#' output grid cannot recover a ridge between rays, steps or DEM cells. Increase
#' the radius and refine both sampling intervals to check convergence for a site.
#' At 5 km the default 0.25-degree rays are approximately 22 metres apart.
#' The supplied altitude uncertainty is an assumption in degrees.
#'
#' \strong{Curvature and terrestrial refraction.} By default the apparent altitude at each sampled distance d (metres) follows
#' Patat (2011), equations 1 and 2. The geometric angle is atan((z-observer)/d - d/(2R)).
#' Add (d/R) times 0.252 P (34.2 + tempgrad) / T squared in radians, using
#' temperature T in kelvins and the latitude-dependent radius R from Appendix A1.
#' The same corrected angles determine the numerical horizon and the visible terrain in the PNG.
#' An explicit numeric refraction.k selects a spherical vertical section of effective radius
#' 6371008.8 / (1 - refraction.k) metres. At ground distance d and elevation z,
#' the forward component is (R + z) sin(d/R) and the vertical component is
#' z minus observer elevation minus 2 (R + z) sin(d/(2R)) squared, where the
#' square applies to the sine. The angle is obtained with atan2. With curvature
#' disabled, use the height difference and horizontal distance directly.
#' Terrestrial refraction applies to the terrain
#' sightline. Atmospheric refraction is controlled separately in astronomical calls.
#'
#' \strong{Panorama and reconstructed features.} By default a second set of rays paints the first visible terrain at each
#' azimuth and elevation, using local surface normals, elevation colours and
#' distance haze. The transparent sky edge follows the numerical horizon.
#' With features, each PNG column is traced directly through the reconstructed
#' terrain; its edge can differ from the interpolated numerical horizon between
#' rays. Refine az.step and width when examining narrow structures.
#' The panorama spans the full 360 degrees from nadir to zenith. The horizon object
#' stores its bounds, checksum, attribution and licence in \code{metadata$image}.
#' Save the external PNG alongside the RDS. If you move them, use \code{\link{relinkHor}}
#' to reconnect the image. A rendering failure returns the numerical horizon with a warning.
#' Plotting also falls back to the numerical horizon if the image is missing or has changed.
#' @param terrain.palette Named land/water pair from \code{\link{demPalettes}}.
#' Options are temperate (default), standalone, arid, alpine, coastal and relief.
#' @param land_cols,sea_cols Opaque R colour names or hex colours in increasing
#' elevation order. NULL inherits the respective ramp from terrain.palette.
#' Supply either vector to override just that ramp. One colour gives a uniform base. Colours interpolate in Lab space.
#' @param sea_level DEM elevation in metres at or below which the sea colour
#' ramp is used. Default 0. A level water surface uses the shallow colour. This is an elevation threshold. Use NULL to apply
#' the land palette over the whole visible elevation range, including depressions.
#' The temperate land ramp starts at the lowest visible land level or sea threshold, whichever is higher. Other palettes with a threshold start at sea_level. The ramp ends at the highest
#' visible elevation. The sea ramp starts at the lowest visible elevation.
#' @param palette_n Total number of palette entries, from 4 to 65536. Default
#' 256.
#' @param n_sea Number of entries reserved for the sea ramp, leaving at least
#' two land colours. Default 3.
#' @param light.azimuth,light.altitude Direction towards the illustrative light
#' source in degrees clockwise from true north and above the horizon. Defaults
#' are 315 and 45 degrees. Hillshading uses upward surface normals in metric
#' east, north and vertical coordinates. This is local illumination without
#' cast shadows from other terrain.
#' @param lighting Relief lighting (default) blends a light at the observer (75 percent)
#' with the directional source (25 percent). It shades both sides of a slope and
#' lifts deep shadows. The standalone renderer brightness response is
#' min(1, ambient + contrast * illumination), with illumination floored at 0.15.
#' Terrain and visible horFeature surfaces use the same lighting. Choose "directional"
#' for one-sided hillshade. Lighting affects the PNG only.
#' @param ambient Fraction of illumination retained on slopes facing away from
#' the light (0 to 1). Default 0.40. Use 1 for unshaded palette colours.
#' @param contrast Relief-light contrast (0 to 1), default 0.70, following the standalone renderer. It controls slope detail independently of ambient brightness. Ignored for directional lighting, which retains ambient + (1 - ambient) * illumination. Regenerate existing panoramas to apply new settings.
#' @param haze Maximum fraction of distance haze (0 to 1). Default 0.65.
#' Set haze=0 and ambient=1 to inspect the colour mapping alone.
#' @param features Optional \code{\link{horFeature}} object or list of such
#' objects describing reconstructed or highlighted structures. Requires dem.
#' Elevations are additive, with the largest positive and deepest negative
#' contributions combined at overlaps. Later footprints take colour precedence.
#' Features affect both numerical horizon and visible terrain. Their geometry
#' is kept out of the returned horizon; save the feature objects separately.
#' @param dem Optional single-band terra SpatRaster or local raster filename.
#' Supply dem in place of az and alt. Requires the optional terra package.
#' @param loc.crs CRS for numeric DEM observer coordinates in x,y order (easting, northing or longitude, latitude). Without it, numeric loc uses latitude, longitude in WGS84. A single terra SpatVector point carries its own CRS. An optional third numeric coordinate is ground elevation in metres.
#' @param ncores Maximum DEM tracing workers. Defaults to half the available logical cores, rounded down with a minimum of one. The startup message reports used and available cores. Workers share file-backed DEM data and process separate bearings.
#' @param verbose Announce the DEM radius, sampling and rendering stages.
#' @param progress Show a progress bar for each tracing stage. Default on in interactive verbose sessions.
#' @param max.distance DEM search radius in metres, default 10000. Include any
#' distant terrain that could form the horizon within the DEM coverage.
#' @param eye.height Observer height above the selected ground elevation in metres.
#' @param az.step Maximum spacing of numerical DEM rays in degrees. The default
#' is 0.25. An equal number of steps fills the circle.
#' @param panorama.distance.step The default, "adaptive", starts at one local DEM cell and increases radial
#' spacing with distance using the larger of the requested angular acuity and
#' one horizontal PNG pixel. This transverse footprint is a proxy for radial
#' detail. A narrow ridge may fall between samples. Use a positive number for
#' fixed spacing in metres or NULL to reuse distance.step. The numerical horizon
#' uses distance.step independently. Inspect metadata$image$sampling for the
#' effective angle, sample count and actual step range. On large DEM-only panoramas, the renderer checks near terrain against a conservative upper bound for the remaining terrain on each ray. It skips distant samples only where all remaining elevations are hidden and keeps neighbours for surface shading. Terrain with added archaeological features uses the full sampling path.
#' @param panorama.acuity Angular detail target in arcminutes for adaptive PNG sampling. Defaults to 1 (nominal 20/20 acuity). Larger values coarsen distant terrain.
#' @param distance.step Maximum spacing of samples along each ray in metres.
#' By default use half the smaller native cell dimension estimated at the site.
#' @param curvature Apply Earth curvature to DEM sightlines. Default TRUE.
#' @param refraction.k Default NULL applies Patat when curvature is enabled. Set zero for geometric curvature without terrestrial refraction, or supply a coefficient below one for the effective-radius model.
#' @param atm Ground atmospheric pressure in millibars for Patat terrestrial refraction.
#' @param temp Ground air temperature in degrees Celsius.
#' @param tempgrad Vertical temperature gradient in kelvins per kilometre (default -10).
#' @param panorama Render an external terrain PNG for a DEM. Default TRUE.
#' Set FALSE for the numerical horizon alone.
#' @param image.dir Directory for the DEM PNG. Defaults to the working directory. Site names form the filename; different images with the same name receive a checksum suffix.
#' @param width DEM panorama width in pixels, an even integer from 32 to 81920.
#' Default 4096 by 2048 pixels. Height is half the width. Widths above 8192 are assembled from temporary image tiles. Both dimensions grow with width, so a tenfold increase generates 100 times as many pixels.
#' @param attribution Credit for the DEM, included in plots and exports. The generic default credit is omitted from plots.
#' @param license Source licence or terms for the DEM, included in exports.
#' @return A skyscapeR.horizon list with metadata (site name and georeference) and a data frame named data. Its az, alt and alt.unc columns contain azimuth, horizon altitude and altitude uncertainty, all in degrees.
#' @references Patat, F. (2011). Horizon synthesis for archaeo-astronomical purposes. Astronomische Nachrichten 332, 743-749. \doi{10.1002/asna.201111570}.
#' @param az Array of azimuth values
#' @param alt Array of altitude values.
#' @param alt.unc Either a single value or an array of altitude
#' uncertainty.
#' @param loc Observer location. For DEMs, pass a single CRS-bearing terra SpatVector point directly as loc, or a numeric latitude and longitude. Ground elevation is read from the DEM unless supplied as a third numeric value. Measured horizons take latitude, longitude and elevation.
#' @param name Name of site.
#' @param smooth Logical to control whether to smooth horizon profile using rolling mean.
#' Defaults to FALSE
#' @param .scale Number of samples in the rolling mean when smooth = TRUE. The grid spacing is 0.01 degrees, so 1000 samples span approximately 10 degrees.
#' @seealso \code{\link{importPeakFinder}}, \code{\link{plot.skyscapeR.horizon}}
#' @import zoo
#' @export
#' @seealso \code{\link{createHWT}}, \code{\link{downloadHWT}}
#' @examples
#' \dontrun{
#' observer <- terra::vect("observer.gpkg") # one point with a CRS
#' h <- createHor(dem="terrain.tif", loc=observer, eye.height=1.6)
#' }
#' # Create a skyscapeR.horizon from 5 measurements:
#' az <- c(0,90,180,270,360)
#' alt <- c(0,5,5,0,0)
#' hor <- createHor(az, alt, alt.unc=0.1, loc=c(40.1,-8), name='Test')
#' plot(hor)
#' \dontrun{
#' # loc is latitude, longitude. Distances and elevations are metres.
#' terrain <- createHor(dem="site-dtm.tif", loc=c(51.19,-1.75),
#'   name="Survey site", max.distance=5000, eye.height=1.6,
#'   image.dir="horizons", attribution="Survey DTM", license="Survey terms")
#' plot(terrain)
#' exportStellarium(terrain, name="Survey site")
#' }
#' @section Accuracy, uncertainty and reproducibility:
#' H01 compared 334 bearings at 31 observing setups.
#' At the surveyed bearings, \code{createHor()} reproduced the observed terrain skyline to within half a degree in 91.9\% of cases.
#' See the \href{https://stonexsky.github.io/skyscapeR/manual/#h01-field-validation}{H01 field validation in the skyscapeR technical manual} for the survey protocol, station results and limitations.
#'
#' With the same 1 m terrain, \code{createHor()} reproduced Horizon’s skylines to a mean absolute difference of 0.004 degrees at each of two English sites. Against field measurements at 344 Irish bearings, \code{createHor()} using SRTM had a mean absolute error of 0.328 degrees, lower than either HeyWhatsThat or PeakFinder on that sample. Field accuracy varied by site. See the \href{https://stonexsky.github.io/skyscapeR/manual/#h02-skyline-comparison}{H02 comparison in the skyscapeR technical manual}.
#'
#' The numerical skyline is a maximum over sampled terrain, not an accuracy
#' certificate for the DEM. Horizontal and vertical datum differences, observer
#' height, missing raster values, vegetation and buildings can exceed the sampling
#' error. DEM elevations and a supplied observer ground elevation must use a
#' compatible vertical datum. No geoid conversion or vertical-datum correction is
#' performed automatically. Refine the angular and distance sampling and expand
#' the radius until the relevant sightlines are stable. Compare with surveyed
#' bearings and altitudes where available.
#'
#' \code{alt.unc} is supplied by the user and is not estimated from DEM resolution,
#' control residuals or atmospheric variability. Smoothing a measured profile can
#' remove narrow notches or peaks. Adaptive panorama sampling and illustrative
#' lighting concern the image; they do not increase numerical horizon precision.
#' The terrestrial refraction approximation assumes the supplied atmosphere and
#' cannot reproduce all near-ground gradients or mirage conditions.
#'
#' Save the DEM provenance, observer location and datum, constructor arguments,
#' feature objects, returned horizon and external PNG. Rebuilding from the same
#' source is distinct from merely replotting the saved horizon.
createHor <- function(az, alt, loc, dem=NULL, loc.crs=NULL, name='', alt.unc=0.5, smooth=F,
                      .scale=1000, eye.height=1.6, max.distance=10000, az.step=.25,
                      distance.step=NULL, features=NULL, curvature=TRUE, refraction.k=NULL,
                      atm=skyscapeR.env$atm, temp=skyscapeR.env$temp, tempgrad=-10, panorama=TRUE,
                      width=4096, panorama.distance.step="adaptive", panorama.acuity=1,
                      terrain.palette="temperate", land_cols=NULL, sea_cols=NULL, sea_level=0,
                      n_sea=3L, palette_n=256L, lighting=c("relief","directional"),
                      light.azimuth=315, light.altitude=45, ambient=.40, haze=.65,
                      image.dir=getwd(), attribution="User-supplied DEM", license="See source DEM terms", ncores=NULL, verbose=TRUE, progress=interactive() && verbose, contrast=.70) {
  if (!is.null(dem)) {
    if (!missing(az) || !missing(alt)) stop("Use either az/alt measurements or dem.")
    if (!identical(smooth,FALSE)) stop("DEM horizons use the sampled terrain profile. Set smooth=FALSE.")
    style <- if (isTRUE(panorama)) .dem_style(terrain.palette,land_cols,sea_cols,sea_level,n_sea,
                                             light.azimuth,light.altitude,ambient,haze,palette_n,contrast)
    if (isTRUE(panorama)) style$lighting <- match.arg(lighting)
    return(.dem_horizon(dem,loc,name,alt.unc,max.distance,eye.height,az.step,distance.step,
                        curvature,refraction.k,panorama,image.dir,width,attribution,license,style,features,loc.crs,ncores,verbose,progress,atm,temp,tempgrad,panorama.distance.step,panorama.acuity))
  }
  if (!is.null(features)) stop("features requires a DEM.")
  # return result
  hor <- c()
  hor$metadata$name <- name
  if (length(loc)<3) { .sky_message('No site elevation given. Assuming 0 metres above sea level.'); loc <- c(loc,0) }
  if (length(loc)<2) { stop('Site latitude, longitude and elevation needed.') }
  hor$metadata$georef <- loc; names(hor$metadata$georef) <- c('Lat','Lon','Elev'); dim(hor$metadata$georef) <- c(1,3)
  if (length(alt.unc) > 1 & length(alt.unc) < length(az)) { stop('Altitude uncertainty must be either a single value or have the same length as "az"') }
  if (length(alt.unc) == 1 ) { alt.unc <- rep(alt.unc, length(az)) }

  az <- c(az-360, az, az+360); alt <- rep(alt, 3); alt.unc <- rep(alt.unc, 3)
  ind <- which(duplicated(az)); if (length(ind)>0) { az <- az[-ind]; alt <- alt[-ind]; alt.unc <- alt.unc[-ind] }
  xx <- seq(0-90, 360+90, 0.01)
  alt <- approx(az, alt, xx)$y
  alt.unc <- approx(az, alt.unc, xx)$y
  hor$data <- data.frame(az = xx, alt = alt, alt.unc = alt.unc)

  if (smooth) {
    yy <- zoo::rollmean(alt, .scale, fill=0)
    hor$data$alt <- yy

    yy <- zoo::rollmean(alt.unc, .scale, fill=0)
    hor$data$alt.unc <- yy
  }
  class(hor) <- "skyscapeR.horizon"
  return(hor)
}



#' Export a site horizon as a Stellarium landscape
#' @description Reproduce a measured or modelled archaeological horizon as a Stellarium landscape for
#' visual comparison.
#' @details A landscape lets the same terrain obstruction be used in a planetarium reconstruction.
#' The export pairs the profile's azimuth and altitude samples in a polygonal horizon file
#' and records the site, appearance and attribution in a Stellarium configuration. The
#' export describes the horizon geometry; the altitude-uncertainty series remains in the R
#' object.
#'
#' DEM and SVG panoramas use a spherical landscape with the PNG rotated -90
#' degrees for Stellarium's east-origin mapping. A new PeakFinder import also
#' records the original ZIP location. If native tiles are requested and that archive still matches its checksum
#' and contains the supported calibrated eight-tile landscape, Stellarium uses
#' those original textures and their angular settings instead. The export also
#' includes the skyscapeR numerical profile and SVG-derived PNG.
#' exportStellarium is an alias for exportHor.
#'
#' The ZIP is saved as name-horizon.zip in the calling working directory.
#' Choose a new archive name each time. Names cannot contain path separators or line breaks.
#' Export requires a working system ZIP command.
#' @param peakfinder.layout Single spherical PNG (default) or original native tiles. Textured DEM and PeakFinder exports use their PNG silhouette.
#' @param max.texture.width Maximum exported DEM or photographic panorama width in pixels. The default is 4096. Use max.texture.width=2048 for a smaller texture, or Inf to keep the original dimensions. Stellarium supports textures only up to the graphics hardware limit. PNG compression does not reduce the decoded texture dimensions. Larger exports are at your own risk and may fail to load or crash Stellarium. This setting does not resize PeakFinder textures.
#' @param use.image Include a valid external panorama PNG when available. Default TRUE.
#' Missing or changed files warn and use a polygonal horizon instead.
#' @param hor Horizon data in \emph{skyscapeR.horizon} format.
#' @param name Horizon name to be displayed in \emph{Stellarium}, if different
#' from one in \emph{skyscapeR.horizon} object.
#' @param author Author, to be included in \emph{landscape.ini} file.
#' @param description Description, to be included in \emph{landscape.ini} file.
#' @param ground_col Color of ground. Defaults to \emph{Stellarium}'s default.
#' @param hor_col Colour of a polygonal horizon line. Textured DEM and PeakFinder landscapes omit this line.
#' @return Invisibly, the absolute path of the created ZIP archive.
#' @seealso \code{\link{createHor}}, \code{\link{downloadHWT}}, \code{\link{plot.skyscapeR.horizon}}
#' @references \href{https://stellarium.org/guide/}{Stellarium User Guide}, chapter on landscapes and texture sizes.
#' @section Validation:
#' H07 checked that a partial photograph is rejected for both image and numerical landscape export. See the \href{https://stonexsky.github.io/skyscapeR/manual/#h07-panorama-coverage}{H07 check in the technical manual}.
#' @export
#' @import utils
#' @examples
#' # Export a measured horizon from the example dataset.
#' \dontrun{
#' hor <- malta_horizons[[1]]
#' exportStellarium(hor, name = "Malta survey")
#' # For a DEM panorama, reduce the texture width on a memory-limited machine.
#' terrain <- createHor(dem = "terrain.tif", loc = c(51.19, -1.75))
#' exportStellarium(terrain, name = "Survey terrain", max.texture.width = 2048)
#' }
exportHor = function(hor, name, author="skyscapeR", description, use.image=TRUE,
                      peakfinder.layout=c("single","native"), max.texture.width=4096, ground_col,
                      hor_col) {
  if (!is.numeric(max.texture.width) || length(max.texture.width)!=1L || is.na(max.texture.width) || max.texture.width<2 || (is.finite(max.texture.width) && max.texture.width!=floor(max.texture.width))) stop("max.texture.width must be an integer of at least 2 or Inf.")
  peakfinder.layout <- match.arg(peakfinder.layout)
  single <- identical(hor$metadata$source$provider,"PeakFinder") && peakfinder.layout=="single" && use.image
  if (!is(hor,'skyscapeR.horizon')) { stop('No skyscapeR.horizon object found.') }

  if (missing(name)) { name <- hor$metadata$name }
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(trimws(name)) || grepl("[/\\\\\r\n]", name)) {
    stop('Provide a nonempty horizon name without path separators or line breaks.')
  }
  destination <- file.path(getwd(), paste0(name, '-horizon.zip'))
  if (file.exists(destination)) stop('Export file already exists: ', destination)
  work <- tempfile('skyscapeR-export-')
  if (!dir.create(work)) stop('Could not create export temporary directory.')
  original.wd <- getwd()
  on.exit({ setwd(original.wd); unlink(work, recursive=TRUE) }, add=TRUE)
  setwd(work)
  if (missing(description)) { description <- paste0("Horizon created using skyscapeR ", packageVersion('skyscapeR'), ".") }
  if (substr(description,nchar(description),nchar(description)) != ".") { description <- paste0(description,". Horizon created using skyscapeR ", packageVersion('skyscapeR'), ".") }
  if (missing(ground_col)) {ground_col = ".15,.45,.15"}
  if (missing(hor_col)) hor_col <- ".75,.45,.45"

  image <- if (use.image) .horizon_image(hor) else NULL

  photo <- identical(hor$metadata$source$provider,"Panorama")
  if (photo && diff(hor$metadata$source$coverage)<360-1e-8)
    stop("Stellarium landscape export requires full 360-degree coverage.")
  # Horizon data
  data <- data.frame(x = hor$data$az, y = hor$data$alt)
  # Prefer the central turn when a stored horizon includes interpolation copies.
  central <- data$x >= 0 & data$x < 360
  data <- data[order(!central), , drop=FALSE]
  data$x <- round(data$x %% 360,8) %% 360
  data <- data[!duplicated(data$x), , drop=FALSE]
  data <- data[order(data$x), , drop=FALSE]
  write.table(sprintf("%.8f %.8f", data$x, data$y), file="horizon.txt", quote=FALSE, sep = " ", row.names=F, col.names=F, fileEncoding = "UTF-8")

  # Landscape.ini file
  string.text <- c("[landscape]", paste0("name = ",name), "type = polygonal", paste0("author = ",author),
                   paste0("description = ", description), "polygonal_horizon_list = horizon.txt",
                   "polygonal_angle_rotatez = 0.00001", paste0("ground_color = ",ground_col),
                   paste0("horizon_line_color =  ",hor_col), "",
                   "[location]", "planet = Earth", paste0("latitude = ",hor$metadata$georef[1],"d"),
                   paste0("longitude = ",hor$metadata$georef[2],"d"), paste0("altitude = ",hor$metadata$georef[3]))
  if (photo && !is.finite(hor$metadata$georef[3])) string.text <- string.text[!grepl("^altitude =",string.text)]
  files <- c("landscape.ini", "horizon.txt")
  if (!is.null(image)) {
    if (photo && !is.null(image$source.rows)) {
      png::writePNG(.read_panorama_pixels(image),"calibrated.png")
      image$path <- "calibrated.png"
    }
    if ((photo || identical(hor$metadata$source$provider,"DEM")) && image$width>max.texture.width) {
      grid <- suppressWarnings(terra::rast(image$path))
      height <- max(1,floor(nrow(grid)*max.texture.width/ncol(grid)))
      reduced <- terra::spatSample(grid,max.texture.width*height,method="regular",as.raster=TRUE)
      terra::writeRaster(reduced,"panorama.png",datatype="INT1U",NAflag=NA)
    } else if (!file.copy(image$path,"panorama.png")) stop("Could not copy panorama PNG.")
    string.text[string.text=="type = polygonal"] <- "type = spherical"
    string.text <- append(string.text,c("maptex = panorama.png",
      paste0("maptex_top = ",image$altitude[2]),paste0("maptex_bottom = ",image$altitude[1]),
      paste0("angle_rotatez = ",image$azimuth[1]-90), "minimal_altitude = -90",
      paste0("license = ",image$license)),after=2L)
    i <- grep("^author =",string.text)
    string.text[i] <- paste0(string.text[i],"; ",image$attribution)
    # The calibrated texture edge defines the visible terrain silhouette.
    # Keep the numerical profile in the archive for independent use.
    if (isTRUE(hor$metadata$source$provider %in% c("DEM","PeakFinder","Panorama")))
      string.text <- string.text[!grepl("^polygonal_",string.text)]
    files <- c(files,"panorama.png")
  }
  # Preserve the original eight calibrated textures if the import archive is available.
  source <- hor$metadata$source
  if (!single && !is.null(image) && identical(source$provider,"PeakFinder") &&
      length(source$archive)==1L && file.exists(source$archive) &&
      identical(unname(tools::md5sum(source$archive)),source$zip.md5)) {
    members <- utils::unzip(source$archive,list=TRUE)$Name
    ini <- members[basename(members)=="landscape.ini"]
    assets <- c(sprintf("peakfinder%d.png",0:7),"peakfinderbg.png")
    prefix <- if (dirname(ini)==".") "" else paste0(dirname(ini),"/")
    cfg <- local({ con <- unz(source$archive,ini); on.exit(close(con)); readLines(con,warn=FALSE) })
    if (all(paste0(prefix,assets) %in% members) &&
        any(trimws(cfg)=="type = old_style") && any(trimws(cfg)=="calibrated = true")) {
      utils::unzip(source$archive,files=paste0(prefix,assets),junkpaths=TRUE)
      geom <- grep("^(decor_|nbside|nbsidetex|nb_decor_repeat|side[0-7] |tex[0-7] |groundtex |ground |ground_angle_|calibrated |tan_mode )",trimws(cfg),value=TRUE)
      string.text <- string.text[!grepl("^(maptex|angle_rotatez|minimal_altitude|polygonal_)",string.text)]
      string.text[string.text=="type = spherical"] <- "type = old_style"
      string.text <- append(string.text,geom,after=2L)
      files <- c(files,assets)
    }
  }
  writeLines(enc2utf8(string.text), "landscape.ini", useBytes=TRUE)

  # Save as zip file
  status <- utils::zip(zipfile='export.zip', files=files, flags='-q -r9X')
  if (!identical(as.integer(status), 0L) || !file.exists('export.zip')) {
    stop('Could not create horizon ZIP archive.')
  }
  members <- utils::unzip('export.zip', list=TRUE)$Name
  if (!setequal(members, files)) {
    stop('Horizon ZIP archive has unexpected contents.')
  }
  if (!file.copy('export.zip', destination, overwrite=FALSE)) {
    stop('Could not save horizon ZIP archive: ', destination)
  }
  invisible(destination)
}



.hwt_text <- function(path, ..., timeout=30, retry=FALSE) {
  if (!requireNamespace("httr2",quietly=TRUE))
    stop("Install httr2 to use the HeyWhatsThat service.")
  req <- httr2::request(paste0("https://www.heywhatsthat.com/",path))
  req <- httr2::req_url_query(req,...)
  req <- httr2::req_timeout(req,timeout)
  if (retry) req <- httr2::req_retry(req,max_tries=3,max_seconds=60)
  if (path=="api/query") req <- httr2::req_options(req,followlocation=FALSE)
  resp <- httr2::req_perform(req)
  if (httr2::resp_status(resp) >= 300L) stop("Unexpected HeyWhatsThat HTTP status: ",httr2::resp_status(resp))
  status <- httr2::resp_header(resp,"x-app-status")
  if (!is.null(status) && !grepl("^200( |$)",status))
    stop("HeyWhatsThat: ",status)
  httr2::resp_body_string(resp)
}

.hwt_capture <- function(html,pattern,field) {
  match <- regmatches(html,regexec(pattern,html,perl=TRUE))[[1]]
  if (length(match)<2L) stop("Missing or invalid HeyWhatsThat ",field,".")
  match[-1L]
}

#' Download a horizon from a HeyWhatsThat panorama
#' @description Recover an existing terrain-derived horizon for use in site-specific orientation
#' calculations.
#' @details A panorama describes the horizon seen from a particular observer location. The
#' optional httr2 package retrieves site metadata and 0.125-degree horizon bins over
#' HTTPS. Site elevation above sea level is reported in metres, converting feet
#' when necessary. Validated profiles are cached by panorama identifier for the R
#' session. Read-only downloads allow up to three attempts for transient HTTP
#' failures, with a 30-second limit per request. For each bin, the distance and
#' angular altitude determine a height difference. The mean angular change for
#' height offsets of plus and minus 9.73 metres estimates altitude uncertainty.
#' @return A skyscapeR.horizon list containing site metadata and horizon samples. The data columns az, alt and alt.unc are in degrees. Retrieval requires access to the HeyWhatsThat service.
#' @param HWTID This is the 8 character ID attributed by
#' \emph{HeyWhatsThat.com}
#' @export
#' @import utils
#' @references \href{https://www.heywhatsthat.com/}{HeyWhatsThat.com}
#' @seealso \code{\link{importPeakFinder}}, \code{\link{createHWT}}
#' @examples
#' \dontrun{
#' # Retrieve the Liverpool Cathedral horizon by its panorama ID.
#' hor <- downloadHWT('HIFVTBGK')
#' }
downloadHWT <- function(HWTID) {
  if (!is.character(HWTID) || length(HWTID)!=1L || is.na(HWTID) ||
      !grepl("^[A-Za-z0-9]{8}$",HWTID)) stop("Incorrect HeyWhatsThat ID.")
  html <- .hwt_text("iphone/pan.cgi",id=HWTID,retry=TRUE)
  visible <- gsub("(?is)<(script|style)\\b[^>]*>.*?</\\1>", " ", html, perl=TRUE)
  visible <- trimws(gsub("[[:space:]]+", " ", gsub("<[^>]+>", " ", visible)))
  if (grepl("panorama not found",visible,ignore.case=TRUE))
    stop("HeyWhatsThat reports 'Panorama not found' for ",HWTID,
         ". Check the ID and whether the panorama is available.",call.=FALSE)
  visible <- gsub("(?i)&(?:nbsp|#0*160|#x0*a0);|\\x{00a0}", " ", visible, perl=TRUE)
  visible <- gsub("(?i)&(?:deg|#0*176|#x0*b0);", "\u00b0", visible, perl=TRUE)
  coord <- function(label,hemispheres) {
    pattern <- paste0("(?i)\\b",label,"\\s*:?\\s*([0-9]+(?:\\.[0-9]+)?)\\s*",
                      "\\x{00b0}\\s*([",hemispheres,"])\\b")
    value <- .hwt_capture(visible,pattern,paste(label,"for",HWTID))
    as.numeric(value[1L]) * if (toupper(value[2L]) %in% c("S","W")) -1 else 1
  }
  lat <- coord("Latitude","NS")
  lon <- coord("Longitude","EW")
  elev <- .hwt_capture(visible,"(?i)\\(\\s*([-+]?[0-9]+(?:\\.[0-9]+)?)\\s*(ft|m)\\s+above sea level",paste("elevation for",HWTID))
  elevation <- as.numeric(elev[1L]) * if (tolower(elev[2L])=="ft") 0.3048 else 1
  name <- .hwt_capture(html,"(?is)<div\\b[^>]*\\bid\\s*=\\s*['\"]pan_top_title['\"][^>]*>(.*?)</div>",paste("site name for",HWTID))
  name <- trimws(gsub("<[^>]+>","",name))
  if (any(!is.finite(c(lat,lon,elevation))) || abs(lat)>90 || abs(lon)>180)
    stop("Invalid HeyWhatsThat site coordinates.")
  cache <- file.path(tempdir(),paste0("skyscapeR-HWT-",HWTID,".rds"))
  if (file.exists(cache)) {
    profile <- readRDS(cache)
  } else {
    csv <- .hwt_text("api/horizon.csv",id=HWTID,resolution=.125,retry=TRUE)
    samples <- utils::read.csv(text=csv)
    fields <- c("bin.bottom","altitude","distance..m.")
    if (!all(fields %in% names(samples)) || !nrow(samples) ||
        !all(vapply(samples[fields],is.numeric,logical(1))) ||
        any(!is.finite(as.matrix(samples[fields]))) || any(samples$distance..m.<=0))
      stop("Invalid HeyWhatsThat horizon CSV.")
    distance <- samples$distance..m.
    height <- distance*tan(samples$altitude*pi/180)
    error <- (atan((height+9.73)/distance)-atan((height-9.73)/distance))*90/pi
    profile <- data.frame(az=samples$bin.bottom,alt=samples$altitude,alt.unc=error)
    tmp <- tempfile(tmpdir=tempdir())
    on.exit(unlink(tmp),add=TRUE)
    saveRDS(profile,tmp)
    file.rename(tmp,cache)
  }
  site <- matrix(c(lat,lon,elevation),nrow=1L)
  names(site) <- c("Lat","Lon","Elev")
  structure(list(metadata=list(ID=HWTID,name=name,georef=site,elevation=elevation),
                 data=profile),class="skyscapeR.horizon")
}


#' Request a HeyWhatsThat panorama and retrieve its horizon
#' @description Request a terrain-derived horizon for an observing site when comparing orientations with
#' possible rising or setting directions.
#' @details A digital horizon supplies estimated angular obstruction around the site. The
#' optional httr2 package submits one HTTPS request for a panorama. Once the
#' service returns its identifier, the function waits two minutes before the first readiness check, then polls at intervals
#' of up to 15 seconds. The polling deadline defaults to 600 seconds. If the wait
#' expires, use the identifier with downloadHWT once the panorama is ready.
#' Individual requests have a timeout of at most 30 seconds. The initial submission
#' and final download are separate from the polling time allowance. Once ready,
#' the profile is downloaded and assigned the supplied name. Terrain-derived
#' profiles need field verification where vegetation or structures affect visibility.
#' @return A skyscapeR.horizon object for the requested panorama, as returned by \code{\link{downloadHWT}}. This function submits a remote request and waits for the calculated profile.
#' @param max.wait Maximum seconds to wait for panorama readiness after submission. Default 600.
#' @param lat The latitude of the location.
#' @param lon The longitude of the location.
#' @param elevation The elevation of the observer above
#' ground level in meters. Default is 1.6 meters (eye level).
#' @param name Name assigned to the returned horizon profile.
#' @param src Request source ID for \emph{HeyWhatsThat}. Default is
#'  'skyscapeR'. Only change this if you have been given a source ID by the
#'  creator of \emph{HeyWhatsThat}.
#' @param verbose Logical switch to control output. Default is \emph{TRUE}.
#' @export
#' @import utils
#' @references \href{https://www.heywhatsthat.com/}{HeyWhatsThat.com}
#' @seealso \code{\link{downloadHWT}}
#' @examples
#' \dontrun{
#' # Create and retrieve horizon data for the London Mithraeum:
#' hor <- createHWT(lat=ten(51,30,45), lon=ten("-0:5:26.1"), name='London Mithraeum')
#' }
createHWT <- function(lat, lon, name, elevation=1.6, src='skyscapeR', max.wait=600, verbose=TRUE) {
  if (length(max.wait)!=1L || !is.numeric(max.wait) || !is.finite(max.wait) || max.wait<=0)
    stop("max.wait must be a positive number of seconds.")
  if (!is.numeric(lat) || length(lat)!=1L || !is.finite(lat) || abs(lat)>90 ||
      !is.numeric(lon) || length(lon)!=1L || !is.finite(lon) || abs(lon)>180 ||
      !is.numeric(elevation) || length(elevation)!=1L || !is.finite(elevation))
    stop("Supply valid site coordinates and a finite observer height in metres.")
  if (verbose) .sky_message("Submitting panorama request to HeyWhatsThat.")
  reply <- .hwt_text("api/query",src=src,lat=lat,lon=lon,elev_agl=elevation)
  id <- strsplit(trimws(reply),"[[:space:]]+")[[1L]][1L]
  if (is.na(id) || !grepl("^[A-Za-z0-9]{8}$",id))
    stop("HeyWhatsThat did not return a valid panorama ID.")
  if (verbose) .sky_message("Waiting for HeyWhatsThat panorama ",id,".")
  deadline <- proc.time()[["elapsed"]] + max.wait
  Sys.sleep(min(120,max.wait))
  repeat {
    remaining <- deadline - proc.time()[["elapsed"]]
    if (remaining<=0.001) stop("Timed out waiting for panorama ",id,
                          ". Once ready, use downloadHWT('",id,"').")
    state <- tryCatch(trimws(.hwt_text("api/ready",src=src,id=id,timeout=min(30,remaining))),
                      error=function(e) stop("Could not check panorama ",id,": ",conditionMessage(e)))
    if (state=="1") {
      # A ready flag can precede publication of the HTML or CSV.
      hor <- tryCatch(downloadHWT(id),error=identity)
      if (!inherits(hor,"error")) break
      if (verbose) .sky_message("Panorama ",id," could not yet be retrieved: ",conditionMessage(hor))
    } else if (state!="0") stop("Unexpected HeyWhatsThat readiness response for ",id,": ",state)
    Sys.sleep(min(15,max(0,deadline-proc.time()[["elapsed"]])))
  }
  if (!missing(name)) hor$metadata$name <- name
  hor
}


#' Interpolate horizon altitude at a bearing
#' @description Estimate the altitude of the local horizon, and optionally its uncertainty, at one or
#' more true bearings.
#' @details A measured axis rarely falls exactly on a stored horizon sample. Separate cubic splines
#' interpolate altitude and its recorded uncertainty as functions of azimuth, and results
#' are rounded to two decimal places. Supply bearings within the profile coverage; a spline
#' beyond the measured range is an extrapolation. Distinguish this lookup from the linear
#' interpolation used inside coordtrans.
#' @param hor A \emph{skyscapeR.horizon} object from which to retrieve horizon altitude.
#' @param az Array of azimuth(s) for which to retrieve horizon altitude(s).
#' @param return.unc Logical controlling whether to return altitude uncertainty.
#' Default is \emph{FALSE}.
#' @return With \code{return.unc = FALSE}, a numeric vector of altitudes.
#' With \code{return.unc = TRUE} and one azimuth, an unnamed numeric vector
#' containing altitude followed by uncertainty. With multiple azimuths, a
#' numeric matrix with columns \code{alt} and \code{alt.unc}, one row per
#' input azimuth in the supplied order. Values are rounded to two decimal places.
#' @section Validation:
#' In H07, lookup at each edge of a north-crossing panorama returned its measured altitude and lookup beyond either edge returned NA. See the \href{https://stonexsky.github.io/skyscapeR/manual/#h07-panorama-coverage}{H07 check in the technical manual}.
#' @export
#' @import stats
#' @seealso \code{\link{createHor}}, \code{\link{downloadHWT}}
#' @examples
#' hor <- createHor(c(0,90,180,270,360), c(0,5,3,1,0),
#'                  alt.unc=c(0.1,0.2,0.3,0.2,0.1), loc=c(40,-8,100))
#' hor2alt(hor, 90)
#' hor2alt(hor, 90, return.unc=TRUE)
#' hor2alt(hor, c(90,180,270), return.unc=TRUE)
hor2alt = function(hor, az, return.unc=F) {
  if (identical(hor$metadata$source$provider,"Panorama")) {
    limits <- hor$metadata$source$coverage
    relative <- (az-limits[1])%%360
    angle <- (limits[1]+relative)%%360
    vals <- stats::approx(hor$data$az,hor$data$alt,xout=angle,rule=2)$y
    vals[relative>diff(limits)+1e-8] <- NA_real_
    if (!return.unc) return(round(vals,2))
    ans <- cbind(alt=round(vals,2),alt.unc=ifelse(is.na(vals),NA_real_,hor$data$alt.unc[1]))
    return(if (length(az)==1L) as.numeric(ans) else ans)
  }
  hh <- splinefun(hor$data$az, hor$data$alt)
  alt <- round(hh(az), 2)
  if (return.unc) {
    hh <- splinefun(hor$data$az, hor$data$alt.unc)
    unc <- round(hh(az), 2)
    if (length(az) == 1L) {
      return(c(alt, unc))
    }
    return(cbind(alt=alt, alt.unc=unc))
  }

  return(alt)
}

#' Find the northern declination limit of a horizon
#' @description Estimate the largest declination indicated by the horizon near true north.
#' @details Terrain relief can shift the limiting sightline away from the meridian. The function
#' samples bearings at 0.01-degree spacing within the specified half-width around 0
#' degrees. It calls az2dec at each bearing with the stored horizon, or with the supplied
#' latitude and constant altitude, then returns the maximum sampled declination. A limited
#' search sector may miss a more extreme sightline elsewhere.
#' @seealso \code{\link{createHor}}, \code{\link{downloadHWT}}, \code{\link{createHWT}}, \code{\link{exportHor}}, \code{\link{hor2alt}}, \code{\link{hor2min.dec}}
#' @return One numeric value in degrees, the largest declination sampled across the specified azimuth range around north.
#' @param hor Either a \emph{skyscapeR.horizon} object or the latitude of the place.
#' @param alt Altitude of the horizon at northernmost point. Only necessary if no
#' \emph{skyscapeR.horizon} object is provided. Default is 0.
#' @param range Range of azimuths to check for minimum declination. Default is 25
#' degrees either side of North.
#' @section Validation:
#' A05 recovered a limiting sightline within the northern search sector on a synthetic skyline with an interior crest. See the \href{https://stonexsky.github.io/skyscapeR/manual/#a05-solar-limits}{A05 check in the technical manual}.
#' @export
#' @examples
#' hor <- malta_horizons[[1]]
#' hor2max.dec(hor)
hor2max.dec = function(hor, alt=0, range=25) {
  xx <- seq(-range, range, 0.01)
  if (is(hor,'skyscapeR.horizon')) { decs <- az2dec(xx, hor) } else { decs <- az2dec(xx, hor[1], alt) }
  return(max(decs))
}

#' Find the southern declination limit of a horizon
#' @description Estimate the smallest declination indicated by the horizon near true south.
#' @details Terrain relief can shift the limiting sightline away from the meridian. The function
#' samples bearings at 0.01-degree spacing within the specified half-width around 180
#' degrees. It calls az2dec at each bearing with the stored horizon, or with the supplied
#' latitude and constant altitude, then returns the minimum sampled declination. A limited
#' search sector may miss a more extreme sightline elsewhere.
#' @seealso \code{\link{createHor}}, \code{\link{downloadHWT}}, \code{\link{createHWT}}, \code{\link{exportHor}}, \code{\link{hor2alt}}, \code{\link{hor2max.dec}}
#' @return One numeric value in degrees, the smallest declination sampled across the specified azimuth range around south.
#' @param hor Either a \emph{skyscapeR.horizon} object or the latitude of the place.
#' @param alt Altitude of the horizon at southernmost point. Only necessary if no
#' \emph{skyscapeR.horizon} object is provided. Default is 0.
#' @param range Range of azimuths to check for minimum declination. Default is 25
#' degrees either side of South
#' @section Validation:
#' A05 recovered a limiting sightline within the southern search sector on a synthetic skyline with an interior trough. See the \href{https://stonexsky.github.io/skyscapeR/manual/#a05-solar-limits}{A05 check in the technical manual}.
#' @export
#' @examples
#' hor <- malta_horizons[[1]]
#' hor2min.dec(hor)
hor2min.dec = function(hor, alt=0, range=25) {
  xx <- seq(-range, range, 0.01)+180
  if (is(hor,'skyscapeR.horizon')) { decs <- az2dec(xx, hor) } else { decs <- az2dec(xx, hor[1], alt) }
  return(min(decs))
}

#' @rdname exportHor
#' @export
exportStellarium <- exportHor
