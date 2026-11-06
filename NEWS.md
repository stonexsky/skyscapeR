# skyscapeR 2.0 (development)

This development version is being tested ahead of the 2.0 release. For detailed
migration instructions and positional-argument changes, see the README.

## Vignette and help

* Other tools now have separate subsections for sky styles, constellation
  patterns, the Milky Way, twilight and Stellarium captures. The planisphere
  example also shows the Chinese lunar mansions.
* `planisphere()` accepts bundled constellation sets, selected figures and
  imported line tables through `constellations`, with `col.constellations`
  controlling their colour. Lines follow spherical arcs and are clipped to the
  accessible sky; faint endpoints are retained independently of `max.mag`.
* Help clarifies methods, defaults, return values and limitations of the revised
  astronomical and plotting functions. The vignette uses named `CC()` arguments.
* EFM moonset pairing and some CC event selections remain under investigation.
  Their help pages now describe these known limitations; this update does not
  change either calculation to incorporate the experimental corrections.

## Changes requiring attention

* Structured astronomical results now carry a `provenance` attribute recording
  loaded skyscapeR and swephR versions. Access it with `attr(result, "provenance")`.
  List entries, table columns and numeric-vector return types are unchanged.
  Preserve original results, since downstream transformations may drop attributes.

* `star.phases()` results now retain the remaining visibility settings, the
  numerical horizon profile and loaded package versions. Existing fields and
  phase calculations are unchanged; older saved objects lack the added record.

* `planisphere()` now uses a south-polar projection at negative latitudes,
  with consistent star selection, RA labels, ecliptic and Milky Way layers.
  The exact latitude-dependent outer radius is retained. The ecliptic now uses
  geocentric coordinates without local refraction. Regenerate affected charts.
* `sketchSky()` help explicitly distinguishes cluster and twilight rendering
  from calibrated naked-eye detection. The approved rendering is unchanged.

* `spatial.equinox()` now accepts an epoch and consistent atmospheric settings.
  It selects first emergence and final disappearance from terrain crossings,
  inverts refraction and solar parallax consistently, and reports undefined
  crossings as `NA` with a warning. Existing two-argument calls remain valid.
  `findTargets()` uses the requested period endpoints for this candidate.
  Recalculate spatial-equinox candidates and downstream interpretations.

* Single `star.phases()` results now plot as a circular seasonal cycle by default.
  Comparisons retain linear strips; use `layout = "linear"` for the previous
  single-result view or `layout = "circular"` for concentric comparisons.
  Event labels respect seasons crossing December solstice, and the event table
  correctly names entry into and exit from the curtailed-passage season.

* Replace `sky.sketch()` with `sketchSky()` followed by a horizon `plot()`.
  Choose constellation figures through `display$constellations`.
* Replace `sigTest()` with `randomTest()` or `modelTest()`, according to the
  null hypothesis. Replace `curvigram()` with `spd()` for uncertain measurements
  or `stats::density()` for a kernel density estimate of point measurements.
* Remove `swapStars()` calls. A single bundled stellar catalogue is configured
  automatically without overwriting files belonging to `swephR` or `swephRdata`.
* Development helpers `phaseAlignment()` and `plotStarPhases()` have been
  removed. Use the event date from `star.phases()`, then `riseset()` and
  `sketchSky()`. Compare phase objects with `plot(first, list(second, third))`.
* Supply `az` and `hor` explicitly when testing transformed declination PDFs or
  plotting their transformation. Objects are no longer recovered by stored names.
  The deprecated PDF plotting argument `show.az` is ignored with a warning.
* Use `xrange` and `yrange` for package plots. Horizon plots retain compatibility
  with `xlim` and `ylim`. Several functions have reordered arguments, so review
  positional calls against the README migration table.
* Rerun significance tests, affected horizon events and downstream analyses
  after upgrading. The statistical calibration, coordinate transformation,
  calendar handling and crossing calculations have changed as described below.

## Field measurements and probability distributions

* Magnetic corrections use IGRF-14 with explicitly paired survey dates and site
  locations. `reduct.compass()` accepts a supplied magnetic correction without
  a date and evaluates the horizon at each corrected bearing.
* `ten()` preserves the sign of numeric angles and character negative-zero
  degrees. `az2dec()` pairs measurements with their supplied horizons.
* `hor2alt(return.unc = TRUE)` returns altitude and uncertainty columns for
  multiple bearings. HeyWhatsThat interpolation handles north and the final bin.
* `coordtrans()` uses linear interpolation, avoids negative spline undershoot,
  handles north-crossing azimuth distributions, and validates transformed widths
  and densities. `.prec` controls integration spacing and `.output.res` controls
  the returned grid. Compact metadata records transformation settings.
* `dec.pdf()` constructs declination densities from measurements and uncertainty.
  `spd(normalise = TRUE)` now gives unit trapezoidal area. `hpdi()` validates
  densities and correctly handles multiple highest-density intervals.

## Statistical tests

* `randomTest()` tests the random-orientation null. `modelTest()` supports
  weighted, user-supplied models in azimuth or declination. For transformed
  densities, its default simulation shifts each observed density to the sampled
  mean. Supplying `unc` retains normal-density simulation.
* Both tests pool the observation once with the simulations when calculating
  envelopes, means and scales. Upper-tail Monte Carlo ranks include ties and use
  `(b + 1) / (B + 1)`. Exact envelope ties contribute zero departure, and
  simulation-constant departures participate in pooled scoring.
* Local intervals must span at least 0.5 degrees and receive a maximum-cluster
  correction using the same interval selection and both-sign search in each
  simulation. Outputs carry an inference tag. These changes affect inference,
  including results from earlier development versions, and require fresh fits.
* Parallel simulations are reproducible and use batched temporary storage.
  The default worker count uses half the available logical cores, rounded down
  with a minimum of one, subject to simulation and memory limits. Interactive
  runs report progress through simulation, envelope and scoring stages.
* `randomTest()` honours normalisation after azimuth-to-declination conversion.
  `bernoulli.trial()` uses stable binomial probabilities, including zero hits.
* Significance plots place the global p-value in the upper-right corner and
  support celestial overlays. Empty local-result tables print cleanly.

## Astronomical positions, dates and events

* `star.phases()` now determines horizon clearance from the stellar track and
  supplied skyline with the selected refraction setting in both hemispheres.
  `metadata$horizon.status` separates always-above, always-below and crossing
  tracks from seasonal visibility. Completely obstructed tracks return an
  all-invisible result; non-setting tracks are no longer rejected. Recalculate
  affected southern, circumpolar and raised-horizon cases.

* `time2jd()` and `jd2time()` handle date rollover, actual-year timezone rules,
  fixed offsets and whole-second rounding. Ambiguous daylight-saving times are
  rejected. Named timezones require Gregorian years 1900–9999. Use UTC or a
  known fixed offset for ancient dates.
* Historical solar searches refine turning points and solve against the
  ephemeris. Annual stellar sampling recovers short event-search windows, and
  returned dates use the requested epoch and calendar.
* `riseset()` supports stellar and solar-system skyline crossings. Terrain-corner
  sampling, rapid-motion refinement and interior-extremum searches detect brief
  emergence and reoccultation. Pure tangencies are not crossings. The refraction
  flag is corrected, and missing events return `NA`.
* `star.phases()` reuses the terrain crossing solver on a fixed stellar track,
  classifies positive-duration visibility windows and retains separated seasons
  and cyclic transitions. `solar.alt` separates twilight from terrain height.
  Event metadata provide Julian-calendar dates for subsequent rise/set searches.
* `EFM()` selects the first sufficiently illuminated Moon with the required
  rising or setting bearing on the solar event's local date. Calendar and
  timezone settings apply throughout, and selection uses topocentric bearings.
* `CC()` calculates the four crescent crossover events following Silva and
  Pimenta (2012), with Caldwell–Laney visibility at sunset or sunrise. Visibility
  uses geocentric arc of light and topocentric lunar altitude. Horizon profiles,
  diagnostic output and explicit missing results are supported.
* Fixed horizon metadata handling in `body.position()` and `antizenith()`, quiet
  geocentric calculations in `eq()`, and equinox overlays in `skyObject()`.
  `findTargets()` handles empty results and scalar dating bounds consistently.
* Catalogue lookups recognise identifiers such as `M45` and `M 45`. Background
  star plots omit integrated clusters and other non-stellar entries, while target
  searches and phase calculations retain them.

## Horizons and landscape images

* `createHor()` accepts a user-supplied DEM, with georeferenced ray sampling,
  optional terrestrial refraction and an external shaded panorama. Image
  metadata preserve render settings and the supplied DEM credit.
* `demPalettes()` supplies six terrain presets. Land and sea ramps, lighting,
  ambient fill and distance haze are configurable. Observer-assisted relief
  lighting is the default. Image shading does not alter horizon calculations.
* `importPanorama()` imports photographic horizons with calibration facilities.
  `importPeakFinder()` imports Stellarium ZIP horizons and optional SVG ridge
  drawings, with calibration validation and a horizon-only fallback.
* `createHWT()` and `downloadHWT()` use HTTPS, checked responses and bounded
  waits through optional `httr2`. Downloaded elevations are converted to metres.
* Horizon plots load attached PNGs on demand. `exportHor()` includes images in
  spherical landscapes, isolates temporary files, checks the archive and refuses
  to overwrite existing exports. `exportStellarium()` is an alias.
  `exportPanorama()` exports a transparent-sky, full-sphere PNG.

## Sky visualisation and plotting

* `sketchSky()` returns a reusable sky background with natural, publication
  colour and monochrome styles. `stellariumSky()` captures an equirectangular
  Stellarium view, and `skyStep()` advances or reverses the selected time.
* Natural skies use shared twilight adaptation for stars with local moonlight
  and directional contrast. The Milky Way is attenuated differentially with
  altitude. These displays remain approximations of visual perception.
* Unresolved clusters can appear as diffuse light over their catalogue footprint,
  fading as individual stars become visible. This is an illustrative rendering,
  not a calibrated prediction of extended-object detection. The Pleiades' member
  stars are plotted individually, without an integrated M45 point symbol.
* Constellation segments are drawn even when their endpoint stars are too faint
  to display. Figures can be selected from modern Western, full Chinese and
  Chinese lunar-mansion sets. `importConstellations()` accepts user-supplied
  figures in Stellarium's legacy line format. Single-star figures use rings and
  overlapping selections are deduplicated.
* Single and multiple stellar-phase plots share a layout and show only applicable
  seasons in the legend. Labels include the star name, Bayer designation and
  epoch. Single-object plots retain event annotations.
* Improved `plotBars()` labels, compass spokes in `plotAzimuth()`, planisphere
  hour labels, celestial range overlays and PDF panel layouts. Plot coordinates
  remain available for annotations and user-defined layouts.

## Data, installation and documentation

* Compressed the stellar catalogue and lazy-data database with lossless xz
  compression. The stellar catalogue is expanded into session-temporary storage.
* Added Stellarium's Chinese figures and the 28 lunar mansions as separate sets.
  A documented SIMBAD/Hipparcos supplement adds 59 stellar records and links one
  existing HIP identifier to complete constellation endpoint coverage.
* Removed `RugglesCKR` and `RugglesRSC` because redistribution permission has not
  been obtained. The vignette uses field examples and Maltese temple data.
* Consolidated stellar catalogue provenance and retained separate source-specific
  notices for bundled datasets, constellation lines, cluster sizes and colours.
* `swephRdata` is optional at installation. Astronomical calculations require
  verified Swiss Ephemeris files covering the requested epoch, including stellar
  calculations. Existing `SE_EPHE_PATH` settings are preserved.
* Reduced dependencies and removed obsolete parallel backends. Regression
  scripts are not bundled in the package. Routine progress uses R messages and
  can be silenced with `options(skyscapeR.verbose = FALSE)`.
* Rebuilt the vignette around data input, single-site analysis, multiple-site
  analysis and other tools. The README includes installation and migration
  guidance and links to the illustrated online vignette. Help pages describe
  methods, assumptions, references and workflows by research question.

# skyscapeR 1.2 development history

The following summarises work recorded during the earlier development cycle.
Later 2.0 changes supersede implementation details where noted above.

* Added optional date reporting to `findTargets()` and corrected solar and lunar
  target reporting, sunrise/set dates and `az2dec()` calculations.
* Added `modelTest()` and `planet.extremes()`, and expanded the stellar catalogue
  to include stars down to visual magnitude 6.
* Extended `coordtrans()` to SPD objects and improved speed through vectorisation.
* Added SNOW-based parallel support and progress reporting to `randomTest()`,
  with faster handling of repeated measurements at a shared site. This backend
  has since been replaced in 2.0.
* Reworked `star.phases()` for wider epoch coverage and corrected horizon and
  circumpolar handling. Added comparison plots for multiple phase objects.
* Added `hor2min.dec()` and `hor2max.dec()`, horizon sky colours and significance
  plot colours. Corrected `plotAzimuth()` and horizon export behaviour.
* Updated stellar catalogue lookup, daily date labels, rise/set reporting,
  documentation and vignette examples.
