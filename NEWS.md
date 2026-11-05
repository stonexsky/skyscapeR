# skyscapeR 2.0 (release candidate)

* Natural sky sketches use a shared twilight adaptation level for stars, retaining local moonlight and directional sky contrast. Natural views attenuate the Milky Way with altitude. Set `display = list(constellations = "Orion")` to select figures. Dubhe and Zeta Herculis are drawn once each using their combined system entries.


* Natural sky sketches draw unresolved cluster light over its catalogue footprint, with a smooth transition to visible stars.
* Sky sketches and planispheres draw the Pleiades from their member stars. Integrated cluster entries, galaxies, transients and reference points are omitted from the background star layer. Target searches and phase calculations retain the full catalogue.

* Follow a stellar phase event to a dated horizon sketch using the Julian date from `star.phases()` and the stellar crossing time from `riseset()`. Compare phase results with `plot(first, list(second, third))`. Significance plots place the global p-value in the upper right.

* Constellation lines in sketchSky remain visible when endpoint stars are faint or hidden; star-symbol visibility is unchanged.

- Stellarium sky captures use a single equirectangular image. Supply a calendar time, choose display layers and move forward or backward with `skyStep()`.
- Pass a colour, a Stellarium capture or a reusable vector background from `sketchSky()` to horizon plots via `sky`. Image-backed plots keep a 2:1 frame. Use `xrange` and `yrange` for package plot ranges.
- Export any site horizon as a transparent-sky, full-sphere PNG with `exportPanorama()`.
- Existing scripts with positional arguments should consult the migration table in the README.
- Import and calibrate photographic horizons with `importPanorama()` and capture a Stellarium sky for horizon plots with `stellariumSky()`.
- DEM panoramas use observer-assisted relief lighting by default. Set `lighting = "directional"` for the previous shading. This changes the image only.
- For transformed declinations, `modelTest()` now shifts each observed density to its sampled mean by default. Explicit `unc` retains the normal-density simulation.
- Interactive runs of `randomTest()` and `modelTest()` show progress through simulations, envelope calculation and departure scoring, including multi-core runs.
- Routine console updates use R messages; `options(skyscapeR.verbose = FALSE)` silences them for the session.

- North-crossing azimuth PDFs use contiguous bearings in coordinate transformation.
  HeyWhatsThat altitude and uncertainty now interpolate across the final bin
  and north. Regenerate affected declination PDFs, SPDs and simulations.
- SPD plots accept celestial target groups. PDF labels fit the plot, combined
  panels put the version above the figure, and azimuth horizon axes show
  their major ticks. Plots keep their data coordinates active for annotations and advance through user-defined layouts.
- Star catalogue lookups accept identifiers such as M45 and M 45, returning
  the common name with the coordinates. The optional ephemeris installation
  notice follows the welcome banner.


- swephRdata is optional for installation and strongly recommended. Astronomical
  positions require verified Swiss files, including stellar calculations.

- The package manual and README now index the workflows by research question.
  The reference, vignette and messages use the term horizon throughout.

- Reorganised the README and migration guide, refreshed the package reference
  and consolidated related regression scripts by feature.

- `EFM()` selects the first sufficiently illuminated Moon with the required
  rising or setting bearing on the same local date as the solar event. Explicit
  timezone and calendar settings now apply throughout. Topocentric bearings
  govern selection in both declination output frames. Retry failed historical
  calls with UTC or the intended fixed offset. Rerun previous EFM series with
  their original site, threshold, refraction and time conventions.

- star.phases() now reuses the riseset() terrain solver for all crossings on a
  fixed stellar track. Every model day is classified using positive-duration
  visibility windows. Separated seasons and cyclic transitions are retained.
  The yearly calculation reuses the track and evaluates solar samples on demand.

- Added demPalettes() and six named terrain.palette presets with separate
  land and water ramps and per-ramp overrides.
- DEM panoramas accept the standalone land and sea palettes, elevation threshold
  and number of sea colours. Lighting direction, ambient fill and distance haze
  are configurable. The image metadata records the render settings.

- Added a DEM mode to createHor() with georeferenced ray sampling, optional
  terrestrial refraction and an external shaded panorama shared by plotting
  and Stellarium export. Terrain images carry the supplied DEM credit.

- Added maximum-cluster adjustment for local intervals in both significance tests,
  using the same selection, 0.5-degree filter and both-sign search in simulations.
  Global pooled scores are unchanged. New local inference tag distinguishes old fits.
- Replaced time-grid-only horizon bracketing with terrain-corner visits, azimuth
  turning/rapid-motion refinement and interior-extremum searches. Detects brief
  emergence/reoccultation; pure tangencies are not crossings. Profiles must close
  at north. Rerun affected horizon events and dependent crossover calculations.

- Corrected global Monte Carlo calibration in `randomTest()` and `modelTest()`:
  pool the observation once with simulations for envelope, mean and scale;
  retain inclusive ranks, scientific nulls, batched disk storage and core defaults.
  Simulation-constant departures now participate in pooled scoring. Local
  intervals receive the maximum-cluster correction described above. Rerun previous fits; outputs carry an inference tag.
  Document the change from the earlier published-method implementation and
  the exchangeability rationale (Dufour 2006).

- Corrected `CC()` visibility to use geocentric arc of light, as explicitly
  defined by Caldwell and Laney, while retaining topocentric lunar altitude.
  Rerun earlier development-version CC results. Independent scientific
  validation also identified small-simulation significance calibration (corrected above) and
  brief horizon-crossing limitations (addressed by the new search above).


- Added `CC()` for the four crescent crossover events, following Silva and
  Pimenta (2012) with Caldwell–Laney visibility at sunset/sunrise. Includes
  SIM1/SIM2, geocentric declinations, numerical horizon support, diagnostic
  output and explicit missing results and selection diagnostics.

- Added `importPeakFinder()` for Stellarium ZIP horizons and optional matching
  Vector SVG ridge drawings. Validated tangent-projection calibration creates
  one external PNG; warning/fallback preserves horizon-only use on failure.
- Horizon plots and `sky.sketch()` use attached PNGs on demand. `exportHor()`
  includes PNGs in spherical landscapes; `exportStellarium()` is its new alias.
- Added optional `xml2` (SVG import) and `png` (image plotting) dependencies,
  provenance, checksum validation, explicit uncertainty/elevation caveats and
  regression tests. Images never enter coordinate or simulation calculations.

* Reorganised the vignette around Maltese entrance and terrain data. A shared field record leads to separate single-site visibility and group-pattern investigations. The specified-null case study is deferred.

* Expanded the existing four vignette sections to cover survey reduction, horizon data, both significance tests, astronomical targets, events and sky plots. The offline examples run during rendering. Fixed graphics-state restoration in plotAzimuth.

* Pair survey dates and sites explicitly for IGRF-14. Both HeyWhatsThat functions use optional httr2, HTTPS, checked responses and bounded waits. Downloaded site elevations are in metres. See the README for rerun guidance.

* Added method and data references to the help pages. The stellar visibility threshold is attributed to Purrington (1988, equation 1), with its magnitude and horizon limits documented. Schaefer (1986, 1987) and Brady (2015) are cited for extinction, visibility modelling and stellar phases. These documentation changes require no rerun.

* Annual stellar sampling recovers short search windows. Solar date searches refine turning points and solve directly against the ephemeris. Named seasons now report dates for the requested epochs and calendar. Updated the migration guide and removed the unused rootSolve dependency.

* Revised all reference topics around the scientific question and the implemented method. Added a workflow guide to the README and removed generated vignette copies and local session files from the source tree.

* Rewrote the significance-test help to distinguish the random-azimuth null (Silva 2020) from the specified-model extension (Silva and Lomsdalen 2025).
* Fixed historical daily sampling and date labels, planisphere hour labels, and installation with long ephemeris search paths. Compressed the lazy-data database with xz and updated three regression tests.

* Global and local p-values in randomTest() and modelTest() use the upper-tail Monte Carlo rank (b + 1)/(B + 1), including ties. Both help pages compare the convention with Silva (2020), equation 2.

* randomTest() and modelTest() compare densities directly with the reported envelope, then standardise the departures. Exact envelope ties contribute zero.

* randomTest() and modelTest() flag departures from constant simulated density in data$constant.departures and exclude those grid points from significance scores.

* Both significance tests report local intervals at least 0.5 degrees wide. Constant null envelopes retain their simulated density and exact matches contribute zero departure.

* randomTest() defaults to one core and supports reproducible parallel runs in azimuth and declination. It shares batched grid storage with modelTest(). Removed snow, doSNOW and foreach.

* modelTest() supports reproducible parallel simulations with ncores. Large runs process the simulation grid in batches using temporary storage.

* Replaced MESS with local grid integration and grDevices colour handling.

* The significance and stellar-phase plots now include their default colours directly.

* Removed three mandatory external dependencies (doParallel, pointr, stringr) and the unused suggested package testthat.

* B02a: modelTest() tests user-supplied weighted models in azimuth or declination. dec.pdf() constructs declination PDFs from measurements and their uncertainties.

* Revised the package reference, bundled data descriptions and introductory vignette.
* Removed the bundled RugglesCKR and RugglesRSC datasets and their help pages because redistribution permission has not been obtained. Vignette examples now use explicitly synthetic declinations.

* plotAzimuth() aligns its eight compass spokes and labels with bearings measured clockwise from north.

* B19d: coordtrans() uses absolute transformed uncertainty widths and clearly
  rejects unsupported zero widths and invalid densities; see help for limits.

* B19c: coordtrans() activates .prec as azimuth integration spacing and adds
  .output.res for returned grid spacing. Defaults preserve numerical grids;
  new transformed PDFs supply all three spacings to randomTest simulations.

* B19b: randomTest() now honours normalise for azimuth input plus a horizon,
  scaling empirical and simulated curves after transformation to declination.

* B19a: spd(normalise=TRUE) now gives unit trapezoidal area, rather than
  unit sum of heights. Normalised scales change; see spd help.

* B18c: star.phases() respects refraction and accepts horizon obstruction
  from horizon objects. New solar.alt separates twilight from terrain; see
  help for defaults, retained approximations and minute-sampling limits.

* B18b: riseset() uses horizon profiles and corrects the inverted refraction
  flag. Profile events are first crossings within 48 hours; missing events
  return NA. Explicit alt conflicts with a profile. See help for scan limits.

* B18a follow-up: sky.sketch() retains qualifying stars when no catalogue
  magnitudes are missing; missing/faint-only selections remain safely empty.

* B18a: sky.sketch() derives stellar years from observation JD for both stars
  and constellation lines, fixing numeric-JD and truncated BCE input years.

* B17: findTargets() preserves tables when no rows need removing, handles
  empty stellar selections, and passes scalar lower/upper dating bounds.

* B16: exportHor() isolates temporary files, restores the working directory,
  fixes the default name and checks the ZIP. Existing archives are refused;
  the saved path is returned invisibly.

* B14: significance-test results with an empty local p-value table print
  "No local regions detected." instead of failing on nonexistent rows.
  Global results and nonempty local-row formatting are preserved.

* B13: antizenith() correctly reads latitude from horizon metadata, restoring
  documented horizon-object support. Formula and range check are unchanged;
  an offline horizon example and regression checks cover both input types.

* B12: eq(verbose=FALSE) correctly returns geocentric zero declination when
  no location is supplied. Verbosity no longer changes the calculation path;
  topocentric corrections are unchanged. Help and quiet-call example updated.

* B11: bernoulli.trial() uses stable binomial probability functions, handles
  zero hits correctly, and avoids factorial overflow and upper-tail
  cancellation. Invalid inputs now error; finite scalar inputs and the
  independent-trial/common-probability assumptions are documented.

* B10: hor2alt(return.unc=TRUE) returns an alt/alt.unc matrix for multiple
  azimuths instead of overwriting an altitude and losing uncertainties.
  Single-azimuth pairs and altitude-only output are unchanged. Return shapes
  are documented with offline examples.

* B09 follow-up: coordtrans() uses linear interpolation to prevent negative
  density undershoot, evaluates the upper bound to avoid extrapolation, and
  documents calculation versus output spacing. PDF normalisation and SPD
  scaling are preserved. Numerical results can change; regenerate saved
  spline-based transformations. Coarse grids can still miss narrow peaks.

* B09: hpdi() handles PDF objects as lists of ordered interval matrices;
  repairs multiple-region bounds and validates mass, densities and grids.
  Negative densities (including coordtrans spline undershoot) now error.
  Integer-grid gaps remain supported; see hpdi help for limitations.

* B08: reduct.compass() accepts magdec without a date (Date is NA), requires
  at least one of those inputs, reports per-measurement horizon altitudes
  consistently with declination calculations, and handles single locations.

* B07: repaired year-range celestial overlays in plotBars() and horizon plots.
  Use stored minimum/maximum declinations, full-height target bands, preserved
  line styles and consistent refraction for both bounding tracks.

* B06c: randomTest() requires explicit az and hor for declination input,
  validates sources/settings and no longer retrieves objects by stored names.
  Compact simulation specifications are exported to workers. Fixed parallel
  verbose=FALSE options and cluster cleanup (related B15). B06 source-input
  changes are now implemented; full-package validation remains pending.

* B06b: PDF plots now use explicit named az/hor inputs for combined panels;
  omit both to plot distributions alone. show.az is deprecated and ignored
  with a warning. Global-name lookup removed; legacy objects follow the same
  explicit-input rules. Declination tests take the same explicit inputs.

* B06a: coordtrans() now returns compact version-2 metadata
  with descriptive source labels and no copies of the source objects.

* B05: fixed body.position() horizon metadata lookup in both coordinate frames
  and preserved horizontal output tables with refraction disabled. Clarified
  location inputs and output in help; added regression tests and an example.

* B04: repaired time2jd()/jd2time() date rollover and actual-year timezone
  handling; added explicit fixed offsets, date validation, DST ambiguity errors
  and nearest-second output rounding. Named zones now require Gregorian years
  1900–9999. For ancient dates, use timezone="UTC", a fixed offset, or set the
  package timezone to UTC. Updated conversion and timezone-setting help.

* B03: ten() now preserves the sign of numeric angles, using the first non-zero
  component. Later negative components are rejected. Character input retains
  a leading minus, including negative zero degrees. Added negative-angle help
  examples and regression tests. This changes results for negative numeric
  inputs that previously returned positive magnitudes.

* B01: use one bundled stellar catalogue at inst/ephemeris/sefstars.txt,
  configured on namespace load without overwriting swephR or swephRdata files.
  All R-side readers use the same catalogue. Preserve SE_EPHE_PATH and clear
  any previously cached stellar catalogue on initialisation.
* Retired swapStars() and its help page; removed the bundled original catalogue
  and the stars selection field from skyscapeR.vars(). Remove swapStars() calls
  from existing scripts: no replacement call is needed.

# skyscapeR 1.2 development history

* added dates optional facility to findTargets()
* fixed issues with sun and moon reporting on findtargets()
* Added modelTest()
* Added all stars up to magV 6, i.e. all visible stars
* Added planetary.extremes()
* Fixed issue with reporting of sunrise/set dates in findTargets()
* Fixed issue with az2dec()
* Fixed manual page for randomTest()
* Fixed plotting issue with plotAzimuth()
* Fixed typos in vignette
* Added parameter col to plot.skyscapeR.sigtest()
* coordtrans() can now transform skyscapeR.spd objects
* coordtrans() has been sped up through vectorisation
* Changed randomTest() to use SNOW, now allowing for progress bar
* randomTest() can now perform the test on multiple measurements at the same site/horizon much faster
* Small fixes to riseset()
* Corrected issues with exportHor()
* Changed star() to use internal version of stellar ephemeris and naming (replacing swephR's version)
* Complete rewrite of the engine behind star.phases() to allow for seasonality calculations outside of swephR's time interval
* Fixed issue in star.phases() to do with alt.hor and circumpolar stars
* Added facility to plot multiple star.phases objects
* Added hor2min.dec() and hor2max.dec() to calculate min/max decs for given horizon profile or location
* Added sky colour option to plot.skyscapeR.horizon()
* Fixed issue with jd2time where day was not correct for change of date
* Fixed issue with riseset where date was not displayed accurately when the day changed

* Added the full Stellarium 0.22.2 Chinese constellation set and a separate
  28 lunar-mansion asterism set, with provenance and coverage notes.
* Constellation imports retain zero-segment entries; single-star figures draw
  as rings. Overlapping named sets are deduplicated when selected together.

* Completed Chinese endpoint coverage using SIMBAD astrometry and Hipparcos
  V magnitudes: 59 records appended and one existing HIP alias linked, with
  per-field provenance and explicit unknown-radial-velocity placeholders.
