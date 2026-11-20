# skyscapeR 2.0 (in development)

## Updating existing analyses

* Refit orientation tests after upgrading. The treatment of circular bearings,
  simulation grids, Monte Carlo ranks and local significance has changed.
  Saved fits retain their original probabilities. For results derived from
  azimuths, first rebuild the declination densities and summed probability
  distributions (SPDs), particularly for steep skylines and near-polar
  declinations.
* Recalculate skyline crossings, spatial-equinox targets and stellar first-
  and last-visibility dates used in earlier interpretations. Regenerate
  southern-hemisphere planispheres. Analyses based on a development version
  of the Milky Way outline should use the revised boundary.
* Replace `sky.sketch()` with `sketchSky()` followed by `plot(hor, sky = sky)`.
  Replace `sigTest()` with `randomTest()` or `modelTest()` according to the
  null hypothesis. For uncertain measurements use `spd()` in place of
  `curvigram()`. Remove `swapStars()` calls.
* Pass source azimuths and horizons explicitly when working with transformed
  densities. The [README migration guide](https://github.com/stonexsky/skyscapeR/blob/development-2.0/README.md#migrating-from-skyscaper-12-to-20)
  lists reordered arguments, retired functions and changes affecting earlier
  development builds.

## Field measurements and orientation tests

* Compass reduction uses IGRF-14 and pairs each site's survey date and
  location. Theodolite reduction preserves site rows and returns true
  bearings on the interval from 0 degrees inclusive to 360 degrees exclusive.
* `az.pdf()` and `dec.pdf()` support normal, uniform, Epanechnikov, triangular,
  tapered uniform and uniform openings convolved with normal survey error.
  For the last family, the default survey-error standard deviation is 0.5 degrees.
* `coordtrans()` retains full internal precision, refines integration on steep
  skylines and propagates altitude uncertainty through the spherical
  transformation near the celestial poles. `spd(normalise = TRUE)` gives unit
  area. `hpdi()` handles disjoint highest-density intervals.
* `randomTest()` tests random orientations. `modelTest()` accepts a specified
  null in azimuth or declination. An azimuth-defined null can also be tested
  in declination by transforming each observed and simulated measurement
  through its site horizon. Simulations preserve measurement shape.
* Monte Carlo ranks include ties and use the observed sample once in the
  reference ensemble. Local excesses and deficits of at least 0.5 degrees
  are assessed against the largest cluster found in each simulation.
  Circular tests count north once and can return intervals spanning north.
  Both tests support reproducible parallel simulation.

## Horizons, sky reconstructions and plots

* `createHor()` builds a skyline and shaded panorama from a supplied elevation
  model. `horFeature()` adds reconstructed features such as banks, stones and
  ditches. Lighting and surface colours can be adjusted independently of
  the skyline calculation.
* `importPanorama()` calibrates photographs and `importPeakFinder()` reads
  Stellarium ZIP landscapes from PeakFinder. Partial photographic skylines
  support analyses within their azimuth coverage. `exportPanorama()` writes a
  full-sphere PNG and `exportHor()` writes a Stellarium landscape.
* `sketchSky()` reconstructs the sky above a site. Natural views combine
  directional twilight, atmospheric extinction, stars, planets, lunar light,
  clusters and the Milky Way. `skyRadiance()` returns the brightness field
  before display compression. `skyStep()` moves through time.
* Sketches and planispheres support Western and Chinese constellation figures,
  the 28 lunar mansions and custom line tables. The Milky Way is traced on a
  0.1-degree grid. Its envelope and internal regions are illustrative. The
  seven-square-degree minimum dark-patch area was chosen visually. Publication
  sketches apply extinction to the band by default.
* Compass plots show epoch ranges as azimuth sectors. Stellar-phase plots
  compare stars or epochs in rings or strips. Significance plots distinguish
  local excesses and deficits and display the global probability.

## Astronomical calculations, data and documentation

* `star.phases()` uses Schaefer's twilight visibility model for heliacal rising
  and acronycal setting. The other passage events retain Purrington's rule.
  The event table records the method and date. Terrain clearance, circumpolar
  tracks and seasons spanning the December solstice are handled explicitly.
* `riseset()` resolves brief emergence and reoccultation against terrain.
  `spatial.equinox()` accepts an epoch and accounts for refraction and solar
  parallax. Calendar conversions handle day rollover and fixed time offsets.
  `EFM()` and `CC()` calculate lunar events against supplied horizons.
* The bundled stellar catalogue is unpacked in temporary storage. Ephemeris
  calculations require files covering the requested epoch, supplied through
  `swephRdata` or `SE_EPHE_PATH`. Structured astronomical results record
  the loaded skyscapeR and swephR versions.
* The vignette now follows field-survey and Maltese temple examples.
  Removed the RugglesCKR and RugglesRSC datasets. The
  [technical manual](https://stonexsky.github.io/skyscapeR/manual/) documents
  the methods, validation studies and limitations.

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
