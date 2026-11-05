# skyscapeR

**Development release ahead of version 2.0.** We are still testing the package and fixing bugs before the public release. Please [report any problems](https://github.com/stonexsky/skyscapeR/issues) with the code that reproduces them.

skyscapeR is an R package for analysing the orientations of archaeological structures. It takes measured bearings and their uncertainties, combines them with site locations and horizon profiles, and expresses the results as declination distributions. You can compare the observed orientations with random orientations or with a model of your own.

Upgrading from an earlier release? [Check which analyses and figures to rerun](#migrating-from-skyscaper-12-to-20).

## Install the development version

In R 4.1 or later, install the 2.0 development branch:

```r
install.packages("remotes")
remotes::install_github("stonexsky/skyscapeR", ref = "development-2.0")
```

For stellar positions, solar sights and other ephemeris calculations, install the Swiss Ephemeris data once and restart R. Alternatively, set `SE_EPHE_PATH` to your existing ephemeris directory before loading the package. The files must cover your epoch.

```r
install.packages("swephRdata", repos = "https://rstub.r-universe.dev", type = "source")
```

## Start with your observations

The [worked vignette](vignettes/skyscapeR.Rmd) starts with an Excel field survey. It follows compass and theodolite readings through reduction, declination and candidate selection, then uses Maltese temple entrances for single-site reconstruction and a multi-site significance test.

For a short example using three true bearings, a level horizon and an illustrative site at 35° N, 8° W:

```r
library(skyscapeR)
az <- az.pdf(az = c(85, 95, 112), unc = 2)
hor <- createHor(az = c(0, 360), alt = c(0, 0), loc = c(35, -8, 100))
dec <- coordtrans(az, hor)
plot(spd(dec))
```

For the illustrated vignette in RStudio, install `knitr`, `rmarkdown` and `readxl`, then repeat the GitHub installation with `build_vignettes = TRUE, force = TRUE`. Open it with `vignette("skyscapeR", package = "skyscapeR")`. Use `?skyscapeR` for a guide to the help pages.

## Explore a site or compare sites

| Question | Where to start |
|---|---|
| What could rise or set in an opening? | Use `az2dec()` and `findTargets()` to screen candidates. For stars, use `star.phases()` for the season and `riseset()` for the skyline crossing on the selected date. |
| What would the view look like? | Use `sketchSky()` with the site's horizon and time, then `plot(hor, sky = sky)`. Choose natural, publication colour or monochrome views, and use `skyStep()` to advance time. |
| How do constellation traditions compare? | In `sketchSky()`, pass `display = list(constellations = "modern")`, `"chinese"` or `"chinese_mansions"` as the constellation selection. For your own line figures, use `importConstellations()` with a Stellarium legacy line file. See `?constellations` for sources and identifiers. |
| Do multiple sites share an orientation pattern? | Start with `plotAzimuth()` and `plotBars()`. Combine declination distributions with `spd()`, then use `randomTest()` or compare a specified hypothesis with `modelTest()`. |

For changes by version, see [NEWS](NEWS.md). For bundled dataset credits and terms, see [COPYRIGHTS](inst/COPYRIGHTS).

## Migrating from skyscapeR 1.2 to 2.0

Start at the **earliest affected step** below and rerun everything that used its output. For example, if a west longitude changed sign, repeat the field reduction before rebuilding the declination curves. If only a compass grid was drawn incorrectly, redraw the figure from the existing bearings. Keep the old results alongside the rerun and record the original inputs, package and ephemeris versions, calendar, timezone, refraction settings, uncertainty model, grid spacings and simulation seed.

### Removed functions and deprecated interfaces

Version 2.0 removes obsolete entry points rather than retaining defunct wrappers. Update scripts before loading the new version:

| Removed interface | Migration action |
|---|---|
| `sky.sketch()` | Construct a georeferenced horizon, then call `sky <- sketchSky(hor, time = "2026/12/21 20:00", timezone = "UTC", style = "natural", display = list(sun = TRUE, moon = TRUE))` and `plot(hor, sky = sky, xrange = c(30, 150))`. Unlike the old function, `sketchSky()` returns a reusable sky background; plotting is a separate step. The old `loc` vector must be used to construct a horizon first. Move the old `const` selection into `display$constellations`. Translate `max.mag` and visibility switches explicitly; there is no one-to-one positional argument replacement. |
| `sigTest()` | Use `randomTest()` to test against random orientations, or `modelTest()` to compare against a specified model. These are not identical statistical tests; re-run the analysis and review the null hypothesis and settings. The `skyscapeR.sigTest` *result class* and its print/plot methods remain supported. |
| `curvigram()` | Use `spd()` for summed probability densities of uncertain measurements, or `stats::density()` for a kernel density estimate of point measurements. These are different estimators; choose according to the data and research question. |

The deprecated `show.az` argument of `plot.skyscapeR.pdf()` is ignored. Replace it with `plot(dec, index = i, az = az, hor = hor)` for a combined transformation view, or omit both `az` and `hor` for a density-only plot. This compatibility argument is retained temporarily with a warning; do not use it in new scripts.

### Positional arguments

This table compares the public v1.2.0 package with 2.0. To migrate a positional call, use the v1.2 column to identify what each value meant, then pass those values by name. Arguments after `...` must be named. For the surveyed-profile constructor, the last column gives a named call so you can reuse the original survey inputs directly.

| Function | v1.2 arguments in order | 2.0 arguments in order or named survey call |
|---|---|---|
| `az.pdf()` | `pdf, az, unc, name, verbose, .cutoff, .res` | `az, unc, pdf, name, .res, .cutoff, verbose` |
| `body.position()` | `obj, time, timezone, calendar, dec, loc, refraction, atm, temp, verbose` | `time, obj, loc, timezone, calendar, dec, refraction, atm, temp, verbose` |
| `coordtrans()` | `pdf, hor, xrange, refraction, atm, temp, verbose, .res, .prec` | `pdf, hor, xrange, .res, .prec, .output.res, verbose` |
| `createHor()` | `az, alt, alt.unc, loc, name, smooth, .scale` | `createHor(az = bearings, alt = heights, alt.unc = uncertainty, loc = site, name = label, smooth = FALSE, .scale = 1000)` |
| `createHWT()` | `lat, lon, elevation, name, src, verbose` | `lat, lon, name, elevation, src, max.wait, verbose` |
| `EFM()` | `season, rise, year, loc, min.phase, refraction, atm, temp, timezone, calendar` | `year, loc, season, rise, min.phase, timezone, calendar, refraction, atm, temp` |
| `exportHor()` | `hor, name, author, description, ground_col, hor_col` | `hor, name, author, description, use.image, peakfinder.layout, max.texture.width, ground_col, hor_col` |
| `findTargets()` | `decrange, timerange, max.mag, loc, dates, calendar` | `decrange, year, loc, max.mag, dates, calendar, timerange` |
| `modelTest()` | `pdf, model, nsims, conf, tails, normalise, ncores, save.sim, verbose, .res` | `pdf, model, unc, nsims, conf, tails, normalise, .res, ncores, verbose` |
| `orbit()` | `dec, loc, res, refraction, atm, temp` | `dec, loc, refraction, atm, temp, res` |
| `planet.extremes()` | `obj, epoch, loc` | `obj, year, loc, epoch` |
| `planisphere()` | `epoch, lat, max.mag, show.ecliptic, col.ecliptic, show.equator, col.equator, star.col` | `year, lat, max.mag, show.ecliptic, col.ecliptic, show.equator, col.equator, star.col, epoch, show.milky.way, show.galactic.equator` |
| `plot.skyscapeR.horizon()` | `x, show.az, xlim, ylim, obj, refraction, col.ground, col.sky, ...` | `x, show.az, xrange, yrange, obj, col.ground, col.sky, use.image, ..., show.skyline, sky, xlim, ylim, fastLoad` |
| `plot.skyscapeR.pdf()` | `x, index, hdr, show.az, xlim, col, ...` | `x, index, hdr, show.az, xrange, col, ..., az, hor` |
| `plot.skyscapeR.sigTest()` | `x, xlim, title, col, show.pval, show.local, pal, ...` | `x, xrange, title, col, show.pval, show.local, pal, ..., obj` |
| `plot.skyscapeR.spd()` | `x, xlim, ylim, title, col, shading, ...` | `x, xrange, yrange, title, col, shading, ..., obj` |
| `plotAzimuth()` | `az, col, lwd, lty, obj, show.obj.labels` | `az, obj, show.obj.labels, col, lty, lwd` |
| `plotBars()` | `val, unc, names, unit, col, shade, mark, sort, xrange, yrange, obj, show.obj.label` | `val, unc, names, unit, obj, show.obj.label, xrange, yrange, sort, shade, mark, col` |
| `randomTest()` | `pdf, .res, nsims, conf, tails, normalise, ncores, verbose, hor, refraction, atm, temp` | `pdf, az, hor, nsims, conf, tails, normalise, .res, ncores, verbose` |
| `riseset()` | `obj, date, jd, alt, loc, calendar, timezone, dec, refraction, atm, temp, verbose` | `loc, date, jd, obj, alt, timezone, calendar, dec, refraction, atm, temp, verbose` |
| `sky.objects()` | `names, epoch, loc, col, lty, lwd` | `names, year, loc, col, lty, lwd, epoch` |
| `sky.sketch()` (removed) | Direct plot with `time, loc, ...` | Use `sketchSky(hor, time = ..., display = ..., style = ...)` followed by `plot(hor, sky = sky)`; see the removed-functions table above. |
| `spd()` | `pdf, normalise, xrange, .cutoff, .res` | `pdf, xrange, normalise, .res, .cutoff` |
| `star.phases()` | `star, year, loc, alt.hor, k, limit, alt.rs, res, refraction, atm, temp` | `star, year, loc, alt.hor, alt.rs, solar.alt, k, limit, refraction, atm, temp, res` |

Use `xrange` and `yrange` for package plots. Existing `xlim` and `ylim` calls to horizon plots still work. Use `sky = "darkblue"` for a solid colour or pass a sky returned by `sketchSky()` or `stellariumSky()`.

For example, replace `randomTest(pdf, 0.5)` with `randomTest(pdf, .res = 0.5)`. Remove explicit atmosphere arguments to `randomTest()`, `coordtrans()` and horizon plots. Set session conditions with `skyscapeR.vars()` for a fresh transformation, or regenerate the observed declination PDF with the intended conditions. `modelTest()` no longer accepts the developer diagnostic `save.sim` or adds `data$simulated` to new fits.

### Field measurements, dates and astronomical targets

| If your analysis used… | What to check or rerun |
|---|---|
| Automatic magnetic correction or downloaded horizons | Rerun `mag.dec(loc, dates)` and `reduct.compass()` from the survey dates and original readings, especially for multi-date or multi-site calls. Use one row per site and align dates with bearings. IGRF-14 uses the survey date at midnight UTC, including leap years. For saved HeyWhatsThat profiles, download the original panorama ID again if site metadata affect downstream results: reported feet are now converted to metres and hemisphere signs are parsed independently. Install `httr2` to create or retrieve online panoramas. |
| Signed angles from `ten()` | Compare with the field record, especially west longitudes, south latitudes and negative altitudes. `ten(-8, 30)` now returns −8.5; `ten("-0:30")` returns −0.5. An existing `-ten(8, 30)` is still valid. Use one sign for the whole angle. Recalculate site-dependent results if a coordinate changed. |
| Local dates, Julian dates or historical calendar labels | Reconvert original timestamps with `time2jd()` and check them with `jd2time()` across midnight, year boundaries and daylight-saving changes. Specify `calendar = "G"` or `"J"` and the intended timezone. Before 1900, set `skyscapeR.vars(timezone = "UTC")` for UT input or supply the known fixed offset for local time. Resolve ambiguous daylight-saving times with the recorded offset. Compare rounded whole seconds as well as the date; astronomical year 0 is 1 BCE. Recompute positions and events derived from changed instants. `planet.extremes()` and `solar.date()` now use UT for year-based sampling; rerun historical tables, including short-year and BCE date labels. |
| Compass bearings or horizon lookups | Rerun `reduct.compass(loc, mag.az, magdec = correction)` when you know the correction, or pass the survey `date` to estimate it. Check altitude and declination at **each corrected bearing**, particularly on a varying horizon. For multiple bearings, `hor2alt(hor, bearings, return.unc = TRUE)` returns a matrix with `alt` and `alt.unc` columns; one bearing returns a pair. In `az2dec()`, match each bearing to its own horizon when supplying a list. |
| Body positions, target lists or epoch ranges | Rerun `body.position()` calls that used horizon objects or disabled refraction. Match the original site, geocentric/topocentric frame and intended geometric/apparent altitude. Retry previously failing `eq(verbose = FALSE)` and `antizenith(hor)` calls. Recheck the resulting solar, lunar and stellar candidate tables. |
| Stellar search windows and seasonal dates | Rerun `findTargets(decrange, timerange, dates = TRUE)` with the original inputs, including empty results and searches shorter than a century. Recreate stellar range overlays with `sky.objects()`. Stellar dates now use annual positions and rounded boundary years. Rerun near-solstice `solar.date(dec, year, calendar)` calls, including formerly rejected values or numeric text in date cells. Keep both crossings even when their day labels coincide. Code reading all dates should use `result[-1, , drop = FALSE]` to include a possible third passage. Rebuild named equinox and solstice tables for the original year(s) and calendar. A slash separates differing endpoint-year dates. |
| Fixed stars or a previously switched catalogue | Remove `swapStars()` and the old `stars` setting. Reload the package in a fresh session and rerun stellar identifications, positions, target lists and visibility calculations if they relied on a different catalogue. If an older installation overwrote a dependency's star file, reinstall that dependency to restore its distributed copy. Version 2.0 selects its bundled catalogue automatically and uses covering Swiss Ephemeris files from `swephRdata` or an explicitly configured `SE_EPHE_PATH`. Apply manual topocentric or sidereal engine settings after loading skyscapeR. |
| Equinoctial full moons | Rerun `EFM(year = year, loc = loc, season = season, rise = rise, timezone = "UTC", calendar = "G")` with your original site, illumination threshold and refraction settings. Substitute the intended timezone and calendar. For historical dates, passing UTC or a fixed offset on the call is sufficient. The chronological search recovers three missing moonsets at the Silva–Pimenta reference site. At other sites, thresholds or timezones, dates can also change. Bearings now always use apparent topocentric positions, including when geocentric declination is requested. |
| Rising, setting or stellar visibility | Rerun `riseset()` even for a level horizon: the old constant-altitude route inverted the refraction flag. If an older call used the horizon only for its location, extract `site <- as.numeric(hor$metadata$georef)` and use `loc = site, alt = h` for the original fixed obstruction, or `loc = hor` alone for a horizon. Supply a profile that closes at 0°/360°. The search now includes terrain corners and brief crossings between minute samples. A profile gives the first upward and downward crossings separately within 48 hours; a missing crossing returns `NA`. Rerun `star.phases()` with the original site, star and atmospheric settings. All 365 model days are now evaluated, including brief terrain openings. Season lengths, phase codes and transition dates can change even on a level horizon. Check for `visibility.search = "terrain-windows-v1"` in the metadata. Use coordinates with `alt.hor` for a constant obstruction, or a horizon object alone for the actual horizon. The twilight threshold defaults to −6° with a profile and to `alt.hor - 6` with coordinates. Set `solar.alt` explicitly when comparing them. |

### Rebuilding distributions and significance tests

**If you retained the original bearings and horizons**, rebuild the azimuth PDF if the bearings or their uncertainties changed. Otherwise, reuse that PDF. Always regenerate an older spline-based declination curve before reanalysing it. A saved declination PDF's metadata labels cannot reconstruct the bearings or horizon. Recover the source objects or recreate them from the original records; preserve the measurement order when using one horizon per site.

The following rerun uses the `az` and `hor` from the first example. Use the actual survey's refraction and atmosphere settings, measurement uncertainties and simulation count when repeating a research analysis.

```r
invisible(skyscapeR.vars(refraction = FALSE))
dec <- coordtrans(az, hor, .res = 0.1,
                  .prec = 0.1, .output.res = 0.01, verbose = FALSE)
summed <- spd(dec, normalise = TRUE)
regions <- hpdi(dec)                # one interval matrix per measurement
plot(dec, index = 1, az = az, hor = hor)
set.seed(42)
fit <- randomTest(dec, az = az, hor = hor, nsims = 10, .res = 2,
                  normalise = TRUE, ncores = 1, verbose = FALSE)
```

- **Recheck peaks, areas and interval boundaries.** `coordtrans()` now interpolates linearly. Its `.res` controls declination calculations, `.prec` now controls azimuth integration, and `.output.res` controls the stored grid. Refine the calculation and integration grids when checking narrow peaks. Check horizon coverage and measured altitude uncertainties if a rerun rejects a zero or nonfinite transformed width. A genuinely zero width is unsupported by this Gaussian approximation. `spd(normalise = TRUE)` now produces unit area on the returned grid. Earlier releases scaled the sampled heights to sum to one, so recompute normalised curves and anything based on their heights. `hpdi(dec)[[i]]` returns ordered intervals for measurement `i`; `hpdi(summed)` returns one matrix. Preserve disjoint regions, including those split at north.
- **Choose the original random-orientation hypothesis and route.** Use `randomTest(az, ...)` for azimuth alone. Use `randomTest(az, hor = hor, ...)` for the shared-horizon shortcut, which transforms the summed curve and now honours `normalise`. To test individually transformed observations, use the call in the example. Its `.res` sets the comparison grid. To change the stored transformation settings, rebuild `dec` first. These routes can produce different curves even at one site.
- **For a specified null model**, supply a two-column table of candidate coordinates and nonnegative weights, or one such table per observation in order. Use `modelTest(az, az.model, ...)` for an azimuth model. For transformed declinations use `modelTest(dec, dec.model, ...)` to shift each observed curve to its sampled mean. Supply `unc` explicitly to reproduce simulations using normal declination errors. A random-azimuth null transformed through a horizon generally differs from a uniform-declination model. Set `.res` explicitly when comparing grids with an earlier analysis.
- **Rerun both significance tests from inputs before reporting p-values.** Saved fits retain the old numbers. Recheck the envelope, global p-value and local intervals. Crossings are now decided in density units; exact ties contribute zero. Local runs must span at least 0.5° inclusive. Examine `fit$data$constant.departures` for differences at zero simulated variance; these now enter the pooled score when the observation differs. The new rank includes ties and the observation itself: the minimum p-value is `1 / (nsims + 1)`. Local values now use maximum-cluster adjustment across qualifying intervals and both signs under the complete null. Cluster endpoints are descriptive, not confidence limits; rerun older local results. New fits identify the methods as `metadata$inference = "pooled-exchangeable-v1"` and `metadata$local.inference = "maximum-cluster-v1"`. The default is 1,000 simulations using half the available logical cores (at least one). Changing the null, uncertainty, normalisation or grid can change the scores as well as the ranks. Keep the seed, `nsims`, `conf` and `tails` for reproducibility. The same seed reproduces draws across 2.0 worker counts; older parallel runs used a different sampling arrangement.
- **For discrete hit counts**, rerun `bernoulli.trial(n, p, r, type = "tail")` for at least `r` hits, or use `type = "single"` for exactly `r`. Check previously undefined or unstable results, especially zero hits, large samples and small tails.

For a HeyWhatsThat horizon, rebuild declination PDFs for north-crossing bearings with the original azimuth PDFs and horizon, then rerun dependent SPDs and significance tests. Saved horizons can be reused.

### Saved figures, datasets and exports

| Output | Migration action |
|---|---|
| Compass, horizon and target-overlay figures | Redraw `plotAzimuth()` to align the cardinal spokes and labels. Recreate year-range targets with `sky.objects()` and redraw `plotBars()` and horizon plots to show both declination bounds or tracks. Regenerate `planisphere()` to update the right-ascension hour labels. Recreate former `sky.sketch()` figures with `sketchSky()` to apply the observation epoch to stars and constellations, especially with numeric Julian dates or BCE dates, and for catalogue magnitude selections that formerly dropped qualifying stars. |
| PDF and test plots or printed tables | Use `plot(dec)` for densities alone; for transformation panels use `plot(dec, index = i, az = az, hor = hor)`. Replace the old PDF plotting switch `show.az` with these named inputs. Print or plot a newly computed significance fit for updated statistics. Empty local tables now print correctly and azimuth intervals are labelled in azimuth. |
| Removed datasets | `RugglesCKR` and `RugglesRSC` are no longer distributed. Scripts using them must import a separately obtained, appropriately authorised copy, or use their own observations. |
| Stellarium landscapes | To rebuild a landscape use `path <- exportHor(hor, name = "My-site")`. The return is the archive path. Choose a fresh name or move the previous ZIP before rerunning. Reimport the archive to update the landscape. |

The switch to built-in plotting colours and the dependency cleanup do not, by themselves, require rerunning scientific analyses. The numerical changes above and the input choices in the help pages are the relevant comparisons.
