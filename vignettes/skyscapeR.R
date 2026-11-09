## ----setup, include=FALSE-----------------------------------------------------
spreadsheet_examples <- requireNamespace("readxl", quietly = TRUE)
# Accept file-backed positions from either swephRdata or SE_EPHE_PATH.
covers_years <- function(years) {
  tryCatch({
    skyscapeR::star("Aldebaran", year = years)
    TRUE
  }, skyscapeR_ephemeris_error = function(e) FALSE)
}
astronomy_examples <- covers_years(c(-3399L, -3099L, 2016L, 2026L))
epoch_examples <- astronomy_examples && covers_years(c(-3999L, -1999L))
knitr::opts_chunk$set(error = FALSE, collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4, dpi = 300, fig.retina = 1,
                      out.width = "100%", fig.align = "center")

## ----own-file, eval=FALSE-----------------------------------------------------
# library(skyscapeR)
# survey <- readxl::read_excel("/Users/yourname/Research/field-survey.xlsx", sheet = "Compass")

## ----load-example, include=FALSE----------------------------------------------
library(skyscapeR)
if (spreadsheet_examples) survey <- readxl::read_excel(system.file("extdata", "field-survey.xlsx", package="skyscapeR"), sheet="Compass")

## ----show-survey, eval=spreadsheet_examples-----------------------------------
survey

## ----reduce-survey, eval=spreadsheet_examples, message=FALSE------------------
locations <- as.matrix(survey[c("latitude", "longitude", "elevation_m")])
compass <- reduct.compass(loc = locations, mag.az = survey$mag_az_deg,
                          date = as.Date(survey$survey_date), ID = survey$id)
compass[c("ID", "Magnetic.Azimuth", "Mag.Dec", "True.Azimuth")]

## ----solar-sight, eval=astronomy_examples, message=FALSE----------------------
library(skyscapeR)
site_location <- c(51.18, -1.83, 105)
theo <- reduct.theodolite(loc = site_location, az = c(298.40, 302.35),
                          date = "2016/02/20", time = "11:07:17",
                          tz = "Europe/London", az.sun = ten(327, 29, 50),
                          limb = "right")
theo$Altitude <- c(2, 5)  # Horizon altitudes measured along the two sightlines
theo[c("Uncorrected.Azimuth", "Sun.Az", "True.Azimuth", "Altitude")]

## ----compass-decs, eval=spreadsheet_examples----------------------------------
survey$true_az_deg <- compass$True.Azimuth %% 360
survey$dec_deg <- az2dec(survey$true_az_deg, loc = survey$latitude,
                         alt = survey$horizon_alt_deg)
survey[c("id", "true_az_deg", "horizon_alt_deg", "dec_deg")]

## ----theo-decs, eval=astronomy_examples---------------------------------------
theo$Declination <- az2dec(theo$True.Azimuth, loc = theo$Latitude, alt = theo$Altitude)
theo[c("True.Azimuth", "Altitude", "Declination")]

## ----own-targets, eval=spreadsheet_examples && astronomy_examples-------------
candidates <- findTargets(survey$dec_deg[1] + c(-2, 2), c("3400 BC", "3200 BC"))
candidates$stellar

## ----entrance-----------------------------------------------------------------
door <- malta_entrances[1, ]
door[, c("Name", "Latitude", "Longitude", "Az L", "Az R")]
edges <- c(door[["Az L"]], door[["Az R"]])
door_frame <- c(edges[1], mean(edges), edges[2])

## ----surveyed, eval=FALSE-----------------------------------------------------
# profile <- read.csv("/Users/yourname/Research/surveyed-horizon.csv")
# survey_location <- c(51.18, -1.83, 105)  # latitude, longitude, elevation in metres
# surveyed <- createHor(az = profile$azimuth, alt = profile$altitude,
#                       loc = survey_location)
# plot(surveyed)

## ----photo, eval=FALSE--------------------------------------------------------
# photo <- importPanorama("doorway-east.png", loc = survey_location,
#                         az.lim = c(75, 115), alt.lim = c(-5, 10), calibrate = TRUE)
# plot(photo, show.skyline = TRUE)

## ----hwt, eval=FALSE----------------------------------------------------------
# hwt <- createHWT(lat = door$Latitude, lon = door$Longitude,
#                  elevation = 1.6, name = "Mnajdra South")

## ----peakfinder, eval=FALSE---------------------------------------------------
# pf <- importPeakFinder("landscape.zip", svg = "panorama.svg")

## ----local-dem, eval=FALSE----------------------------------------------------
# terrain <- createHor(dem = "survey-dem.tif", loc = c(door$Latitude, door$Longitude),
#                      eye.height = 1.6, max.distance = 10000, image.dir = tempdir())
# plot(terrain)

## ----door-horizon-------------------------------------------------------------
hor <- malta_horizons[[1]]
hor2alt(hor, edges)
plot(hor, xrange = c(75, 115))
abline(v = door_frame, col = "firebrick", lty = c(1, 2, 1))
door_dec <- az2dec(edges, loc = hor)

## ----door-targets, eval=astronomy_examples------------------------------------
period <- c("3400 BC", "3100 BC")
targets <- findTargets(range(door_dec), period)
targets$solar
targets$stellar

## ----equinox, eval=astronomy_examples-----------------------------------------
year <- "3250 BC"
sun <- riseset(hor, date = "3250 BC/04/17", obj = "sun",
               calendar = "J", timezone = "UTC")$rise
sun
sky <- sketchSky(hor, time = paste(sun$date, sun$time), calendar = "J", timezone = "UTC")
equinox <- skyObject("Equinox", year = year, loc = hor, col = "steelblue", lty = 3)
sunrise_track <- skyObject(sun$declination, year = year, loc = hor, col = "orange")
plot(hor, sky = sky, obj = list(equinox, sunrise_track), xrange = c(85, 105))
abline(v = door_frame, col = "firebrick", lty = c(1, 2, 1))

## ----phase-circle, eval=astronomy_examples, fig.height=7----------------------
phases <- star.phases("M 45", year, hor)
plot(phases)
events <- phases$metadata$events
hr_date <- events$date.julian[events$event == "heliacal rising"]
hr_date

## ----phase-linear, eval=astronomy_examples------------------------------------
plot(phases, layout = "linear")

## ----cluster-rise, eval=astronomy_examples------------------------------------
predicted_jd <- time2jd(hr_date, calendar = "J", timezone = "UTC")
rise <- riseset(hor, jd = predicted_jd + 5, obj = "M 45", calendar = "J", timezone = "UTC")$rise
rise
rise_jd <- time2jd(paste(rise$date, rise$time), timezone = "UTC", calendar = "J")
sky <- sketchSky(hor, jd = rise_jd + 5 / 1440)
plot(hor, sky = sky, xrange = c(85, 105))
abline(v = door_frame, col = "firebrick", lty = c(1, 2, 1))
target <- as.numeric(sky$clusters[sky$clusters$identifier == "M 45", c("az", "alt")])
segments(target[1] + 0.18, target[2] + 0.14,
         target[1] + 0.70, target[2] + 0.58, col = "grey85", lwd = 0.7)
text(target[1] + 0.75, target[2] + 0.60, "Pleiades", adj = c(0, 0), col = "grey85", cex = 0.78)

## ----alcyone-dawn, eval=astronomy_examples------------------------------------
alcyone <- star.phases("Alcyone", year, hor)
events <- alcyone$metadata$events
date <- events$date.julian[events$event == "heliacal rising"]
alcyone_day <- time2jd(date, calendar = "J", timezone = "UTC")
alcyone_rise <- riseset(hor, jd = alcyone_day + 5, obj = "Alcyone", calendar = "J", timezone = "UTC")$rise
jd <- time2jd(paste(alcyone_rise$date, alcyone_rise$time), timezone = "UTC", calendar = "J")
sky <- sketchSky(hor, jd = jd + 16 / 1440)
plot(hor, sky = sky, xrange = c(85, 105), yrange = c(0, 11))
abline(v = door_frame, col = "firebrick", lty = c(1, 2, 1))
target <- as.numeric(sky$stars["eta Tau", c("az", "alt")])
segments(target[1] + 0.18, target[2] + 0.14,
         target[1] + 0.70, target[2] + 0.58, col = "grey85", lwd = 0.7)
text(target[1] + 0.75, target[2] + 0.60, "Alcyone", adj = c(0, 0), col = "grey85", cex = 0.78)

## ----later-rise, eval=astronomy_examples--------------------------------------
later <- riseset(hor, jd = predicted_jd + 60, obj = "M 45", calendar = "J", timezone = "UTC")$rise
sky <- sketchSky(hor, time = paste(later$date, later$time), calendar = "J", timezone = "UTC")
sector <- c(85, 105)
plot(hor, sky = sky, xrange = sector, yrange = c(0, 11))
abline(v = door_frame, col = "firebrick", lty = c(1, 2, 1))

## ----later-step, eval=astronomy_examples--------------------------------------
sky <- skyStep(hor, sky, minutes = 10, xrange = sector, yrange = c(0, 11))
abline(v = door_frame, col = "firebrick", lty = c(1, 2, 1))

## ----group-bearings, fig.width=5.75, fig.height=5.75, out.width="82%"---------
rows <- which(malta_entrances$Phase == "Ggantija")
sites <- malta_entrances[rows, ]
profiles <- malta_horizons[rows]
width <- (sites[["Az R"]] - sites[["Az L"]]) %% 360
centre <- (sites[["Az L"]] + width / 2) %% 360
plotAzimuth(centre)

## ----compass-targets, eval=astronomy_examples, fig.width=5.75, fig.height=5.75, out.width="82%"----
solstices <- skyObject("solar extremes", year = year, loc = hor)
plotAzimuth(centre, obj = solstices, marks = "rim", obj.position = "outside")

## ----group-opening-widths-----------------------------------------------------
plotBars(centre, unc = width / 2, names = sites$Name, unit = "Azimuth")

## ----group-decs---------------------------------------------------------------
bearings <- az.pdf(pdf = "uniform", az = centre, unc = width / 2, name = sites$Name)
decs <- coordtrans(bearings, profiles)
plot(spd(decs))

## ----linked-doorway-----------------------------------------------------------
plot(decs, index = 1, az = bearings, hor = profiles, hdr = 0.68)

## ----random-null, eval=FALSE--------------------------------------------------
# set.seed(42)
# fit <- randomTest(decs, az = bearings, hor = profiles, nsims = 500)

## ----random-result------------------------------------------------------------
fit <- readRDS(system.file("extdata", "ggantija-random-500.rds", package = "skyscapeR"))
plot(fit, show.local = TRUE)
fit$metadata$local.pval

## ----peak-targets, eval=astronomy_examples------------------------------------
bands <- fit$metadata$local.pval
peak <- subset(bands, type == "+" & p.value < 0.05)
candidates <- findTargets(c(peak$start, peak$end), period)
candidates$stellar

## ----gacrux, eval=astronomy_examples------------------------------------------
gacrux <- skyObject("Gacrux", year = period)
plot(fit, obj = gacrux, show.local = TRUE)

## ----gacrux-rise, eval=astronomy_examples-------------------------------------
i <- which(sites$Name == "Skorba West")
site_hor <- profiles[[i]]
opening <- c(sites[["Az L"]][i], sites[["Az R"]][i])
g_phases <- star.phases("Gacrux", year, site_hor)
g_events <- g_phases$metadata$events
g_date <- g_events$date.julian[g_events$event == "heliacal rising"]
g_day <- time2jd(g_date, calendar = "J", timezone = "UTC")
g_rise <- riseset(site_hor, jd = g_day + 5, obj = "Gacrux",
                  calendar = "J", timezone = "UTC")$rise
g_rise
g_jd <- time2jd(paste(g_rise$date, g_rise$time), timezone = "UTC", calendar = "J")
sky <- sketchSky(site_hor, jd = g_jd + 20 / 1440)
plot(site_hor, sky = sky, xrange = c(122, 144), yrange = c(-1, 10))
abline(v = c(opening[1], mean(opening), opening[2]), col = "firebrick", lty = c(1, 2, 1))
target <- as.numeric(sky$stars["gam Cru", c("az", "alt")])
segments(target[1] + 0.18, target[2] + 0.14,
         target[1] + 0.70, target[2] + 0.58, col = "grey85", lwd = 0.7)
text(target[1] + 0.75, target[2] + 0.60, "Gacrux", adj = c(0, 0), col = "grey85", cex = 0.78)

## ----sky-styles, eval=astronomy_examples, fig.show="hold", fig.cap=c("Natural", "Publication colour", "Publication monochrome")----
moon_sector <- c(119, 157)
moon_heights <- c(72.7, 77.4)
when <- "2026/12/21 20:00"

natural <- sketchSky(hor, time = when, timezone = "UTC", style = "natural")
plot(hor, sky = natural, xrange = moon_sector, yrange = moon_heights)
cluster <- subset(natural$clusters, identifier == "M 45")
moon <- natural$bodies$moon
segments(cluster$az - 4, cluster$alt - 0.62, cluster$az - 1.7, cluster$alt - 0.28,
         col = "grey70", lwd = 0.7)
text(cluster$az - 4, cluster$alt - 0.62, "Pleiades", pos = 2, col = "grey85", cex = 0.8)
segments(moon[1] + 1.35, moon[2] + 0.22, moon[1] + 3.4, moon[2] + 0.58,
         col = "grey70", lwd = 0.7)
text(moon[1] + 3.4, moon[2] + 0.58, "Moon", pos = 4, col = "grey85", cex = 0.8)

colour <- sketchSky(hor, time = when, timezone = "UTC", style = "publication-colour")
plot(hor, sky = colour, xrange = moon_sector, yrange = moon_heights)

mono <- sketchSky(hor, time = when, timezone = "UTC", style = "publication-mono")
plot(hor, sky = mono, xrange = moon_sector, yrange = moon_heights)

## ----compare-epochs, eval=epoch_examples, fig.height=7------------------------
ald <- lapply(c("4000 BC", "3250 BC", "2000 BC"),
              function(y) star.phases("Aldebaran", y, hor))
plot(ald[[1]], ald[-1])

## ----western-sky, eval=astronomy_examples, fig.cap="Modern Western constellations above the eastern horizon at Mnajdra South."----
sector <- c(0, 180)
heights <- c(-5, 85)
sky <- sketchSky(hor, time = when, timezone = "UTC", style = "publication-colour",
                 display = list(constellations = "modern"))
plot(hor, sky = sky, xrange = sector, yrange = heights)

## ----chinese-constellations, eval=astronomy_examples, fig.cap="Chinese lunar-mansion asterisms visible in the same eastern sky."----
sky <- sketchSky(hor, time = when, timezone = "UTC", style = "publication-colour",
                 display = list(constellations = "chinese_mansions"))
plot(hor, sky = sky, xrange = sector, yrange = heights)

## ----custom-constellations, eval=FALSE----------------------------------------
# my_figures <- importConstellations("my-constellationship.fab")
# # Subset rows here if only some figures should be drawn.
# sky <- sketchSky(hor, time = "2026/12/21 20:00", timezone = "UTC",
#                  display = list(constellations = my_figures))
# plot(hor, sky = sky, xrange = sector)
# abline(v = door_frame, col = "firebrick", lty = c(1, 2, 1))

## ----milky-way-natural, eval=astronomy_examples, fig.cap="The summer Milky Way above the southern horizon at Mnajdra South, 15 July 2026 at 23:00 UTC."----
mw_time <- "2026/07/15 23:00"
mw <- sketchSky(hor, time = mw_time, timezone = "UTC",
                display = list(moon = FALSE))
plot(hor, sky = mw, xrange = c(115, 235), yrange = c(-2, 65))

## ----milky-way-guides, eval=astronomy_examples, fig.cap="The same field in publication colour. The Galactic equator is dot-dashed, the ecliptic dashed and the celestial equator dotted."----
guides <- sketchSky(hor, time = mw_time, timezone = "UTC",
                    style = "publication-colour", display = list(moon = FALSE,
                    galactic.equator = TRUE, equator = TRUE, ecliptic = TRUE))
plot(hor, sky = guides, xrange = c(115, 235), yrange = c(-2, 65))

## ----milky-way-off, eval=astronomy_examples, fig.cap="The same natural view with the Milky Way switched off; the stars and atmosphere are still drawn."----
without_band <- sketchSky(hor, time = mw_time, timezone = "UTC",
                          display = list(moon = FALSE, milky.way = FALSE))
plot(hor, sky = without_band, xrange = c(115, 235), yrange = c(-2, 65))

## ----efm, eval=astronomy_examples---------------------------------------------
full_moon <- EFM(year = 2026, loc = c(40, 0, 400), season = "spring", timezone = "UTC")
full_moon$Moon
full_moon$Sun

## ----cc, eval=astronomy_examples----------------------------------------------
crescent <- CC(year = 2026, loc = c(40, 0, 400), season = "winter", crescent = "first",
               timezone = "UTC", simulation = "SIM1")
crescent$Moon
crescent$diagnostics[, c("status", "visibility.margin", "day.offset")]

## ----planisphere, eval=astronomy_examples, fig.height=7-----------------------
planisphere(year, lat = door$Latitude, max.mag = 3, show.galactic.equator = TRUE)

## ----planisphere-mansions, eval=astronomy_examples, fig.height=7, fig.cap="Chinese lunar-mansion asterisms on the planisphere for Mnajdra South's latitude and the chosen epoch."----
planisphere(year, lat = door$Latitude, max.mag = 3,
            show.galactic.equator = TRUE, constellations = "chinese_mansions")

