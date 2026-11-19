# Rebuilding the illustrative Milky Way model

Milky Way illustrative image classification — revised 19 November 2026
Source: diffuse-light component of Deep Star Maps 2020, Ernie Wright,
NASA/Goddard SVS and ESA/Gaia/DPAC, https://svs.gsfc.nasa.gov/4851/.
The exact 1024 x 512 RGB source image used in visual selection is retained in
 data-raw/milky-way/source-map.png (SHA256 in metrics.json).
It is an image-tone model, not calibrated surface-brightness photometry.

Reproduction: python data-raw/milky-way/rebuild.py, then
Rscript data-raw/milky-way/save-grid.R, from the repository root.
Python dependencies: numpy, scipy, Pillow. R uses base functions only.

RGB tone weights 0.2126, 0.7152, 0.0722; values on a 0–255 scale.
Bilinear resampling into Galactic coordinates; 1-degree Gaussian presmoothing.
Outer envelope: further 4-degree Gaussian smoothing, tone >=40; largest
connected region near the Galactic plane, with enclosed holes filled.
Dark regions: presmoothed tone more than 8% below an 8-degree Gaussian local
reference. Connected dark regions must occupy at least 7 square degrees.
Bright and brightest: tone >72.47510740164273 and >115.30130722988615,
respectively; connected regions of at least 10 square degrees, excluding dark
regions. These fixed values preserve the reviewed brightness classes: originally
1.10 and 1.75 times a reference mean of 65.88646127422066 in the pilot model.
They are not recomputed after a change to the dark mask.

All Gaussian widths are FWHM. Smoothing approximately corrects longitude width
by cos(latitude); components use 8-neighbour connectivity with periodic longitude.
Areas use spherical cell solid angles expressed in square degrees.
Thresholds were selected by visual comparison. The 7-square-degree cutoff is
illustrative, not a physiological or independently validated visibility limit.
Dark patches remain within the outer Milky Way envelope for geometric alignment.

Grid: 0.1-degree sampling; longitude -180 to <180; latitude -60 to +60.
R raw matrix [longitude, latitude]: 3600 x 1201.
Class codes: 0 outside, 1 ordinary, 2 dark, 3 bright, 4 brightest.
Fine sampling does not imply 0.1-degree independent detail in the source image.
Rendering weights, surface-brightness normalisation and atmospheric transmission
are unchanged. Dependent geometric results must be recalculated for this outline.
The CSV gives extremal occupied-cell edges per longitude for legacy line plotting.
