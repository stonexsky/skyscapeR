# skyscapeR technical manual

Open `docs/index.html` for the rendered book. To rebuild, install Quarto and R, change to this directory and run `quarto render`. Rendering uses base R's tools package to convert both complete help pages and calculate the heliacal-rising summary from archived CSVs. The Rd files come from the 2.0.9004 source archive. The archives contain the heliacal-rising study and its 2.0.9003 source, and the H01 international field validation for 2.0.9004.

To publish this draft, commit the source and the rendered `docs/` directory to the intended GitHub repository, and select that branch's `/docs` folder as the GitHub Pages source. A `.nojekyll` marker is included. No service credentials are needed to read the manual locally.

The draft includes a complete function index, nine annotated chapter outlines, worked chapters on horizons, stellar phases, sky reconstruction and significance testing, a linked specialist methods catalogue and versioned individual help pages, the H01 and heliacal-rising studies, a register of sky-rendering studies, and an introduction to the astronomical dependency. H01 includes the international field comparison, its downloadable report and reproducible archive. The reproducibility archives preserve the data used for each study. The later methods chapters include executed examples for the documented development build. The ephemeris page records upstream accuracy claims and the data generation of the inspected files. It contains no new ephemeris benchmark.

The study register links to original review archives and a downloadable extract. Before public deployment, replace the access-controlled review links with permanent public archival links.

The density and event-search chapters include executed examples with synthetic horizons, two figures and an archive of R scripts and outputs.

The lunar crossover and candidate-screening chapters include executed examples, comparison with the 2012 simulation series and two additional evidence downloads.

This draft adds worked chapters for astronomical conventions, field reduction and panorama calibration, with executable examples and their outputs.

The vignette guide covers the archived package tutorial section by section. With the matching package installed, run `Rscript scripts/check-coverage.R` from the project root to compare the catalogue with the namespace, formal argument lists and explicit calls in every vignette R chunk. The check covers documentation structure. Source review is still needed when an implementation changes behaviour without changing a function signature.
