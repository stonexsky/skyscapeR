# skyscapeR technical manual

Open `docs/index.html` for the prebuilt Quarto book. The methods guide groups chapters by research question. The function index links to individual methodological entries and versioned R help. The validation register links to study reports, results and reproduction archives.

## Editing and rebuilding

Edit the relevant chapter at the project root. Function pages use the `fn-` prefix. Validation pages use their study ID in the filename. Register new chapters in `_quarto.yml`. Keep existing filenames and heading identifiers stable so published links continue to work.

From this directory, run `quarto render` with Quarto, Rscript and Python 3 available. The pre-render step regenerates the versioned help excerpts and the stellar visibility results table. The post-render step preserves links from the former single-page edition. The prebuilt edition includes the archived excerpts and table.

To check the function catalogue against an installed version of skyscapeR, run `Rscript scripts/check-coverage.R`.

## GitHub Pages

Copy the contents of this project's rendered `docs/` into the repository's `docs/manual/` directory. The repository's `docs/index.html` is the vignette. See the upload instructions distributed with this edition.
