# skyscapeR technical manual

Open `docs/index.html` for the chapter-based edition. Run `quarto render` from this directory to rebuild it. The pre-render script converts the versioned R help pages and calculates the heliacal-rising summary from its archived data.

For the single-page edition currently used on GitHub Pages, render the book first and run `python scripts/standalone/build.py /path/to/manual/index.html` with Python and lxml installed. The output embeds its figures, help pages and study downloads. The hosted vignette remains at the site root; the technical manual is in the `manual/` subdirectory.

The methods guide organises the book by archaeological question. The function reference gives the implementation files, formal arguments and links to individual help pages. The validation register includes H01, H02, A02, the sky-rendering studies and the validation programme. Some study links lead to access-controlled review copies; a public release should use archival links available to its readers.

Run `Rscript scripts/check-coverage.R` with skyscapeR 2.0.9004 installed to check the catalogue against the package namespace, formal argument lists and calls in the vignette. Scientific illustrations, their data sources and reproduction instructions are under `figures/` and `scripts/illustrations/`.
