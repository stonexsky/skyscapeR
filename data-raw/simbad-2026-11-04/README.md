# Stellar supplement, retrieved 2026-11-04

Run `python3 data-raw/simbad-2026-11-04/append.py` from the package root (or use
an absolute script path). Python's standard library is sufficient. The script
uses the archived responses, performs no downloads, and is idempotent. It writes
`inst/ephemeris/sefstars.txt.xz` and the installed provenance CSV.

SIMBAD TAP endpoint: https://simbad.cfa.harvard.edu/simbad/sim-tap/sync
Parameters: REQUEST=doQuery, LANG=ADQL, FORMAT=json, QUERY=<saved ADQL>.
The base query returns ICRS/J2000 astrometry and bibliographic references.
Separate queries supply V photometry, aliases, PM and parallax measurement
records, and the empty radial-velocity measurement result for the five objects
without a basic radial velocity. `plx29.json.adql` is the parallax query.
The exact VizieR request and catalogue metadata are preserved in the TSV header.
Only Hipparcos Johnson Vmag is imported from that table, not its J1991.25 positions.

The input records and the derived supplement retain the SIMBAD ODbL 1.0 notice:
https://opendatacommons.org/licenses/odbl/1-0/
See `inst/COPYRIGHTS` for acknowledgements, scientific sources and source-data
limitations. Unknown radial velocity is represented by a documented numeric
zero in sefstars because the engine needs a number; this is not a measurement.

No original numeric catalogue values are changed. HIP 63948 links to the
already-present 39 Com A entry rather than adding a second copy of that object.
