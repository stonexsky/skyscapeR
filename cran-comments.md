# skyscapeR submission comments (working draft)

## Submission history

I request restoration of skyscapeR to CRAN. CRAN's incoming check records
that version 1.0.0 was archived on 24 October 2025 for a policy violation
concerning Internet access. The 10 October message requested corrections by
24 October. I left Bournemouth University in 2025 and had not updated my
maintainer email address promptly. I remain the maintainer. Please use
fsilva@stonexsky.com for future correspondence.

Installation and package loading require no downloads. Users who need
ephemeris files can install swephRdata explicitly or set SE_EPHE_PATH to
existing files. HeyWhatsThat requests are made by explicit user calls and
use HTTPS with timeouts. Network-dependent examples are skipped during checks.

## Changes

This release revises orientation-distribution tests and the propagation of
measurement uncertainty. It adds terrain and sky reconstruction tools and
updates astronomical event calculations. NEWS.md summarises changes affecting
published analyses. The README contains migration instructions.

The package unpacks its stellar catalogue in session-temporary storage.
Bundled data credits and terms are documented in COPYRIGHTS.

## R CMD check results

[Replace with the release version, R versions and platforms, and counts of
ERRORs, WARNINGs and NOTEs from checks on the final source tarball. Include
the R-devel --as-cran result and explain each remaining NOTE.]

## Reverse dependencies

A search of the CRAN source package index retrieved on 20 November 2026
found no packages listing skyscapeR in Depends, Imports, LinkingTo,
Suggests or Enhances. [Refresh at submission if the index has changed.]
