# Archived sky-rendering checks

These five R scripts accompanied the 2.0.9004 package. They check twilight, colour, zoom and shared lighting. They use internal functions and should be interpreted against the development build used for each study.

Run scripts in a fresh R session with skyscapeR and its astronomy dependencies installed. For the twilight event comparison, run `events.R` from the `twilight/` directory before `compare.R`, which reads the resulting `events.csv`. The twilight holdout script can run separately.
