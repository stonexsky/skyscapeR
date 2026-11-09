# Render the full help from the versioned package sources.
write_help <- function(source, output) {
  help_body <- capture.output(tools::Rd2HTML(
    source, standalone = FALSE, no_links = TRUE, texmath = "text"
  ))
  writeLines(c("```{=html}", help_body, "```"), output)
}
write_help("source/star.phases.Rd", "_help.qmd")
write_help("source/createHor.Rd", "_horizon-help.qmd")

# The archived study stores signed offsets: positive is a later morning.
study <- read.csv("study/results.csv", na.strings = c("", "NA"))
main <- subset(study, sample == "main")
events <- subset(main, is.finite(calibrated_days))
offset <- events$calibrated_days
max_case <- events[which.max(abs(offset)), ]
unmatched <- subset(main, !is.finite(calibrated_days))
stopifnot(nrow(main) == 80L, nrow(unmatched) == 10L, all(unmatched$star == "Capella"))
lines <- c(
  "| Measure | Result |", "|:---|---:|",
  sprintf("| Comparable dates | %d |", nrow(events)),
  sprintf("| Absolute difference of three days or less | %d of %d |", sum(abs(offset) <= 3), length(offset)),
  sprintf("| Absolute difference of seven days or less | %d of %d |", sum(abs(offset) <= 7), length(offset)),
  sprintf("| Median signed difference | %+g days |", median(offset)),
  sprintf("| Largest absolute difference | %g days |", abs(max_case$calibrated_days)),
  "",
  sprintf("The largest difference occurred for %s at %s in %s. The comparator placed its first qualifying morning %g days later than the phase date.",
          max_case$star, max_case$site, max_case$epoch, max_case$calibrated_days)
)
writeLines(lines, "_results.qmd")
