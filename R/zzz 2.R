.onLoad <- function(libname, pkgname) {
  .configure_ephemeris(lib.loc=libname)
}

.onAttach <- function(...) {
  label <- paste("  skyscapeR", utils::packageVersion("skyscapeR"))
  packageStartupMessage(paste(c(
    "         *",
    "*              *                     *",
    "",
    "                      *",
    "  Welcome to",
    sprintf("%-36s*", label),
    "                        *"
  ), collapse="\n"))
  if (interactive() && !.ephemeris_state$available && !.ephemeris_state$notified) {
    packageStartupMessage("For solar, lunar, planetary and stellar positions, install swephRdata: ", .ephemeris_install)
    .ephemeris_state$notified <- TRUE
  }
}
