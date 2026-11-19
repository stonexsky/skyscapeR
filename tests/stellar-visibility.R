library(skyscapeR)
if (tryCatch({star("Sirius",2000); TRUE}, error=function(e) FALSE)) {
  loc <- c(36,14,0)
  physical <- star.phases("Sirius",2000,loc)
  passage <- star.phases("Sirius",2000,loc,method="purrington")
  e <- physical$metadata$events
  other <- !e$event %in% c("heliacal rising","acronycal setting")
  old <- passage$metadata$events
  stopifnot(nrow(physical$metadata$visibility.calendar)==366L,
            identical(physical$data,passage$data),
            all(e$method[other]=="purrington"), all(e$method[!other]=="schaefer"),
            identical(e$date.julian[other],old$date.julian[old$event %in% e$event[other]]),
            identical(skyscapeR:::.phase_plot_events(physical)$day,e$model.day))
  hidden <- star.phases("Sirius",2000,loc,alt.hor=89)
  stopifnot(!nrow(hidden$metadata$events),!any(hidden$metadata$visibility.calendar$morning))
}
