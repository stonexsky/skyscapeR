# Normalize the container without flattening independently dated groups.
.sky_groups <- function(obj) {
  if (is.null(obj)) return(list())
  if (inherits(obj,"skyscapeR.object")) obj <- list(obj)
  if (!is.list(obj)) stop("obj must be a skyObject or a list of skyObject objects.")
  for (i in seq_along(obj)) {
    x <- obj[[i]]
    bad <- function(message) stop("Invalid skyObject group ",i,": ",message,call.=FALSE)
    if (!inherits(x,"skyscapeR.object")) bad("expected a skyObject (nested lists are not supported).")
    if (!is.numeric(x$epoch) || !length(x$epoch) %in% 1:2 || any(!is.finite(x$epoch))) bad("supply one epoch or two endpoints.")
    if (!is.matrix(x$decs) || !is.numeric(x$decs) || nrow(x$decs)!=length(x$epoch) ||
        !ncol(x$decs) || any(!is.finite(x$decs)) || any(abs(x$decs)>90)) bad("invalid declination matrix.")
    if (length(x$n)!=1L || is.na(x$n) || x$n!=ncol(x$decs)) bad("target count does not match declinations.")
    for (field in c("col","lty","lwd")) {
      value <- x[[field]]
      if (is.null(value)) value <- switch(field,col="red",lty=1,lwd=1)
      if (length(value)==1L) value <- rep(value,x$n)
      if (length(value)!=x$n) bad(paste("invalid",field,"settings."))
      x[[field]] <- value
    }
    obj[[i]] <- x
  }
  obj
}
