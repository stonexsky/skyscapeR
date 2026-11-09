#' Represent measured declinations as probability densities
#' @description Represent the uncertainty of directly supplied declinations so they can be summed or
#' tested against a declination model.
#' @details A declination estimate already incorporates the sightline geometry. This constructor
#' works directly in that coordinate. For each centre it evaluates a normal density with
#' the given standard deviation or a uniform density with the given half-width. A single
#' uncertainty is shared, or one can be supplied per observation. Save the original
#' uncertainty interpretation when comparing results from different surveys.
#'
#' Each grid is restricted to physical declinations and trimmed at .cutoff. The retained curve is not renormalised after truncation. Use \code{\link{spd}} to combine measurements.
#'
#' Silva (2020) sets out the probability-distribution framework for uncertain
#' orientations. Here the supplied measurements and uncertainties are already
#' expressed in declination.
#' @param pdf Distribution family, "normal" or "uniform", with the same abbreviations as \code{\link{az.pdf}}.
#' @param dec Measured declinations in degrees, between -90 and 90 inclusive.
#' @param unc Positive uncertainties in degrees of declination. Supply one value or one per measurement. For normal densities this is the standard deviation. For uniform densities it is the half-width.
#' @param name Optional measurement labels in input order.
#' @param verbose Logical. Report completion.
#' @param .cutoff Minimum retained density height. Default 0.0001.
#' @param .res Grid spacing in degrees. Default 0.01.
#' @return A skyscapeR.pdf list with one x/y data frame per observation. Metadata record coord = "dec", the input measurements, family, uncertainties, labels and numerical settings.
#' @section Validation:
#' In D01, sampled declination densities agreed with analytic references. Build 2.0.9006 retained the corresponding declination centres in all four subset checks. See \href{https://stonexsky.github.io/skyscapeR/manual/#d01-density-checks}{D01 validation in the technical manual}.
#' @export
#' @references Silva, F. (2020). A probabilistic framework and significance test for the analysis of structural orientations in skyscape archaeology. Journal of Archaeological Science 118, 105138. \doi{10.1016/j.jas.2020.105138}.
#' @seealso \code{\link{az.pdf}}, \code{\link{coordtrans}}, \code{\link{spd}}, \code{\link{modelTest}}
#' @examples
#' d <- dec.pdf(dec = c(-18, 12), unc = c(0.8, 1.4))
#' head(d$data[[1]])
dec.pdf <- function(dec, unc, pdf="normal", name, .res=0.01, .cutoff=1e-4, verbose=TRUE) {
  if (!is.numeric(dec) || !length(dec) || any(!is.finite(dec) | abs(dec) > 90))
    stop("dec must contain finite declinations between -90 and 90 degrees.")
  unc <- .measurement_unc(unc, length(dec))
  if (missing(name)) name <- paste("Measurement", seq_along(dec))
  if (length(name) != length(dec)) stop("Supply one name per declination.")
  if (!is.numeric(.res) || length(.res)!=1L || !is.finite(.res) || .res<=0 ||
      !is.numeric(.cutoff) || length(.cutoff)!=1L || !is.finite(.cutoff) || .cutoff<0)
    stop("Supply a positive finite spacing and a nonnegative finite cutoff.")
  if (length(pdf)!=1L || !is.character(pdf) ||
      !tolower(pdf) %in% c("n", "normal", "u", "uniform"))
    stop("Choose a normal or uniform density.")
  normal <- tolower(pdf) %in% c("n","normal")
  x <- seq(-90, 90, by=.res)
  data <- lapply(seq_along(dec), function(i) {
    y <- if (normal) stats::dnorm(x, dec[i], unc[i]) else
      stats::dunif(x, dec[i]-unc[i], dec[i]+unc[i])
    keep <- y >= .cutoff
    if (sum(keep)<2L || !any(y[keep]>0))
      stop("Too few density samples for measurement ", i, ". Reduce .res or .cutoff.")
    data.frame(x=x[keep], y=y[keep])
  })
  metadata <- list(coord="dec", pdf=pdf, name=name, dec=dec, unc=unc,
                   param=list(.cutoff=.cutoff, .res=.res))
  if (verbose) .sky_message("Created ", length(data), " declination densities.")
  structure(list(metadata=metadata, data=data), class="skyscapeR.pdf")
}

.measurement_unc <- function(unc, n) {
  if (!is.numeric(unc) || !length(unc) %in% c(1L,n) ||
      any(!is.finite(unc) | unc<=0))
    stop("unc must be positive and finite, with one value or one per observation.")
  rep_len(unc, n)
}
