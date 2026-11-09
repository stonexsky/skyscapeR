#' Represent measured declinations as probability densities
#' @description Represent the uncertainty of directly supplied declinations so they can be summed or
#' tested against a declination model.
#' @details Use this constructor when both the measurement and its uncertainty
#' are already expressed in declination. All six families have the same parameter
#' meaning as in \code{\link{az.pdf}}. The doorway examples in 'Choosing a
#' measurement model' illustrate their shapes. When starting from measured bearings
#' and site horizons, use \code{\link{coordtrans}} to propagate the bearing model.
#' Direct declinations are linear. A triangle must satisfy lower <= dec <= upper.
#' The sampled curve is clipped to physical declinations from -90 to +90 degrees
#' and trimmed at \code{.cutoff}. The retained area can consequently be less than one.
#' @inheritSection az.pdf Choosing a measurement model
#' @inheritParams az.pdf
#' @param dec Measured declinations between -90 and +90 degrees. For triangular, the modal declination.
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
#' dec.pdf(-18, pdf="triangular", lower=-19, upper=-15)
#' dec.pdf(-18, 2, pdf="uniform_normal", survey.sd=0.2)
dec.pdf <- function(dec, unc=NULL, pdf="normal", name, .res=0.01, .cutoff=1e-4, verbose=TRUE,
                    taper=NULL, survey.sd=0.5, lower=NULL, upper=NULL) {
  if (!is.numeric(dec) || !length(dec) || any(!is.finite(dec) | abs(dec) > 90))
    stop("dec must contain finite declinations between -90 and 90 degrees.")
  if (missing(name)) name <- paste("Measurement", seq_along(dec))
  family <- .measurement_family(pdf)
  if (missing(survey.sd) && family != "uniform_normal") survey.sd <- NULL
  spec <- .measurement_spec(dec, "dec", unc, family, taper, survey.sd, lower, upper)
  .measurement_pdf(dec, "dec", spec$unc, pdf, spec$shape, name, .res, .cutoff, verbose)
}

.measurement_unc <- function(unc, n) {
  if (!is.numeric(unc) || !length(unc) %in% c(1L,n) ||
      any(!is.finite(unc) | unc<=0))
    stop("unc must be positive and finite, with one value or one per observation.")
  rep_len(unc, n)
}
