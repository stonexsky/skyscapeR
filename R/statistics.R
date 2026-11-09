utils::globalVariables("hor")
#' @noRd
bernoulli <- function(n,p,r) {
  stats::dbinom(r, size=n, prob=p)
}


#' Calculate binomial probabilities for aligned structures
#' @description Assess how often a specified number of independent structures would face a target by
#' chance, given the chance for one structure.
#' @details Define a hit and the per-structure chance of a hit before counting the orientations. A
#' single-count probability answers how often exactly that count arises; the upper tail
#' answers how often a count at least that large arises.
#'
#' The binomial model assumes independent structures with the same
#' target probability \code{p}. The upper-tail probability is a one-sided
#' p-value under this model. Calculations use \code{stats::dbinom} and
#' \code{stats::pbinom}, avoiding factorial overflow and subtraction from 1.
#' At least zero hits has probability 1. With zero structures, zero hits has
#' probability 1 for either type. Inputs are scalar; invalid inputs error.
#' @param n Total number of structures: a single finite nonnegative integer.
#' @param p Target probability: a single finite number between 0 and 1.
#' @param r Number of hits: a single finite integer between 0 and n.
#' @param type Either \code{"single"} for exactly r hits or \code{"tail"}
#' for at least r hits (the default).
#' @return A numeric probability between 0 and 1.
#' @section Validation:
#' Exact and upper-tail probabilities agreed with independent binomial calculations for the tested counts up to 1,000, including zero hits. See \href{https://stonexsky.github.io/skyscapeR/manual/#s03-scoring-checks}{S03 validation in the technical manual}.
#' @export
#' @references Ruggles, C (1999) \emph{Astronomy in Prehistoric Britain and Ireland}. Yale University Press.
#' @seealso \code{\link{randomTest}}, \code{\link{modelTest}}
#' @examples
#' # At least 10 out of 30 structures, with targets covering 20% of the horizon:
#' bernoulli.trial(30, 0.2, 10)
#' # Exactly 10 hits:
#' bernoulli.trial(30, 0.2, 10, type="single")
#' bernoulli.trial(30, 0.2, 0) # Always 1
#' bernoulli.trial(1000, 0.2, 250) # No factorial overflow
bernoulli.trial <- function(n,p,r, type='tail') {
  valid.scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x)
  }
  if (!valid.scalar(n) || n < 0 || n != floor(n)) {
    stop('n must be a single finite nonnegative integer.')
  }
  if (!valid.scalar(p) || p < 0 || p > 1) {
    stop('p must be a single finite probability between 0 and 1.')
  }
  if (!valid.scalar(r) || r < 0 || r > n || r != floor(r)) {
    stop('r must be a single finite integer between 0 and n.')
  }
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c('single','tail')) {
    stop('type must be "single" or "tail".')
  }
  if (type == 'single') return(bernoulli(n,p,r))
  stats::pbinom(r - 1, size=n, prob=p, lower.tail=FALSE)
}


#' Represent measured azimuths as probability densities
#' @description Represent uncertainty in a measured orientation before comparing it with other bearings
#' or transforming it to declination.
#' @details A field bearing defines a centre and an angular uncertainty. Model the uncertainty as a
#' normal density with the supplied standard deviation, or as a uniform density over the
#' centre plus or minus the supplied half-width. A scalar uncertainty is reused for each
#' measurement. The function evaluates the periodic density on a grid starting at zero
#' and ending before 360 degrees. Normal densities sum wrapped Gaussian contributions;
#' uniform densities count all copies within the inclusive half-width. The cutoff is
#' applied to the combined density at each bearing. The input measurements remain separate and truncated tails
#' are not renormalised. The supported families are normal and uniform.
#' @return A skyscapeR.pdf list. The data component contains one data frame per measurement, in input order, with azimuth x in degrees and density y. The metadata record the measurements, names, uncertainties, distribution type and numerical settings.
#' @param pdf Uncertainty distribution. Choose \code{"normal"} (default) or \code{"uniform"}.
#' @param az True bearings in degrees clockwise from north, one per observation.
#' @param unc Angular uncertainty in degrees, one value shared across observations or one per azimuth. For a normal distribution this is the standard deviation. For a uniform distribution it is the half-width on either side of the measurement.
#' @param name Optional labels in the same order as the measurements.
#' @param verbose Logical. Display progress or informational messages. Default \code{TRUE}.
#' @param .cutoff Minimum density retained on the output grid. Default \code{1e-5}.
#' @param .res Spacing of the output bearing grid in degrees. Default 0.01.
#' @section Validation:
#' Build 2.0.9006 passed all 173 D01 assertions, including 0.005-degree samples across north and the broad-width periodic density cases. See \href{https://stonexsky.github.io/skyscapeR/manual/#d01-density-checks}{D01 validation in the technical manual}.
#' @export
#' @references Silva, F (2020) A probabilistic framework and significance test for the analysis of structural orientations
#'  in skyscape archaeology \emph{Journal of Archaeological Science} 118, 105138. \doi{10.1016/j.jas.2020.105138}
#' @seealso \code{\link{coordtrans}}, \code{\link{spd}}, \code{\link{hpdi}}, \code{\link{randomTest}}, \code{\link{dec.pdf}}, \code{\link{modelTest}}
#' @examples
#' test <- az.pdf(az=c(87,93,90,110), unc=3)
#' plot(test)
az.pdf <- function(az, unc, pdf='normal', name, .res=0.01, .cutoff=1e-5, verbose=T) {

  if (!is.character(pdf) || length(pdf)!=1L || is.na(pdf) ||
      !tolower(pdf) %in% c("n","normal","u","uniform"))
    stop("Choose a normal or uniform density.")
  if (!is.numeric(az) || !length(az) || any(!is.finite(az)))
    stop("Supply finite azimuths.")
  n <- length(az)
  unc <- .measurement_unc(unc,n)
  if (!is.numeric(.res) || length(.res)!=1L || !is.finite(.res) || .res<=0 || .res>=360 ||
      !is.numeric(.cutoff) || length(.cutoff)!=1L || !is.finite(.cutoff) || .cutoff<0)
    stop("Supply a positive grid spacing below 360 degrees and a nonnegative cutoff.")
  if (!missing(name) && length(name)!=n) stop("Supply one name per azimuth.")
  normal <- tolower(pdf) %in% c("n","normal")
  x <- seq(0,360,by=.res)
  x <- x[x<360]
  out <- vector("list",n)
  if (.sky_progress(verbose) && n>1L) {
    pb <- txtProgressBar(max=n,style=3)
    on.exit(close(pb),add=TRUE)
  }
  for (k in seq_len(n)) {
    delta <- (x-(az[k] %% 360)+180) %% 360-180
    width <- unc[k]
    if (!normal) {
      # Count all copies of this bearing in the inclusive uniform interval.
      turns <- floor((width-delta)/360) + floor((width+delta)/360) + 1
      yy <- (turns/width)/2
    } else if (width<360) {
      yy <- stats::dnorm(delta,mean=0,sd=width)
      for (turn in seq_len(ceiling((9*width+180)/360))) {
        yy <- yy + stats::dnorm(delta+360*turn,mean=0,sd=width) + stats::dnorm(delta-360*turn,mean=0,sd=width)
      }
    } else {
      # At sigma >= one turn, Fourier harmonics beyond the first are below double precision.
      yy <- (1+2*exp(-0.5*(width*(pi/180))^2)*cos(delta*pi/180))/360
    }
    ind <- which(yy >= .cutoff)
    out[[k]] <- data.frame(x = x[ind], y = yy[ind])

    if (.sky_progress(verbose) && length(az)>1) { setTxtProgressBar(pb, k) }
  }

  if (missing(name)) { name <- paste('Measurement',seq(1,length(az))) }

  mtdta <- list(
    coord = 'az',
    pdf = pdf,
    name = name,
    az = az,
    unc = unc)
  mtdta$param <- list(.cutoff=.cutoff, .res=.res)
  out <- list(metadata = mtdta, data = out)
  class(out) <- 'skyscapeR.pdf'

  return(out)
}


#' @noRd
coordtrans_unc <- function(pdf, hor, refraction, atm, temp, xrange, .res=0.1, .prec=0.1,
                           .output.res=0.01, verbose=T) {

  if ((inherits(pdf, 'skyscapeR.pdf') & pdf$metadata$coord != 'az') & (!inherits(pdf, 'skyscapeR.spd'))) stop('Azimuthal probability density function(s) not recognized or not in correct format.')
  if (!inherits(hor, 'skyscapeR.horizon') & !inherits(hor, 'list')) stop('Horizon(s) not recognized or not in correct format.')

  valid.spacing <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  for (setting in c('.res', '.prec', '.output.res')) {
    if (!valid.spacing(get(setting))) stop(paste(setting, 'must be a single positive finite spacing in degrees.'))
  }
  if (.prec > 540) stop('.prec must not exceed the 540-degree integration domain.')
  if (!missing(xrange) && (!is.numeric(xrange) || length(xrange) != 2L ||
      any(!is.finite(xrange)) || xrange[1] >= xrange[2])) {
    stop('xrange must contain two finite increasing bounds.')
  }

  ## init parameters
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (is(pdf, 'skyscapeR.spd')) {
    n <- 1
    pdf$data <- list(pdf$data)
  } else {
    n <- length(pdf$metadata$az)
  }
  .cutoff <- pdf$metadata$param$.cutoff


  shared.horizon <- is(hor, "skyscapeR.horizon")
  if (shared.horizon && verbose) {
    .sky_message('Single horizon profile found. Using it for all measurements.')
  }

  xx <- seq(-90, 360+90, .prec)

  out <- list()
  if (.sky_progress(verbose) && n>1) { pb <- txtProgressBar(max = n, style = 3); on.exit(close(pb),add=TRUE) }

  for (k in 1:n) {
    horizon <- if (shared.horizon) hor else hor[[k]]
    curve <- .unwrap_north_pdf(pdf$data[[k]])
    azs <- curve$x
    grid <- xx
    if (!identical(azs,pdf$data[[k]]$x)) {
      grid <- seq(floor(min(azs)/.prec)-1, ceiling(max(azs)/.prec)+1)*.prec
    }

    profile <- .coord_horizon(horizon)
    query <- .profile_bearings(azs, profile$az)
    fhor <- approx(profile$az, profile$alt, xout=query)$y
    fhor.unc <- approx(profile$az, profile$alt.unc, xout=query)$y

    if (any(!is.finite(fhor)) || any(!is.finite(fhor.unc)) || any(fhor.unc < 0)) {
      stop(paste0('Distribution ', k, ': horizon altitude and uncertainty must be finite at every input azimuth; uncertainty must be nonnegative.'))
    }

    aux1 <- az2dec(azs, horizon, fhor-4*fhor.unc, refraction=refraction, atm=atm, temp=temp)
    aux2 <- az2dec(azs, horizon, fhor+4*fhor.unc, refraction=refraction, atm=atm, temp=temp)

    if (any(!is.finite(aux1)) || any(!is.finite(aux2))) {
      stop(paste0('Distribution ', k, ': transformed uncertainty bounds must be finite.'))
    }
    mean.dec <- (aux2+aux1)/2
    sd.dec <- abs(aux2-aux1)/8
    if (any(!is.finite(sd.dec)) || any(sd.dec == 0)) {
      stop(paste0('Distribution ', k, ': zero or nonfinite transformed uncertainty width. ',
                  'The Gaussian approximation cannot represent this case; check altitude uncertainties and transformation geometry. ',
                  'Do not add artificial uncertainty merely to bypass this error.'))
    }

    if (missing(xrange)) {
      min.dec <- min(aux1, aux2, na.rm=T); max.dec <- max(aux1, aux2, na.rm=T)
    } else {
      min.dec <- xrange[1]; max.dec <- xrange[2]
    }

    dec <- seq(min.dec, max.dec, by = .res)
    # Include the upper bound so interpolation never extrapolates the tail.
    if (tail(dec, 1L) < max.dec) dec <- c(dec, max.dec)

    # Interpolation is linear. Reuse the position and mixing fraction of
    # each integration point across all trial declinations.
    inside <- grid >= azs[1] & grid <= tail(azs,1)
    pts <- grid[inside]
    lo <- findInterval(pts,azs,all.inside=TRUE)
    mix <- (pts-azs[lo])/(azs[lo+1L]-azs[lo])
    az.density <- approx(azs,curve$y,xout=pts)$y
    steps <- diff(grid)
    weight <- ((c(0,steps)+c(steps,0))/2)[inside]*az.density
    dens <- vapply(dec,function(value) {
      d <- dnorm(value,mean.dec,sd.dec)
      sum(weight*(d[lo]+mix*(d[lo+1L]-d[lo])))
    },numeric(1))

    if (any(!is.finite(dens)) || any(dens < 0)) {
      stop(paste0('Distribution ', k, ': calculated density is nonfinite or negative; inspect uncertainty widths and numerical spacings.'))
    }

    output.grid <- seq(min.dec, max.dec, .output.res)
    if (length(output.grid) < 2L) stop('Requested output grid has fewer than two points; reduce .output.res or widen xrange.')
    dens <- approx(dec, dens, xout=output.grid, method='linear')
    dec <- dens$x; dens <- dens$y

    if (!is(pdf, 'skyscapeR.spd')) {
      area <- .grid_area(dec,dens)
      if (!is.finite(area) || area <= 0) {
        stop(paste0('Distribution ', k, ': no positive finite density area on the output grid; inspect xrange and numerical spacings.'))
      }
      dens <- dens/area
    }

    out[[k]] <- data.frame(x=dec, y=dens)
    if (.sky_progress(verbose) && n>1) { setTxtProgressBar(pb, k) }
  }

  mtdta <- list(
    coord = 'dec',
    name = pdf$metadata$name
  )
  mtdta$param <- list(refraction = refraction, atm=atm, temp=temp, .cutoff=.cutoff, .res=.res,
                       .prec=.prec, .output.res=.output.res)
  # Labels are descriptive only: never resolve them with get() or eval().
  # Capture no input objects or environments in the returned metadata.
  mtdta$format.version <- 2L
  mtdta$source.labels <- list(
    az.pdf = deparse(substitute(pdf)),
    horizon = deparse(substitute(hor))
  )
  out <- list(metadata = mtdta, data = out)
  if (is(pdf, 'skyscapeR.spd')) {
    class(out) <- 'skyscapeR.spd'
    out$data <- out$data[[1]]
  } else {
    class(out) <- 'skyscapeR.pdf'
  }
  return(out)
}

#' Transform azimuth densities into declination densities
#' @description Estimate which declinations a measured bearing could indicate after allowing for
#' uncertainty in both direction and horizon altitude.
#' @details A sightline has a distribution of possible bearings, and each bearing samples a
#' potentially different part of the horizon. For every sampled bearing the function
#' interpolates horizon altitude and uncertainty, maps an altitude range to declination and
#' approximates the conditional declination distribution by a Gaussian. It weights that
#' distribution by the azimuth density and integrates the product by the trapezoidal rule.
#' A shared horizon can serve several observations, or a supplied list can represent their
#' individual sites.
#'
#' Density values are calculated at declination spacing \code{.res}, then
#' linearly interpolated onto a grid with spacing \code{.output.res}. The upper
#' calculation bound is evaluated explicitly if the step does not land on it,
#' avoiding extrapolation at the tail. Linear interpolation preserves
#' nonnegative calculated densities and avoids cubic-spline overshoot.
#' Individual PDF outputs retain their existing area normalisation; SPD outputs
#' retain their existing unnormalised transformation behaviour. Coarse
#' calculation spacing can still miss narrow peaks; interpolation does not
#' replace a resolution-convergence check. Azimuth integration uses a separate
#' grid from -90 to 450 degrees with step \code{.prec}. For a distribution split across north, integration uses contiguous bearings around zero at the same spacing. The output grid starts
#' at the lower bound and ends at the last full step within the upper bound.
#' Refine the calculation grid to check that narrow features are resolved.
#' Set atmospheric refraction, pressure and temperature with \code{\link{skyscapeR.vars}} before transforming. The values used are recorded in the returned PDF or SPD.
#' All three spacings are stored in \code{metadata$param}.
#'
#' The Gaussian uncertainty approximation uses the midpoint of declinations
#' obtained at horizon altitude plus/minus four altitude standard deviations,
#' with standard deviation equal to the absolute endpoint separation divided
#' by eight. Reversed bounds are supported. Horizon uncertainties must be
#' finite and nonnegative, and the profile must cover all input azimuths.
#' Zero transformed widths (including zero altitude uncertainty or coincident
#' endpoint declinations) are unsupported and produce an error. Do not add
#' arbitrary uncertainty to bypass it: these cases need a different treatment.
#' Nonfinite transformed bounds/densities and PDF outputs with no positive
#' finite area also error. Positive but very narrow widths still require
#' convergence checks; no minimum uncertainty is imposed.
#'
#' \strong{Keeping the inputs.} Save the original azimuth distributions and horizons
#' in measurement order if you plan to draw the transformation curves or run a
#' declination test. The result stores the numerical distributions and settings.
#' The names in \code{metadata$source.labels} describe the inputs but cannot
#' reconstruct them. See \code{\link{plot.skyscapeR.pdf}} and \code{\link{randomTest}}.
#' @param pdf A \emph{skyscapeR.pdf} object created with \code{\link{az.pdf}} or a \emph{skyscapeR.spd} object created with \code{\link{spd}}
#' @param hor A shared \code{skyscapeR.horizon} or, for individual azimuth PDFs, a list of horizons in measurement order.
#' @param xrange Array of values (min and max) for SPD when transforming a \emph{skyscapeR.spd} object.
#' @param verbose Logical. Display progress or informational messages. Default \code{TRUE}.
#' @param .res Spacing of the declination calculation grid, in degrees. Default 0.1.
#' @param .prec Azimuth integration spacing, in degrees. Default 0.1.
#' Previously ignored; smaller values increase integration work. Must not exceed 540.
#' @param .output.res Spacing of the returned declination grid, in degrees.
#' Default 0.01. Coarser spacing reduces output size. All spacings must be
#' positive finite scalars, and the returned grid must contain at least two points.
#' @return A \code{skyscapeR.pdf} object for PDF input, or a \code{skyscapeR.spd}
#' object for SPD input. Both contain \code{metadata} and \code{data}. PDF data
#' are a list of data frames with declination \code{x} and density \code{y};
#' SPD data are a single such data frame. Metadata include the coordinate type,
#' transformation settings, format version and descriptive source labels.
#' @export
#' @references Silva, F (2020) A probabilistic framework and significance test for the analysis of structural orientations
#'  in skyscape archaeology \emph{Journal of Archaeological Science} 118, 105138. \doi{10.1016/j.jas.2020.105138}
#' @seealso \code{\link{az.pdf}}, \code{\link{spd}}, \code{\link{hpdi}}, \code{\link{randomTest}}, \code{\link{dec.pdf}}, \code{\link{modelTest}}
#' @examples
#' Az <- az.pdf(az=c(87,93,90,110), unc=3)
#' hor <- createHor(az=c(0,360), alt=c(0,0), loc=c(35,-8,25)) # flat horizon
#' Dec <- coordtrans(Az, hor)
#' head(Dec$data[[1]])
#' Dec$metadata$source.labels
#' # Reduce stored output size without changing the calculation spacing:
#' Dec.coarse <- coordtrans(Az, hor, .output.res=0.1)
#' Dec.coarse$metadata$param
coordtrans <- function(pdf, hor, xrange, .res=.1, .prec=.1, .output.res=.01, verbose=TRUE) {
  result <- coordtrans_unc(pdf, hor, xrange=xrange, .res=.res, .prec=.prec,
                           .output.res=.output.res, verbose=verbose)
  result$metadata$source.labels <- list(az.pdf=deparse(substitute(pdf)), horizon=deparse(substitute(hor)))
  result
}


# Integrate a piecewise linear density on its ordered sampling grid.
.grid_area <- function(x, y) {
  0.5 * sum(diff(x) * (head(y, -1L) + tail(y, -1L)))
}

#' @noRd
.normalise_spd_area <- function(xx, spd) {
  if (length(xx) < 2L || any(!is.finite(spd)) || any(spd < 0)) {
    stop('Cannot normalise SPD: need at least two grid points and finite nonnegative densities.')
  }
  area <- sum(diff(xx) * (head(spd, -1L) + tail(spd, -1L))/2)
  if (!is.finite(area) || area <= 0) {
    stop('Cannot normalise SPD: area on the returned grid must be positive and finite.')
  }
  spd/area
}

#' Sum orientation probability densities
#' @description Combine uncertain orientations to show where a group of structures concentrates in
#' azimuth or declination.
#' @details A summed probability density gives each input measurement one density curve. The
#' function linearly interpolates these curves onto one angular grid, assigns zero outside
#' each curve's sampled support and adds their heights point by point. With an
#' automatically chosen range it zeroes summed heights below the cutoff. With an explicit
#' range it retains the calculated heights across that range.
#'
#' Normalisation takes place after interpolation and any cutoff, over
#' the returned grid. An explicit restricted \code{xrange} therefore normalises
#' the retained portion only. Zero-area, nonfinite or negative densities and
#' grids with fewer than two points cannot be normalised and produce an error.
#' Earlier versions divided by the sum of sampled heights. Normalised outputs
#' now have a different scale; on a regular grid with zero endpoint densities,
#' new heights equal old heights divided by \code{.res}. Existing saved objects
#' are not converted. Changing resolution can still affect interpolation and
#' truncation accuracy, but the returned trapezoidal area is one.
#' @param pdf A \emph{skyscapeR.pdf} object created with either \code{\link{az.pdf}} or \code{\link{coordtrans}}
#' @param normalise If TRUE, scale the returned density to unit area using
#' trapezoidal integration. Default FALSE leaves the summed density unnormalised.
#' @param xrange Array of values (min and max) for SPD if different from range of \emph{pdf}. If given,
#' it overrides \emph{.cutoff}
#' @param .cutoff When the range is chosen automatically, set summed densities below this threshold to zero. Default \code{1e-5}.
#' @param .res Resolution with which to output SPD. Default is 0.01 degrees.
#' @section Validation:
#' The tested sums agreed with independently constructed density sums, and the normalised curves had unit trapezoidal area. See \href{https://stonexsky.github.io/skyscapeR/manual/#d01-density-checks}{D01 validation in the technical manual}.
#' @export
#' @references Silva, F (2020) A probabilistic framework and significance test for the analysis of structural orientations
#'  in skyscape archaeology \emph{Journal of Archaeological Science} 118, 105138. \doi{10.1016/j.jas.2020.105138}
#' @return A \code{skyscapeR.spd} list with \code{metadata} and a \code{data}
#' data frame containing coordinate \code{x} (degrees) and density \code{y}.
#' Normalised densities have units of inverse degrees.
#' @seealso \code{\link{az.pdf}}, \code{\link{coordtrans}}, \code{\link{hpdi}}, \code{\link{randomTest}}, \code{\link{dec.pdf}}, \code{\link{modelTest}}, \code{\link{plot.skyscapeR.spd}}
#' @examples
#' # SPD of azimuths
#' Az <- az.pdf(az=c(87,93,90,110), unc=3)
#' s1 <- spd(Az)
#' plot(s1)
#' s.norm <- spd(Az, normalise=TRUE)
#' with(s.norm$data, sum(diff(x) * (head(y,-1) + tail(y,-1))/2)) # 1
#'
#' # SPD of declinations
#' hor <- createHor(az=c(0,360), alt=c(0,0), loc=c(35,-8,25)) # flat horizon
#' Dec <- coordtrans(Az, hor)
#' s2 <- spd(Dec)
#' plot(s2)
spd <- function(pdf, xrange, normalise = F, .res = 0.01, .cutoff = 1e-5) {
  if (!inherits(pdf, 'skyscapeR.pdf')) { stop('Please provide a valid skyscapeR.pdf object.') }

  aux <- pdf$data
  x <- lapply(aux, "[[", 'x')
  y <- lapply(aux, "[[", 'y')
  if (missing(xrange)) {
    xrange <- range(x)
    cutoff <- T
  } else { cutoff <- F }

  xx <- seq(xrange[1], xrange[2], .res)
  spd <- rep(0, length(xx))
  for (i in 1:NROW(x)) {
    curve <- aux[[i]]
    query <- xx
    if (identical(pdf$metadata$coord,"az")) {
      # Interpolate periodically, with zero density in the discarded tails.
      step <- pdf$metadata$param$.res
      if (is.null(step)) step <- min(diff(curve$x))
      turn <- rbind(transform(tail(curve,1),x=x-360),curve,transform(head(curve,1),x=x+360))
      gaps <- which(diff(turn$x)>step*(1+1e-7))
      if (length(gaps)) {
        left <- turn$x[gaps]+step
        right <- turn$x[gaps+1L]-step
        missing <- unique(round(c(left,right),12))
        turn <- rbind(turn,data.frame(x=missing,y=0))
      }
      curve <- turn
      query <- xx %% 360
    }
    spd <- spd + approx(curve$x, curve$y, xout=query, yleft=0, yright=0)$y
  }

  if (cutoff) { spd[spd < .cutoff] <- 0 }

  if (normalise) spd <- .normalise_spd_area(xx, spd)
  out <- c()
  out$metadata <- c()
  out$metadata$coord <- pdf$metadata$coord
  out$metadata$xrange <- xrange
  atmospheric <- intersect(c("refraction","atm","temp"), names(pdf$metadata$param))
  if (length(atmospheric)) out$metadata$param <- pdf$metadata$param[atmospheric]
  if (cutoff) { out$metadata$param$.cutoff <- .cutoff }
  out$metadata$param$.res <- .res
  out$data <- data.frame(x = xx, y = spd)

  class(out) <- 'skyscapeR.spd'
  return(out)
}



#' @noRd
randomTest_unc <- function(pdf, az, hor, nsims=1000, conf=.95, tails=2, normalise=F, .res=0.1,
                           ncores=NULL, verbose=T) {
  if (!inherits(pdf, 'skyscapeR.pdf') ||
      !identical(pdf$metadata$coord, 'az') && !identical(pdf$metadata$coord, 'dec')) {
    stop("pdf must be an azimuth or declination skyscapeR.pdf object.")
  }
  if (pdf$metadata$coord == 'az' && !missing(az)) {
    stop("az is only needed when pdf contains declination distributions; supply azimuth data as pdf.")
  }
  if (pdf$metadata$coord == 'dec') {
    if (missing(az) || missing(hor)) {
      stop("For declination input, supply both az (the original azimuth PDF) and hor (the original horizon or horizons).")
    }
    if (!inherits(az, 'skyscapeR.pdf') || !identical(az$metadata$coord, 'az') ||
        length(az$data) != length(pdf$data)) {
      stop("az must contain the original azimuth distributions in the same number and order as pdf.")
    }
    saved <- pdf$metadata$param
    if (!is.logical(saved$refraction) || length(saved$refraction) != 1L || is.na(saved$refraction) ||
        !is.numeric(saved$atm) || length(saved$atm) != 1L || !is.finite(saved$atm) ||
        !is.numeric(saved$temp) || length(saved$temp) != 1L || !is.finite(saved$temp)) {
      stop("pdf lacks valid stored transformation settings; regenerate it with coordtrans().")
    }
    valid.horizon <- function(h) {
      inherits(h, 'skyscapeR.horizon') && is.numeric(h$metadata$georef) &&
        length(h$metadata$georef) == 3L && all(is.finite(h$metadata$georef)) &&
        is.data.frame(h$data) && all(c('az','alt','alt.unc') %in% names(h$data)) &&
        nrow(h$data) >= 2L &&
        all(vapply(h$data[c('az','alt','alt.unc')], is.numeric, logical(1))) &&
        all(is.finite(as.matrix(h$data[c('az','alt','alt.unc')])))
    }
    if (!valid.horizon(hor) &&
        !(is.list(hor) && length(hor) == length(pdf$data) &&
          all(vapply(hor, valid.horizon, logical(1))))) {
      stop("hor must be a valid shared horizon or an ordered list matching pdf.")
    }
  }
  source <- if (pdf$metadata$coord == 'dec') az else pdf
  n <- length(source$data)
  if (!n || length(source$metadata$name) != n || length(source$metadata$unc) != n ||
      !is.numeric(source$metadata$unc) || any(!is.finite(source$metadata$unc)) ||
      any(source$metadata$unc <= 0) ||
      !is.character(source$metadata$pdf) || length(source$metadata$pdf) != 1L ||
      !source$metadata$pdf %in% c('N','Normal','n','normal','U','Uniform','u','uniform') ||
      !is.numeric(source$metadata$param$.res) || length(source$metadata$param$.res) != 1L ||
      !is.finite(source$metadata$param$.res) || source$metadata$param$.res <= 0) {
    stop("Source azimuth PDF must have valid distribution type, uncertainties, names and resolution.")
  }
  if (!is.null(ncores) && (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores < 1 || ncores != floor(ncores))) stop("ncores must be a positive integer.")
  if (!is.numeric(nsims) || length(nsims)!=1L || !is.finite(nsims) ||
      nsims < 2 || nsims != floor(nsims)) stop("nsims must be an integer of at least two.")
  if (length(tails)!=1L || !tails %in% 1:2) stop("Use tails=1 or 2.")
  transform.res <- .res
  transform.prec <- 0.1
  transform.output.res <- 0.01
  if (pdf$metadata$coord == 'dec' && !is.null(saved$.output.res)) {
    settings <- saved[c('.res','.prec','.output.res')]
    if (!all(vapply(settings, function(x) is.numeric(x) && length(x)==1L &&
        is.finite(x) && x>0, logical(1)))) {
      stop('Invalid stored transformation spacings; regenerate pdf with coordtrans().')
    }
    transform.res <- saved$.res
    transform.prec <- saved$.prec
    transform.output.res <- saved$.output.res
  }
  quick <- F
  if (pdf$metadata$coord == 'dec' | !missing(hor) ) {
    if (verbose) .sky_message('Testing uniformly random azimuths in declination (',nsims,' simulations).')
    coord <- 'dec'
    xrange <- c(-90,90)
    if (pdf$metadata$coord == 'dec') {
      refraction <- pdf$metadata$param$refraction
      atm <- pdf$metadata$param$atm
      temp <- pdf$metadata$param$temp
      empirical <- spd(pdf, xrange=xrange, normalise = normalise, .res=.res)
    } else {
      az <- pdf
      if (!inherits(hor,"skyscapeR.horizon")) { stop('Please include a valid skyscapeR.horizon object.') }
      if (verbose) { .sky_message('Single horizon profile found. Using it for all measurements.') }
      refraction <- skyscapeR.env$refraction
      atm <- skyscapeR.env$atm
      temp <- skyscapeR.env$temp
      quick <- T
      empirical <- spd(pdf, .res=.res)
      empirical <- coordtrans_unc(empirical, hor, xrange=xrange, refraction=refraction, atm=atm, temp=temp, verbose=F, .res=.res)
      if (normalise) empirical$data$y <- .normalise_spd_area(empirical$data$x, empirical$data$y)
    }
  } else {
    if (verbose) .sky_message('Testing uniformly random azimuths in azimuth (',nsims,' simulations).')
    coord <- 'az'
    xrange <- c(0,360)
    az <- pdf
    empirical <- spd(pdf, normalise = normalise, xrange=xrange, .res=.res)
  }

  grid <- seq(0, 360, source$metadata$param$.res)
  grid <- grid[grid<360]
  settings <- list(family=source$metadata$pdf, unc=source$metadata$unc,
                   coord=coord, quick=quick, xrange=xrange, resolution=.res,
                   normalise=normalise)
  if (coord == "dec") {
    if (verbose) .sky_message(sprintf("Using atmospheric refraction %s, pressure %g mbar and temperature %g C.",refraction,atm,temp))
    settings$transform <- list(hor=hor, refraction=refraction, atm=atm, temp=temp,
                               .res=transform.res)
    if (!quick) {
      settings$transform$.prec <- transform.prec
      settings$transform$.output.res <- transform.output.res
    }
  }
  result <- .run_simulations(empirical$data, nsims, ncores, n,
    function(i) sample(grid, n, replace=TRUE), .random_curve, settings, conf, tails, verbose)

  ## output
  out <- c()
  out$metadata$coord <- coord
  if (coord=="dec") out$metadata$atmosphere <- settings$transform[c("refraction","atm","temp")]
  out$metadata$nsims <- nsims
  out$metadata$inference <- result$method
  out$metadata$local.inference <- result$local.method
  out$metadata$conf <- conf
  out$metadata$tails <- tails
  out$metadata$normalise <- normalise
  out$metadata$global.pval <- result$global
  out$metadata$local.pval <- result$local

  out$data$empirical <- empirical$data
  out$data$null.hyp <- result$envelope
  out$data$constant.departures <- result$departures
  class(out) <- 'skyscapeR.sigTest'

  return(out)
}


#' Test orientations against the null hypothesis of random azimuths
#' @description Use a Monte Carlo significance test to assess whether the measured orientations of a group of structures are compatible with
#' uniformly random azimuths, taking their uncertainties and, for a declination comparison,
#' their horizons into account.
#' @details Following Silva (2020), the test combines the observed measurement densities into a summed probability density (SPD).
#' For each simulated dataset, the function draws one azimuth per measurement with replacement
#' from the uniform grid from 0 to 360 degrees at the source azimuth spacing. It applies the
#' source density family and measurement uncertainties and constructs a comparable SPD.
#' For a declination test, random azimuths are transformed using the site geometry; the
#' resulting null distribution in declination depends on latitude and horizon.
#'
#' A new transformation takes its atmosphere from \code{\link{skyscapeR.vars}}. A saved declination PDF supplies its original settings. The fit records these in metadata$atmosphere for declination tests.
#'
#' \strong{Input coordinates.} With azimuth input alone, all curves are compared in azimuth. With azimuth input and one
#' shared horizon, the function sums the azimuth densities before transforming the SPD.
#' With a transformed declination PDF, supply its original \code{az} and \code{hor} in measurement
#' order. Simulated measurements are then transformed individually before summation, using the
#' stored refraction and atmosphere settings. Newer PDFs store the
#' transformation grid spacings; older PDFs use the legacy defaults. The comparison grid uses
#' \code{.res}. The two transformation routes can give different results.
#' If \code{normalise=TRUE}, each final curve has unit trapezoidal area. Each parallel worker
#' receives its own copy of the required horizons.
#'
#' \strong{Scoring the curves.} At each angle, the observed summed density and all simulated summed densities jointly define a pointwise envelope at confidence
#' level \code{conf}. An upper-tail test scores excesses; a two-tailed test also scores deficits,
#' using equal tail probabilities for the bounds. For every curve, each departure beyond a bound
#' is divided by the standard deviation of all B+1 curves at that angle. Ties at a bound contribute zero.
#' The global statistic sums the departure magnitudes across the grid. A small global p-value
#' means that few simulated datasets depart from the envelope as much as the observed dataset.
#'
#' \strong{Interpreting local intervals.} Local results identify contiguous clusters of departures on the stored linear grid.
#' At least two consecutive points must cross the same bound and span at least 0.5 degrees,
#' inclusive. The same selection and width filter apply to every simulated curve. Each
#' observed cluster's summed departure magnitude is compared with the largest qualifying
#' cluster in each simulation (zero if none), across both signs in a two-tailed test.
#' The inclusive rank gives a maximum-cluster adjusted p-value. Under the complete specified
#' null and exchangeable curves this controls the probability of reporting any cluster as
#' significant, accounting for the interval search. It is not strong error control for
#' arbitrary mixtures of true and false local nulls, and does not establish that every
#' angle inside a cluster is significant or that its endpoints are confidence limits.
#' Azimuth scoring counts each direction once and joins qualifying excursions across north for both observed and simulated curves. A start greater than the end denotes a wrapped interval. Declination retains its two endpoints.
#' The filter does not change the global test. The maximum-cluster principle is adapted
#' from Maris and Oostenveld (2007); cluster interpretation follows Sassenhagen and Draschkow
#' (2019).
#'
#' If all B+1 curves agree at an angle, its standard deviation is zero and it contributes zero.
#' When only the simulations agree, a differing observation remains in the pooled calculation
#' and can contribute to the score. Such points are also listed in
#' \code{result$data$constant.departures} with coordinate \code{x} (degrees), \code{observed}
#' and \code{simulated} densities and \code{direction} (above or below), for either tail setting.
#'
#' \strong{Monte Carlo ranks.} Global and local Monte Carlo p-values use \eqn{(b+1)/(B+1)}, where \eqn{B} is the number of
#' simulations and \eqn{b} counts simulated scores at least as large as the observed score,
#' including ties (North, Curtis and Sham, 2002). The minimum is \eqn{1/(B+1)}; all ties give 1.
#' The observation enters the pooled bounds, mean and standard deviation once. The envelope
#' describes each grid point separately. The rank test assumes that observed and simulated
#' curves are exchangeable under the specified null (Dufour, 2006). Apply the same
#' uncertainty and transformation procedures to both. If you estimate the null model
#' from these observations, account for that fitting when interpreting the test.
#'
#' \strong{Running simulations.} With a fixed seed, sampling in the calling session gives the same result for different worker
#' counts. For grids larger than 16 MiB, allow about eight bytes of temporary
#' disk space per simulation per grid point. Temporary files are removed on return or error.
#' @param pdf A \emph{skyscapeR.pdf} object created with either \code{\link{az.pdf}} or \code{\link{coordtrans}}
#' @param nsims Number of simulations, at least two. Default is 1000.
#' @param conf Confidence level as a fraction. For example, 0.95 requests a 95 percent envelope.
#' @param tails 1 for upper excesses, 2 for excesses and deficits (default).
#' @param normalise Logical to control whether to normalize SPDs. Default is FALSE
#' @param ncores Worker count. By default use half the detected available logical cores, rounded down with a minimum of one. Set ncores=1 for serial execution. The simulation count and memory batch size can further limit workers. Verbose output reports the effective count.
#' @param verbose Logical. Display progress or informational messages. Default \code{TRUE}.
#' @param .res Resolution with which to output probability distribution(s). Default is 0.1 degrees.
#' @param hor For azimuth input, optionally supply one shared \emph{skyscapeR.horizon}
#' to test in declination. For declination input, required: the original shared
#' horizon or an ordered list matching the distributions in \code{pdf}.
#' @param az Original azimuth \emph{skyscapeR.pdf}; required only when \code{pdf}
#' contains declination distributions. Must match its distribution count and
#' order. Do not supply this argument when \code{pdf} already contains azimuths.
#' @return A \emph{skyscapeR.sigTest} list containing metadata (test settings,
#' global and local p-values) and data (empirical SPD and null confidence envelope).
#' No original PDFs or horizon objects are attached.
#' @section Validation:
#' In S03, independent scoring reproduced the azimuth and declination ranks on their respective grids. Build 2.0.9006 counts north once and reports one interval for a qualifying excess across north. Refit earlier azimuth tests to obtain the revised probabilities. See \href{https://stonexsky.github.io/skyscapeR/manual/#s03-scoring-checks}{S03 validation in the technical manual}.
#' @export
#' @references Silva, F. (2020). A probabilistic framework and significance test for the analysis of structural orientations in skyscape archaeology. \emph{Journal of Archaeological Science}, 118, 105138. \doi{10.1016/j.jas.2020.105138}.
#'
#' Silva, F. and Lomsdalen, T. (2025). ‘No easy way from the earth to the stars’: a new statistical approach to the orientation of the Maltese temples. \emph{Archaeological and Anthropological Sciences}, 17, 96. \doi{10.1007/s12520-025-02208-4}.
#'
#' Maris, E. and Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG-data. \emph{Journal of Neuroscience Methods}. \doi{10.1016/j.jneumeth.2007.03.024}.
#'
#' Sassenhagen, J. and Draschkow, D. (2019). Cluster-based permutation tests of MEG/EEG data do not establish significance of effect latency or location. \emph{Psychophysiology}. \doi{10.1111/psyp.13335}.
#'
#' Dufour, J.-M. (2006). Monte Carlo tests with nuisance parameters: A general approach to finite-sample inference and nonstandard asymptotics. \emph{Journal of Econometrics}, 133, 443-477. \doi{10.1016/j.jeconom.2005.06.007}.
#'
#' North, B. V., Curtis, D. and Sham, P. C. (2002). A note on the calculation of empirical P values from Monte Carlo procedures. \emph{American Journal of Human Genetics}, 71, 439-441. \doi{10.1086/341527}.
#' @seealso \code{\link{az.pdf}}, \code{\link{coordtrans}}, \code{\link{spd}}, \code{\link{hpdi}}, \code{\link{modelTest}}
#' @examples
#' \dontrun{
#' # significance test for azimuth
#' Az <- az.pdf(az=c(87,93,90,110), unc=3)
#' st1 <- randomTest(Az, nsims=1000)
#' plot(st1)
#'
#' # significance test for declination
#' hor <- createHor(az=c(0,360), alt=c(0,0), loc=c(35,-8,25)) # flat horizon
#' Dec <- coordtrans(Az, hor)
#' st2 <- randomTest(Dec, az=Az, hor=hor, nsims=1000)
#' plot(st2)
#' }
#' @section Study design and numerical checks:
#' The simulated null describes independently sampled orientations with the
#' specified site geometry and measurement model. It does not automatically model
#' shared construction episodes, correlated measurement errors, geographical
#' clustering or selection of sites after inspecting their orientations. Such
#' structure needs a justified sampling model or separate sensitivity analysis.
#'
#' The pointwise envelope is not a simultaneous confidence band. A significant
#' global result rejects the specified random-orientation model; it does not
#' identify a celestial target. Refine the comparison grid and increase
#' \code{nsims} to assess borderline results, preserving the seed and inputs.
#' More simulations improve Monte Carlo resolution, not the validity of the
#' archaeological sampling design. Repeated searches over alternative subsets or
#' models are not included in the within-fit local-cluster adjustment.
randomTest <- compiler::cmpfun(randomTest_unc, options=list(enableJIT = 3))



#' @noRd
modelTest_unc <- function(pdf, model, unc=NULL, nsims=1000, conf=.95, tails=2, normalise=FALSE,
                          .res=.01, ncores=NULL, verbose=TRUE) {
  if (!inherits(pdf, "skyscapeR.pdf")) stop("Supply an empirical skyscapeR.pdf.")
  coord <- pdf$metadata$coord
  if (!identical(coord,"az") && !identical(coord,"dec"))
    stop("Supply an azimuth or declination skyscapeR.pdf.")
  n <- length(pdf$data)
  if (!n) stop("At least one observation is required.")
  if (!is.null(ncores) && (!is.numeric(ncores) || length(ncores)!=1L || !is.finite(ncores) ||
      ncores<1 || ncores!=floor(ncores))) stop("ncores must be a positive integer.")
  if (!is.numeric(nsims) || length(nsims)!=1L || !is.finite(nsims) ||
      nsims < 2 || nsims != floor(nsims)) stop("nsims must be an integer of at least two.")
  if (!is.numeric(conf) || length(conf)!=1L || !is.finite(conf) || conf<=0 || conf>=1 ||
      length(tails)!=1L || !tails %in% 1:2) stop("Use 0 < conf < 1 and tails=1 or 2.")
  if (!is.numeric(.res) || length(.res)!=1L || !is.finite(.res) || .res<=0)
    stop(".res must be one positive finite spacing.")
  if (verbose) .sky_message("Testing the specified model in ",if (coord=="az") "azimuth" else "declination"," (",nsims," simulations).")
  if (missing(model)) stop("Supply a model table or an ordered list of model tables.")
  if (is.data.frame(model)) model <- rep(list(model), n)
  if (!is.list(model) || length(model)!=n) stop("Supply one model table per observation.")
  xrange <- if (coord=="az") c(0,360) else c(-90,90)
  for (i in seq_along(model)) {
    tab <- model[[i]]
    if (!is.data.frame(tab) || ncol(tab)<2L || !nrow(tab))
      stop("Model ", i, " must be a nonempty table of coordinates and weights.")
    x <- tab[[1L]]
    w <- tab[[2L]]
    if (!is.numeric(x) || any(!is.finite(x) | x<xrange[1] | x>xrange[2]) ||
        !is.numeric(w) || any(!is.finite(w) | w<0) || !any(w>0))
      stop("Invalid coordinates or sampling weights in model ", i, ".")
  }
  template <- coord=="dec" && is.null(pdf$metadata$unc) && is.null(unc)
  settings <- list(coord=coord, xrange=xrange, normalise=normalise, resolution=.res)
  if (template) {
    family <- "shifted-template"
    settings$templates <- lapply(pdf$data,.declination_template)
    for (i in seq_len(n)) {
      tab <- model[[i]]
      bounds <- range(tab[[1L]][tab[[2L]]>0]) + range(settings$templates[[i]]$x)
      if (bounds[1L] < -90-1e-12 || bounds[2L] > 90+1e-12)
        stop("Model ",i," shifts part of its declination PDF outside [-90, 90]. Review its candidate means.")
    }
  } else {
    if (is.null(unc)) unc <- pdf$metadata$unc
    unc <- .measurement_unc(unc,n)
    family <- pdf$metadata$pdf
    if (is.null(family)) family <- "normal"
    settings$family <- family
    settings$unc <- unc
  }
  empirical <- spd(pdf, xrange=xrange, normalise=normalise, .res=.res)
  result <- .run_simulations(empirical$data, nsims, ncores, n, function(i) sampleList(model),
    if (template) .shifted_model_curve else .model_curve, settings, conf, tails,
    verbose)

  ## output
  out <- c()
  out$metadata$unc <- unc
  out$metadata$pdf <- family
  atmospheric <- intersect(c("refraction","atm","temp"), names(pdf$metadata$param))
  if (length(atmospheric)) out$metadata$source.atmosphere <- pdf$metadata$param[atmospheric]
  out$metadata$.res <- .res
  out$metadata$coord <- pdf$metadata$coord
  out$metadata$nsims <- nsims
  out$metadata$inference <- result$method
  out$metadata$local.inference <- result$local.method
  out$metadata$conf <- conf
  out$metadata$tails <- tails
  out$metadata$normalise <- normalise
  out$metadata$global.pval <- result$global
  out$metadata$local.pval <- result$local

  out$data$empirical <- empirical$data
  out$data$null.hyp <- result$envelope
  out$data$constant.departures <- result$departures
  class(out) <- 'skyscapeR.sigTest'

  return(out)
}

#' Test orientations against a specified null model
#' @description Use a Monte Carlo significance test to assess whether the measured orientations of a group of structures are compatible with a
#' specified null model in azimuth or declination. The model can be shared across sites or
#' specific to each site.
#' @details A proposed topographical or astronomical explanation can predict some orientations more
#' frequently than others. Testing that explanation requires simulating its predicted
#' distribution at the sites under study.
#'
#' Silva and Lomsdalen (2025) extended the random-orientation framework of Silva (2020)
#' to test specific hypotheses about temple orientations in Malta, including topographical and
#' astronomical explanations. Here the hypothesis specifies the relative probability of each
#' candidate orientation. A small global p-value is evidence against that particular model.
#'
#' The input PDF determines the coordinate of the test. Express the model coordinates and
#' uncertainties in that same coordinate. Each simulation draws one centre per observation
#' according to the weights in its model table. The function constructs an azimuth or declination
#' density around each centre and sums these densities on the same comparison grid as the
#' observed SPD. The observed measurement curves are used as supplied.
#'
#' \strong{Measurement uncertainty.} Direct measurement PDFs use their stored normal or uniform family and uncertainties.
#' For declinations from \code{\link{coordtrans}}, omit \code{unc} to use each observation's own density as a template.
#' A model draw represents the expected mean declination for that observation. The function calculates the
#' density-weighted mean of its piecewise-linear curve and adds (draw minus mean) to every declination.
#' For example, a curve with mean -17.6 degrees shifts by -4.4 degrees for a draw of -22 degrees.
#' The shift preserves the curve's width, asymmetry and peaks before interpolation onto the comparison grid.
#' It requires only the stored declination curve and assumes its shape applies at the modelled mean.
#' A positive-weight candidate that would shift density beyond -90 or +90 degrees raises an error.
#' An explicit \code{unc} selects a normal density around each simulated centre for transformed input.
#' With \code{normalise=TRUE}, each final SPD has unit area.
#'
#' \strong{Reading the result.} The global p-value compares the departure of the entire
#' observed curve with departures under the specified model. Local results use the largest
#' qualifying cluster per simulation to adjust for searching the grid. Boundaries describe
#' where departures were observed. The envelope is pointwise. Both tests use inclusive
#' Monte Carlo ranks, with minimum p-value 1/(nsims + 1).
#' See \code{\link{randomTest}} for the shared scoring rules, width threshold, treatment
#' of ties and zero variance, and assumptions behind the adjustment. When a model was
#' fitted to the same observations, its fitting procedure must be accounted for separately.
#'
#' \strong{Running simulations.} With a fixed seed, sampling in the calling session gives the same result for different worker
#' counts. For grids larger than 16 MiB, allow about eight bytes of temporary
#' disk space per simulation per grid point. Temporary files are removed on return or error.
#' @param pdf Empirical skyscapeR.pdf from \code{\link{az.pdf}}, \code{\link{dec.pdf}} or \code{\link{coordtrans}}. Its coordinate selects the simulation path.
#' @param model One shared data frame or an ordered list with one table per observation. The first column contains candidate coordinates in degrees and the second contains nonnegative discrete sampling weights, with at least one positive weight. On equally spaced grids, density heights can serve as weights.
#' @param nsims Number of simulations, at least two. Default is 1000.
#' @param conf Confidence level as a fraction.
#' @param tails 1 for upper excesses, 2 for excesses and deficits.
#' @param normalise Scale each summed density to unit area before comparison.
#' @param ncores Worker count. By default use half the detected available logical cores, rounded down with a minimum of one. Set ncores=1 for serial execution. The simulation count and memory batch size can further limit workers. Verbose output reports the effective count.
#' @param verbose Display simulation progress.
#' @param .res Spacing of the empirical and simulated comparison grids in degrees. Default 0.01.
#' @param unc Positive standard deviation (normal) or half-width (uniform), in degrees of the tested coordinate. Supply a scalar or one per observation. Omit to reuse stored measurement uncertainties or mean-shift transformed declination curves. See Details.
#' @return A skyscapeR.sigTest list with settings and global/local results in metadata, empirical density and null envelope in data.
#' @references Silva, F. (2020). A probabilistic framework and significance test for the analysis of structural orientations in skyscape archaeology. \emph{Journal of Archaeological Science}, 118, 105138. \doi{10.1016/j.jas.2020.105138}.
#'
#' Silva, F. and Lomsdalen, T. (2025). ‘No easy way from the earth to the stars’: a new statistical approach to the orientation of the Maltese temples. \emph{Archaeological and Anthropological Sciences}, 17, 96. \doi{10.1007/s12520-025-02208-4}.
#'
#' Maris, E. and Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG-data. \emph{Journal of Neuroscience Methods}. \doi{10.1016/j.jneumeth.2007.03.024}.
#'
#' Sassenhagen, J. and Draschkow, D. (2019). Cluster-based permutation tests of MEG/EEG data do not establish significance of effect latency or location. \emph{Psychophysiology}. \doi{10.1111/psyp.13335}.
#'
#' Dufour, J.-M. (2006). Monte Carlo tests with nuisance parameters: A general approach to finite-sample inference and nonstandard asymptotics. \emph{Journal of Econometrics}, 133, 443-477. \doi{10.1016/j.jeconom.2005.06.007}.
#'
#' North, B. V., Curtis, D. and Sham, P. C. (2002). A note on the calculation of empirical P values from Monte Carlo procedures. \emph{American Journal of Human Genetics}, 71, 439-441. \doi{10.1086/341527}.
#' @section Validation:
#' In S03, independent scoring reproduced the azimuth and declination ranks on their respective grids. Build 2.0.9006 counts north once and reports one interval for a qualifying excess across north. Refit earlier azimuth tests to obtain the revised probabilities. See \href{https://stonexsky.github.io/skyscapeR/manual/#s03-scoring-checks}{S03 validation in the technical manual}.
#' @export
#' @seealso \code{\link{az.pdf}}, \code{\link{dec.pdf}}, \code{\link{coordtrans}}, \code{\link{spd}}, \code{\link{randomTest}}, \code{\link{plot.skyscapeR.sigTest}}
#' @examples
#' d <- dec.pdf(dec = c(-18, 12), unc = c(0.8, 1.4))
#' m <- list(data.frame(x = -20:-16, w = 1), data.frame(x = 10:14, w = 1))
#' set.seed(42)
#' fit <- modelTest(d, m, nsims = 10, .res = 0.5,
#'                  verbose = FALSE)
#' print(fit)
#' \dontrun{
#' set.seed(42)
#' fit <- modelTest(d, m, nsims = 5000, ncores = 2)
#' # For transformed observations, model coordinates specify mean declinations.
#' Az <- az.pdf(az = c(87, 110), unc = 2)
#' hor <- createHor(az = c(0, 360), alt = c(0, 0), alt.unc = 0.1, loc = c(35, -8))
#' Dec <- coordtrans(Az, hor)
#' fit <- modelTest(Dec, data.frame(x = c(-22, 18), w = c(1, 1)))
#' }
#' @section Model weights and interpretation:
#' Weights are discrete sampling probabilities up to proportionality. If a
#' continuous model is evaluated on an uneven coordinate grid, density heights
#' alone do not represent probability mass; construct weights that account for
#' the bin widths before calling the function. Put site-specific model tables in
#' the same order as the observed measurements.
#'
#' A large p-value does not prove the proposed explanation or measure its
#' probability. The test asks whether the specified model is incompatible with
#' the observed departure statistic. Alternative models can remain compatible.
#' The simulations do not automatically represent dependence between sites,
#' selection after viewing the data, or uncertainty in fitted model parameters.
#' See \code{\link{randomTest}} for study-design and Monte Carlo limitations.
modelTest <- compiler::cmpfun(modelTest_unc, options=list(enableJIT = 3))


# HWT stores lower bin edges. Supply the neighbouring bins across true north.
.coord_horizon <- function(horizon) {
  d <- horizon$data
  if (is.null(horizon$metadata$ID)) return(d)
  d <- d[order(d$az), , drop=FALSE]
  first <- d[1L, , drop=FALSE]
  last <- d[nrow(d), , drop=FALSE]
  first$az <- first$az + 360
  last$az <- last$az - 360
  rbind(last,d,first)
}


.unwrap_north_pdf <- function(curve) {
  gap <- which(diff(curve$x)>180)
  if (length(gap)!=1L) return(curve)
  west <- seq.int(gap+1L,nrow(curve))
  curve$x[west] <- curve$x[west]-360
  curve[c(west,seq_len(gap)),,drop=FALSE]
}

# Preserve native coordinates where covered. Try an equivalent bearing only
# outside that interval, leaving actual undercoverage to the caller's checks.
.profile_bearings <- function(x, bearings) {
  bounds <- range(bearings)
  outside <- x < bounds[1] | x > bounds[2]
  x[outside] <- x[outside] %% 360
  x
}
