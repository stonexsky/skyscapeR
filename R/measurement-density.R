# Shared by the public constructors and by simulations which move the same
# measurement model to a different bearing or declination.
.measurement_family <- function(pdf) {
  families <- c("normal", "uniform", "epanechnikov", "tapered_uniform",
                "triangular", "uniform_normal")
  if (is.character(pdf) && length(pdf) == 1L && !is.na(pdf)) {
    family <- switch(tolower(pdf), n = "normal", u = "uniform", tolower(pdf))
    if (family %in% families) return(family)
  }
  stop("Choose a measurement density: ", paste(families, collapse = ", "), ".")
}

.measurement_values <- function(x, n, label) {
  if (!is.numeric(x) || !is.null(dim(x)) || !length(x) %in% c(1L, n) || any(!is.finite(x)))
    stop(label, " must be finite, with one value or one per observation.")
  rep_len(x, n)
}

.measurement_spec <- function(centre, coord, unc, family, taper, survey.sd, lower, upper) {
  n <- length(centre)
  shape <- list()
  if (!is.null(taper) || family == "tapered_uniform") {
    if (family != "tapered_uniform") stop("taper applies only to tapered_uniform.")
    shape$taper <- .measurement_values(taper, n, "taper")
    if (any(shape$taper < 0)) stop("taper must be nonnegative, in degrees.")
  }
  if (!is.null(survey.sd) || family == "uniform_normal") {
    if (family != "uniform_normal") stop("survey.sd applies only to uniform_normal.")
    shape$survey.sd <- .measurement_values(survey.sd, n, "survey.sd")
    if (any(shape$survey.sd < 0)) stop("survey.sd must be nonnegative, in degrees.")
  }
  if (!is.null(lower) || !is.null(upper)) {
    if (family != "triangular" || !is.null(unc))
      stop("Use lower and upper together for triangular, omitting unc.")
    left <- centre - .measurement_values(lower, n, "lower")
    right <- .measurement_values(upper, n, "upper") - centre
    if (coord == "az") {
      left <- left %% 360
      right <- right %% 360
      if (any(left + right >= 360))
        stop("The clockwise arc from lower to upper must contain the modal azimuth and span less than one turn.")
    }
    if (any(left < 0 | right < 0 | left + right <= 0 | !is.finite(left + right)))
      stop("The triangular limits must enclose the mode and have positive separation.")
    unc <- (left + right) / 2
    shape$left.fraction <- left / (left + right)
  }
  list(unc = .measurement_unc(unc, n), shape = shape)
}

# Gaussian smoothing of an interval. Reflection avoids subtracting two CDFs
# close to one. For very small apertures use the interval-average Taylor limit.
.uniform_normal_density <- function(delta, width, sigma) {
  z <- -abs(delta) / sigma
  h <- width / sigma
  if (h < 1e-4) {
    y <- stats::dnorm(delta, sd = sigma)
    keep <- y > 0
    y[keep] <- y[keep] * (1 + (z[keep]^2 - 1) * h^2 / 6)
    return(y)
  }
  (stats::pnorm((width - abs(delta)) / sigma) -
     stats::pnorm((-width - abs(delta)) / sigma)) / (2 * width)
}

.line_density <- function(delta, width, family, shape) {
  switch(family,
    epanechnikov = (0.75 / width) * pmax(0, 1 - (delta / width)^2),
    tapered_uniform = {
      distance <- pmax(0, abs(delta) - width)
      fraction <- pmin(1, distance / shape$taper)
      (1 + cos(pi * fraction)) / (2 * (2 * width + shape$taper))
    },
    triangular = {
      left <- right <- width
      if (!is.null(shape$left.fraction)) {
        left <- (2 * width) * shape$left.fraction
        right <- (2 * width) * (1 - shape$left.fraction)
      }
      height <- as.numeric(delta == 0)
      if (left > 0) height[delta < 0] <- pmax(0, 1 + delta[delta < 0] / left)
      if (right > 0) height[delta > 0] <- pmax(0, 1 - delta[delta > 0] / right)
      height / width
    },
    uniform_normal = .uniform_normal_density(delta, width, shape$survey.sd)
  )
}

.measurement_density <- function(x, centre, width, family, shape, coord) {
  if ((family == "tapered_uniform" && shape$taper == 0) ||
      (family == "uniform_normal" && shape$survey.sd == 0)) family <- "uniform"
  if (coord == "dec") {
    if (family == "normal") return(stats::dnorm(x, centre, width))
    if (family == "uniform") return(stats::dunif(x, centre - width, centre + width))
    return(.line_density(x - centre, width, family, shape))
  }
  delta <- (x - (centre %% 360) + 180) %% 360 - 180
  if (family == "uniform") {
    turns <- floor((width - delta) / 360) + floor((width + delta) / 360) + 1
    return((turns / width) / 2)
  }
  if (family %in% c("normal", "uniform_normal")) {
    sigma <- if (family == "normal") width else shape$survey.sd
    if (sigma >= 360) {
      # At this scale only the first wrapped Gaussian harmonic is measurable.
      # A uniform aperture multiplies that harmonic by sinc(half-width).
      factor <- 1
      if (family == "uniform_normal") {
        angle <- width * (pi / 180)
        factor <- if (angle == 0) 1 else sin(angle) / angle
      }
      return((1 + 2 * factor * exp(-0.5 * (sigma * (pi / 180))^2) * cos(delta * pi / 180)) / 360)
    }
    extent <- 9 * sigma
    kernel <- function(at) stats::dnorm(at, sd = sigma)
    if (family == "uniform_normal") {
      extent <- extent + width
      kernel <- function(at) .uniform_normal_density(at, width, sigma)
    }
    y <- kernel(delta)
    for (turn in seq_len(ceiling((extent + 180) / 360)))
      y <- y + kernel(delta + 360 * turn) + kernel(delta - 360 * turn)
    return(y)
  }
  y <- .line_density(delta, width, family, shape)
  # A skew triangle can extend past the antipode on one side of its mode.
  if (family == "triangular" && !is.null(shape$left.fraction) && width > 90)
    y <- y + .line_density(delta + 360, width, family, shape) +
             .line_density(delta - 360, width, family, shape)
  y
}

.measurement_pdf <- function(centre, coord, unc, pdf, shape, name, .res, .cutoff, verbose) {
  family <- .measurement_family(pdf)
  circular <- coord == "az"
  if (!is.numeric(.res) || length(.res) != 1L || !is.finite(.res) || .res <= 0 ||
      (circular && .res >= 360) || !is.numeric(.cutoff) || length(.cutoff) != 1L ||
      !is.finite(.cutoff) || .cutoff < 0)
    stop("Supply a positive finite grid spacing (below 360 for bearings) and a nonnegative finite cutoff.")
  if (length(name) != length(centre)) stop("Supply one name per measurement.")
  if (circular && !family %in% c("normal", "uniform")) {
    outer <- unc
    if (family == "tapered_uniform") outer <- outer + shape$taper
    if (any(outer > 180)) stop("The azimuth interval, including any cosine tapers, must span at most one turn.")
  }
  x <- if (circular) seq(0, 360, by = .res) else seq(-90, 90, by = .res)
  if (circular) x <- x[x < 360]
  progress <- circular && .sky_progress(verbose) && length(centre) > 1L
  if (progress) {
    pb <- utils::txtProgressBar(max = length(centre), style = 3)
    on.exit(close(pb))
  }
  data <- lapply(seq_along(centre), function(i) {
    params <- lapply(shape, `[[`, i)
    y <- .measurement_density(x, centre[i], unc[i], family, params, coord)
    keep <- y >= .cutoff
    if ((!circular || !family %in% c("normal", "uniform")) &&
        (sum(keep) < 2L || !any(y[keep] > 0)))
      stop("Too few density samples for measurement ", i, ". Reduce .res or .cutoff.")
    if (progress) utils::setTxtProgressBar(pb, i)
    data.frame(x = x[keep], y = y[keep])
  })
  metadata <- list(coord = coord, pdf = pdf, name = name)
  metadata[[coord]] <- centre
  metadata$unc <- unc
  if (length(shape)) metadata$shape <- shape
  metadata$param <- list(.cutoff = .cutoff, .res = .res)
  if (!circular && verbose) .sky_message("Created ", length(data), " declination densities.")
  structure(list(metadata = metadata, data = data), class = "skyscapeR.pdf")
}
