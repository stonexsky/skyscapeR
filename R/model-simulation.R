.simulation_cores <- function() {
  limits <- parallel::detectCores(logical=TRUE)
  if (Sys.info()[["sysname"]] == "Linux") {
    affinity <- tryCatch(parallel::mcaffinity(), error=function(e) NULL)
    if (length(affinity)) limits <- c(limits,length(affinity))
    if (file.exists("/sys/fs/cgroup/cpu.max")) {
      quota <- scan("/sys/fs/cgroup/cpu.max",what="",quiet=TRUE,nmax=2)
      if (length(quota)==2L && quota[1L]!="max")
        limits <- c(limits,max(1,floor(as.numeric(quota[1L])/as.numeric(quota[2L]))))
    }
  }
  check <- Sys.getenv("_R_CHECK_LIMIT_CORES_")
  limits <- c(limits, if (check %in% c("TRUE","true","warn","warned")) 2,
              suppressWarnings(as.numeric(Sys.getenv(c("SLURM_CPUS_PER_TASK","NSLOTS")))))
  limits <- limits[is.finite(limits) & limits>=1]
  if (length(limits)) min(limits) else 1L
}

# Use the same default and report the effective worker count in all parallel entry points.
.worker_count <- function(ncores,limit,verbose) {
  available <- .simulation_cores()
  if (is.null(ncores)) ncores <- max(1L,floor(available/2))
  workers <- min(ncores,available,limit)
  if (verbose) .sky_message(sprintf("Using %d of %d available logical cores.",workers,available))
  workers
}

# Workers evaluate each curve with settings installed once at cluster startup.
.worker_curve <- function(draw) {
  get("curve", envir=.GlobalEnv, inherits=FALSE)(draw,
    get("settings", envir=.GlobalEnv, inherits=FALSE))
}

.declination_template <- function(curve) {
  x <- curve$x
  y <- curve$y
  if (!is.numeric(x) || !is.numeric(y) || length(x)!=length(y) || length(x)<2L ||
      any(!is.finite(x)) || any(!is.finite(y)) || any(diff(x)<=0) || any(y<0) || !any(y>0))
    stop("A declination template needs increasing coordinates and a finite, nonnegative density with positive area.")
  # Retain the neighbouring zero at each edge for linear interpolation.
  positive <- which(y>0)
  keep <- seq.int(max(1L,positive[1L]-1L), min(length(x),tail(positive,1L)+1L))
  x <- x[keep]
  y <- y[keep]
  dx <- diff(x)
  left <- head(y,-1L)
  right <- tail(y,-1L)
  mass <- dx*(left+right)/2
  # Exact first moment of the same piecewise-linear curve used by spd().
  centre <- sum(head(x,-1L)*mass + dx^2*(left+2*right)/6)/sum(mass)
  data.frame(x=x-centre,y=y)
}

.shifted_model_curve <- function(draw, settings) {
  curves <- Map(function(template, centre) {
    template$x <- template$x + centre
    template
  }, settings$templates, draw)
  sim <- structure(list(metadata=list(coord="dec"),data=curves),class="skyscapeR.pdf")
  spd(sim,xrange=settings$xrange,normalise=settings$normalise,.res=settings$resolution)$data$y
}

.model_curve <- function(draw, settings) {
  density <- if (settings$coord == "dec") dec.pdf else az.pdf
  sim <- density(draw, settings$unc, pdf=settings$family, verbose=FALSE)
  spd(sim, xrange=settings$xrange, normalise=settings$normalise, .res=settings$resolution)$data$y
}

.random_curve <- function(draw, settings) {
  sim <- az.pdf(pdf=settings$family, az=draw, unc=settings$unc, verbose=FALSE)
  if (settings$quick) {
    transform <- settings$transform
    transformed <- coordtrans_unc(spd(sim), transform$hor, xrange=settings$xrange,
      refraction=transform$refraction, atm=transform$atm, temp=transform$temp, .res=transform$.res, verbose=FALSE)$data
    if (settings$normalise) return(.normalise_spd_area(transformed$x, transformed$y))
    return(transformed$y)
  }
  if (settings$coord == "dec") {
    transform <- settings$transform
    sim <- coordtrans_unc(sim, transform$hor, refraction=transform$refraction, atm=transform$atm, temp=transform$temp,
      .res=transform$.res, .prec=transform$.prec, .output.res=transform$.output.res, verbose=FALSE)
  }
  spd(sim, xrange=settings$xrange, normalise=settings$normalise, .res=settings$resolution)$data$y
}

.run_simulations <- function(observed, nsims, ncores, n, draw, curve, settings,
                             conf, tails, verbose) {
  width <- nrow(observed)
  disk <- 8 * as.double(nsims) * width > 16*1024^2
  batch.size <- max(1, floor(8*1024^2 / (8 * (as.double(width)+n))))
  workers <- .worker_count(ncores,min(nsims,batch.size),verbose)
  index <- seq_len(nsims)
  batches <- unname(split(index, ceiling(index/batch.size)))
  paths <- character()
  memory <- NULL
  if (disk) {
    work <- tempfile("skyscapeR-model-")
    if (!dir.create(work)) stop("Cannot create simulation temporary directory.")
    on.exit(unlink(work,recursive=TRUE),add=TRUE)
    paths <- file.path(work,as.character(seq_along(batches)))
  } else memory <- matrix(0, nsims, width)
  progress <- .sky_progress(verbose)
  cl <- NULL
  if (workers>1) {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(if (!is.null(cl)) parallel::stopCluster(cl),add=TRUE)
    invisible(parallel::clusterCall(cl, .libPaths, .libPaths()))
    if (!all(unlist(parallel::clusterCall(cl, requireNamespace, "skyscapeR", quietly=TRUE))))
      stop("Socket workers could not load skyscapeR.")
    parallel::clusterExport(cl, c("curve", "settings"), envir=environment())
  }
  if (progress) {
    .sky_message("Running ", nsims, " simulations.")
    pb <- txtProgressBar(max=nsims,style=3)
    on.exit(if (!is.null(pb)) close(pb),add=TRUE)
  }
  # Keep storage batches independent of the frequency of progress updates.
  wave <- if (progress) workers * max(1, ceiling(nsims/(100*workers))) else batch.size
  for (b in seq_along(batches)) {
    rows <- batches[[b]]
    block <- matrix(0,length(rows),width)
    first <- 1L
    while (first <= length(rows)) {
      count <- if (progress && b==1L && first==1L) workers else wave
      last <- min(length(rows), first+count-1L)
      part <- first:last
      draws <- lapply(rows[part],draw)
      curves <- if (is.null(cl)) lapply(draws,curve,settings) else parallel::parLapply(cl,draws,.worker_curve)
      block[part,] <- do.call(rbind,curves)
      rm(draws,curves)
      if (progress) setTxtProgressBar(pb,rows[last])
      first <- last+1L
    }
    if (disk) writeBin(as.double(block),paths[b]) else memory[rows,] <- block
    rm(block)
  }
  if (progress) {
    close(pb)
    pb <- NULL
  }
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  result <- .simulation_inference(observed, nsims, batches, paths, memory, conf, tails, progress,
                                  circular=identical(settings$coord,"az"))
  result
}

# Each file contains a row batch in column-major order for contiguous tile reads.
.model_tile <- function(path, nrows, first, count) {
  con <- file(path, "rb")
  on.exit(close(con))
  seek(con, where=8 * nrows * (first-1), origin="start")
  size <- nrows * as.double(count)
  vals <- readBin(con, double(), n=size)
  if (length(vals) != size) stop("Incomplete simulation temporary file.")
  matrix(vals, nrow=nrows)
}

.simulation_inference <- function(observed, nsims, batches, paths, memory, conf, tails, progress=FALSE, circular=FALSE) {
  width <- nrow(observed)
  # Preserve the closing sample for plotting, but rank only distinct directions.
  duplicate.north <- circular && width>1L && abs(diff(range(observed$x))-360)<1e-10
  last.index <- width-as.integer(duplicate.north)
  disk <- length(paths) > 0L
  stride <- max(1, floor(8*1024^2 / (8 * (as.double(nsims)+1))))
  mean.sim <- sd.sim <- upper <- null.mean <- numeric(width)
  null.constant <- logical(width)
  lower <- if (tails == 2) numeric(width)
  probs <- if (tails == 2) c(1-(1-conf)/2, (1-conf)/2) else conf
  if (progress) {
    .sky_message("Computing the simulation envelope.")
    pb <- txtProgressBar(max=width,style=3)
    on.exit(if (!is.null(pb)) close(pb),add=TRUE)
  }
  col.tick <- max(1L,ceiling(width/100))
  for (first in seq(1, width, by=stride)) {
    cols <- seq(first, min(width, first+stride-1))
    tile <- if (disk) matrix(0, nsims, length(cols)) else memory[,cols,drop=FALSE]
    if (disk) for (b in seq_along(batches)) {
      rows <- batches[[b]]
      tile[rows,] <- .model_tile(paths[b],length(rows),first,length(cols))
    }
    # Disk blocks store doubles. Use the same arithmetic for integer inputs.
    if (!is.double(tile)) storage.mode(tile) <- "double"
    null.mean[cols] <- colMeans(tile)
    for (j in seq_along(cols)) {
      k <- cols[j]
      if (progress && k>1L && (k-1L) %% col.tick==0) setTxtProgressBar(pb,k-1L)
      null.constant[k] <- all(tile[,j] == tile[1L,j])
      # Sorting makes the numerical summaries invariant to row relabelling.
      # Pool only one column: no extra (B+1)-by-grid matrix is retained.
      v <- sort(c(observed$y[k],tile[,j]))
      mean.sim[k] <- mean(v)
      if (all(v == v[1L])) {
        mean.sim[k] <- v[1L]
        upper[k] <- v[1L]
        if (tails==2) lower[k] <- v[1L]
        next
      }
      sd.sim[k] <- stats::sd(v)
      q <- stats::quantile(v, probs, names=FALSE)
      upper[k] <- q[1L]
      if (tails==2) lower[k] <- q[2L]
    }
  }
  if (progress) {
    setTxtProgressBar(pb,width)
    close(pb)
    pb <- NULL
    .sky_message("Scoring the simulated departures.")
    pb <- txtProgressBar(max=nsims,style=3)
  }
  row.tick <- max(1L,ceiling(nsims/100))
  constant <- sd.sim == 0
  flagged <- which(null.constant & observed$y != null.mean)
  flagged <- flagged[flagged<=last.index]
  departures <- data.frame(x=observed$x[flagged], observed=observed$y[flagged],
    simulated=null.mean[flagged],
    direction=c("below","above")[1L + (observed$y[flagged] > null.mean[flagged])])
  # Score the observation and each simulation against the same density bounds.
  departure <- function(y) {
    score <- numeric(width)
    above <- which(!constant & y > upper)
    score[above] <- (y[above]-upper[above])/sd.sim[above]
    if (tails==2) {
      below <- which(!constant & y < lower)
      score[below] <- (y[below]-lower[below])/sd.sim[below]
    }
    if (duplicate.north) score[width] <- 0
    score
  }
  observed.score <- departure(observed$y)
  above <- which(observed.score > 0)
  below <- which(observed.score < 0)
  runs <- function(i) {
    if (!length(i)) return(list())
    groups <- split(i,cumsum(c(TRUE,diff(i)>1)))
    if (circular && length(groups)>1L && i[1L]==1L && tail(i,1L)==last.index) {
      groups <- c(list(c(groups[[length(groups)]],groups[[1L]])),groups[-c(1L,length(groups))])
    }
    groups[vapply(groups,function(g) {
      span <- observed$x[tail(g,1L)]-observed$x[g[1L]]
      if (span<0 && circular) span <- span+360
      length(g)>1L && span>=0.5-1000*.Machine$double.eps
    },logical(1))]
  }
  pos <- runs(above)
  intervals <- c(pos,runs(below))
  types <- rep(c("+","-"),c(length(pos),length(intervals)-length(pos)))
  global <- sum(observed.score[above]) - sum(observed.score[below])
  threshold <- vapply(intervals, function(i) abs(sum(observed.score[i])), numeric(1))
  extreme <- 0
  local.extreme <- numeric(length(intervals))
  for (b in seq_along(batches)) {
    rows <- batches[[b]]
    block <- if (disk) .model_tile(paths[b],length(rows),1,width) else memory[rows,,drop=FALSE]
    for (r in seq_len(nrow(block))) {
      score <- departure(block[r,])
      up <- pmax(score,0)
      down <- pmax(-score,0)
      extreme <- extreme + (sum(up)+sum(down) >= global)
      # Apply the same interval selection to every curve, then compare each
      # observed cluster with the simulated maximum over both directions.
      simulated.intervals <- c(runs(which(score>0)),runs(which(score<0)))
      maximum <- max(c(0,vapply(simulated.intervals,
                                function(i) abs(sum(score[i])),numeric(1))))
      local.extreme <- local.extreme + (maximum >= threshold)
      if (progress && (rows[r] %% row.tick==0 || rows[r]==nsims)) setTxtProgressBar(pb,rows[r])
    }
  }
  pval <- function(b) (b+1)/(nsims+1)
  local <- data.frame(type=if (length(types)) types else logical(),
    start=vapply(intervals,function(i) observed$x[i[1L]],numeric(1)),
    end=vapply(intervals,function(i) observed$x[tail(i,1L)],numeric(1)),
    p.value=pval(local.extreme),row.names=NULL)
  rownames(local) <- NULL
  envelope <- list(x=observed$x, CE.mean=mean.sim, CE.upper=upper)
  if (tails==2) envelope$CE.lower <- lower
  list(global=pval(extreme),local=local,envelope=envelope,departures=departures,
       method="pooled-exchangeable-v1",local.method="maximum-cluster-v1")
}
