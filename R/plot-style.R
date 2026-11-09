# Shared drawing conventions. Angular coordinates always belong to the caller.
.plot_style <- function(style) {
  mono <- match.arg(style, c("publication-colour", "publication-mono")) == "publication-mono"
  list(mono=mono, ink=if (mono) "#282828" else "#23383A",
       muted=if (mono) "#5F5F5F" else "#51676A", grid=if (mono) "#E2E2E2" else "#DEE5E3",
       data=if (mono) "#303030" else "#12696E", target=if (mono) "#686868" else "#A9542E",
       lower=if (mono) "#555555" else "#725786", envelope=if (mono) "#D1D1D1" else "#C7D0CE")
}

.plot_colour <- function(col, theme) {
  if (!theme$mono || !length(col)) return(col)
  rgb <- grDevices::col2rgb(col, alpha=TRUE)/255
  grey <- drop(c(.2126,.7152,.0722) %*% rgb[1:3,,drop=FALSE])
  grDevices::gray(grey, alpha=rgb[4,])
}

.polar_xy <- function(angle, radius=1) {
  a <- angle*pi/180
  cbind(x=radius*sin(a), y=radius*cos(a))
}

.site_caption <- function(loc, name=NULL) {
  if (length(name) && !is.na(name[1]) && nzchar(name[1])) return(as.character(name[1]))
  if (!is.numeric(loc) || length(loc)<2L || any(!is.finite(loc[1:2]))) return("")
  sprintf("%.2f\u00b0%s  %.2f\u00b0%s",abs(loc[1]),if(loc[1]<0) "S" else "N",
          abs(loc[2]),if(loc[2]<0) "W" else "E")
}

.target_id <- function(group,i) {
  ids <- colnames(group$decs)
  if (is.null(ids) || is.na(ids[i]) || !nzchar(ids[i])) return(paste("Target",i))
  ids[i]
}

.target_name <- function(id) {
  ids <- c(jS="June solstice", dS="December solstice", eq="Equinox",
           nMjLX="Northern major lunar extreme", sMjLX="Southern major lunar extreme",
           nmnLX="Northern minor lunar extreme", smnLX="Southern minor lunar extreme")
  ifelse(id %in% names(ids), unname(ids[id]), id)
}

# Allocate one caller-owned plot slot. Each panel uses identical relative bounds
# regardless of its labels, so linked marginal axes coincide physically.
.plot_panels <- function() {
  par(mar=rep(0,4))
  plot.new()
  fin <- par("fin")
  function(left, bottom, right, top) {
    par(mai=c(bottom*fin[2],left*fin[1],(1-top)*fin[2],(1-right)*fin[1]), new=TRUE)
  }
}

.plot_ticks <- function(limits, n=5) {
  ticks <- pretty(limits,n=n)
  ticks[ticks >= min(limits) & ticks <= max(limits)]
}

.plot_x_baseline <- function(pos=NULL, col=par("fg"), lwd=.65) {
  axis(1, at=par("usr")[1:2], labels=FALSE, pos=pos, tcl=0, lwd=lwd, col=col)
}

.plot_axis <- function(side, at, labels=TRUE, theme, pos=NULL, ...) {
  if (side == 1) .plot_x_baseline(pos, theme$muted)
  axis(side, at=at, labels=labels, las=1, col=theme$muted, col.axis=theme$muted,
       pos=pos, lwd=if (side == 1) 0 else .65, lwd.ticks=.65, tcl=-.19, ...)
}

.angle_title <- function(coord) {
  paste0(if (coord=="az") "True azimuth" else "Declination", " (\u00b0)")
}

# Separate nearby external rows without moving the data-bearing endpoints.
.label_rows <- function(y, lower, upper, gap) {
  if (!length(y)) return(y)
  ord <- order(y)
  gap <- min(gap,(upper-lower)/max(1,length(y)-1))
  v <- pmin(upper,pmax(lower,y[ord]))
  if (length(v)>1) {
    for (i in 2:length(v)) v[i] <- max(v[i],v[i-1]+gap)
    v[length(v)] <- min(upper,v[length(v)])
    for (i in rev(seq_len(length(v)-1L))) v[i] <- min(v[i],v[i+1]-gap)
  }
  y[ord] <- v
  y
}
