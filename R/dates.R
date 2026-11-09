.calendar_flag <- function(calendar) {
  if (!is.character(calendar) || length(calendar)!=1L || is.na(calendar) ||
      !calendar %in% c('G','Gregorian','J','Julian'))
    stop('calendar must be G, Gregorian, J or Julian.', call.=FALSE)
  if (calendar %in% c('G','Gregorian')) 1L else 0L
}

.timezone_spec <- function(timezone) {
  if (is.numeric(timezone) && length(timezone)==1L && is.finite(timezone) &&
      abs(timezone)<=14) return(list(offset=timezone*3600))
  if (!is.character(timezone) || length(timezone)!=1L || is.na(timezone))
    stop('Use a timezone name or a fixed offset in hours between -14 and 14.', call.=FALSE)
  if (timezone %in% c('', 'UTC','GMT')) return(list(offset=0))
  m <- regmatches(timezone, regexec('^(UTC|GMT)([+-])([0-9]{1,2})(:([0-9]{2}))?$',timezone))[[1]]
  if (length(m)) {
    h <- as.numeric(m[4]); minutes <- if (nzchar(m[6])) as.numeric(m[6]) else 0
    if (minutes>=60 || h+minutes/60>14) stop('Invalid fixed UTC offset.',call.=FALSE)
    return(list(offset=(if(m[3]=='-') -1 else 1)*(h*3600+minutes*60)))
  }
  if (!timezone %in% OlsonNames()) stop('Unknown timezone. Use an Olson name or UTC offset.',call.=FALSE)
  list(zone=timezone)
}

.named_offset <- function(jd, zone) {
  yr <- swephR::swe_revjul(jd,1L)$year
  if (yr<1900 || yr>9999) stop('Named timezones require Gregorian years 1900-9999. Use UTC or an explicit fixed offset for ancient dates.',call.=FALSE)
  stamp <- as.POSIXct((jd-2440587.5)*86400,origin='1970-01-01',tz='UTC')
  local <- as.POSIXlt(stamp,tz=zone)
  if (is.na(local$year)) stop('Date is not supported by the system timezone database.',call.=FALSE)
  local_jd <- swephR::swe_julday(local$year+1900,local$mon+1,local$mday,
                               local$hour+local$min/60+local$sec/3600,1L)
  # IANA timezone offsets are integral seconds, including historical offsets.
  round((local_jd-jd)*86400)
}

.local_to_ut <- function(local_jd, spec) {
  if (!is.null(spec$offset)) return(local_jd-spec$offset/86400)
  # Collect offsets either side of a transition, then validate each candidate.
  # Use a clamped range near the supported year boundaries.
  lo <- swephR::swe_julday(1900,1,1,0,1L)
  hi <- swephR::swe_julday(10000,1,1,0,1L)-1/86400
  if (local_jd<lo || local_jd>hi) stop('Named timezones require Gregorian years 1900-9999. Use UTC or a fixed offset.',call.=FALSE)
  probes <- pmin(hi,pmax(lo,local_jd+seq(-2,2,by=.25)))
  offsets <- unique(vapply(probes,.named_offset,numeric(1),zone=spec$zone))
  candidates <- local_jd-offsets/86400
  valid <- vapply(seq_along(candidates),function(i) {
    candidates[i]>=lo && candidates[i]<=hi &&
      .named_offset(candidates[i],spec$zone)==offsets[i]
  },logical(1))
  candidates <- candidates[valid]
  if (!length(candidates)) stop('Nonexistent local time in this timezone (clock change). Supply a valid time or a fixed offset.',call.=FALSE)
  if (length(candidates)>1L) stop('Ambiguous local time in this timezone (clock change). Supply an explicit fixed offset.',call.=FALSE)
  candidates[1]
}

#' Convert local calendar date and time to a UT Julian date
#' @description Place a dated field observation or historical event on the continuous UT Julian-date
#' scale used in ephemeris calculations.
#' @details Civil dates require both a calendar and a time convention to specify an instant. The
#' function parses each supplied date and time, verifies that the calendar date exists by
#' converting it to a Julian date and back, then applies the chosen offset or civil
#' timezone to obtain UT. A date without a time uses local midnight. Astronomical year
#' numbering includes year zero. An era label can also precede the month, for example
#' \code{"3700 BCE/06/21"}.
#'
#' Fixed offsets may be numeric hours east of UTC (for example 5.5), or
#' strings such as \code{"UTC+05:30"}, \code{"UTC-02"} or \code{"GMT+01"}.
#' Offsets must lie between -14 and +14 hours. UTC, GMT and the empty string
#' mean zero offset; an empty string does not select the system timezone.
#'
#' Named zones listed by \code{OlsonNames()}, such as \code{"Europe/London"},
#' use the installed timezone database for the actual Gregorian date, not a
#' surrogate year. They are supported only for Gregorian years 1900-9999.
#' Earlier dates require UTC or a fixed offset. For ancient-date workflows,
#' set \code{skyscapeR.vars(timezone="UTC")} or provide timezone explicitly.
#' The package's initial timezone remains Europe/London. Database versions
#' may differ across systems; fixed offsets provide reproducible conversions.
#'
#' Gregorian and Julian calendars are proleptic: no automatic calendar reform
#' is applied. Astronomical year numbering is used: 0 is 1 BCE, -1 is 2 BCE.
#' Fixed-offset input dates support years from -100000 to 100000.
#' A timezone offset represents civil clock time, not local solar time.
#' These functions apply clock offsets to a UT Julian day; they do not apply
#' leap-second, UTC-to-UT1 or TT corrections. Leap-second input is rejected.
#'
#' Ambiguous and nonexistent local times at clock transitions are rejected.
#' For an ambiguous time, specify a fixed offset to identify the intended instant.
#' Calendar dates and clock components are validated before conversion.
#' @param time Character vector in YYYY/MM/DD format, optionally followed by HH:MM or HH:MM:SS. Fractional seconds are accepted. Omitted clock time means midnight.
#' @param timezone Timezone name or fixed offset; defaults to the value in \code{skyscapeR.vars()}. See Details.
#' @param calendar Calendar: G or Gregorian, J or Julian. Defaults to the package calendar setting, initially Gregorian.
#' @param verbose Logical; display a conversion message. Default is FALSE.
#' @return A numeric vector of UT Julian dates, in the same order as time.
#' @import swephR
#' @section Validation:
#' Build 2.0.9005 passed all 499 calendar, clock, era, angle and setting checks in A01. The largest discrepancy across 240 absolute Julian dates was 0.000021 seconds. See \href{https://stonexsky.github.io/skyscapeR/manual/#a01-calendar-checks}{A01 validation in the technical manual}.
#' @export
#' @seealso \code{\link{jd2time}}, \code{\link{timestring}}, \code{\link{BC.AD}}, \code{\link{skyscapeR.vars}}
#' @examples
#' time2jd("2000/01/01 12:00:00", "UTC") # 2451545
#' # A positive offset crosses into the preceding UTC year:
#' jd <- time2jd("2020/01/01 00:30:00", "UTC+02")
#' jd2time(jd, "UTC") # 2019/12/31 22:30:00
#' # Ancient Julian-calendar date, with an explicit fixed offset:
#' time2jd("-4000/06/21 06:00:00", "UTC+02", calendar="J")
#' # A modern date uses that year's timezone rules:
#' time2jd("2007/03/15 12:00:00", "America/New_York")
time2jd <- function(time, timezone, calendar, verbose=F) {
  if (missing(timezone)) timezone <- skyscapeR.env$timezone
  if (missing(calendar)) calendar <- skyscapeR.env$calendar
  cal <- .calendar_flag(calendar); spec <- .timezone_spec(timezone)
  if (!is.character(time) || anyNA(time)) stop('time must be a character vector without missing values.',call.=FALSE)
  time <- .era_date(time)
  out <- numeric(length(time))
  for (i in seq_along(time)) {
    m <- regmatches(trimws(time[i]),regexec('^([+-]?[0-9]+)/([0-9]{1,2})/([0-9]{1,2})([[:space:]]+([0-9]{1,2}):([0-9]{2})(:([0-9]{2}(\\.[0-9]+)?))?)?$',trimws(time[i])))[[1]]
    if (!length(m)) stop('Use YYYY/MM/DD with optional HH:MM or HH:MM:SS.',call.=FALSE)
    y <- as.numeric(m[2]); mo <- as.numeric(m[3]); d <- as.numeric(m[4])
    h <- if(nzchar(m[6])) as.numeric(m[6]) else 0
    mi <- if(nzchar(m[7])) as.numeric(m[7]) else 0
    sec <- if(nzchar(m[9])) as.numeric(m[9]) else 0
    if (!all(is.finite(c(y,mo,d,h,mi,sec))) || abs(y)>100000 || mo<1 || mo>12 || d<1 || d>31 || h>23 || mi>59 || sec>=60)
      stop('Invalid calendar date or clock time; leap seconds are not accepted.',call.=FALSE)
    midnight <- swephR::swe_julday(y,mo,d,0,cal)
    back <- swephR::swe_revjul(midnight,cal)
    if (!identical(as.numeric(c(back$year,back$month,back$day)),c(y,mo,d)))
      stop('Invalid date in the selected calendar.',call.=FALSE)
    out[i] <- .local_to_ut(midnight+(h*3600+mi*60+sec)/86400,spec)
  }
  if (verbose) .sky_message('Converted local calendar time to a UT Julian date.')
  out
}

#' Convert a UT Julian date to local calendar date and time
#' @description Express an ephemeris instant as a calendar date and clock time suitable for comparison
#' with an observation or historical record.
#' @details An instant is specified by its UT Julian date. To label it consistently, choose the
#' calendar and the offset or civil timezone appropriate to the record. The conversion
#' rounds to the nearest second, applies the offset at the rounded instant, handles day
#' rollover and converts the day number with Swiss Ephemeris.
#'
#' Fixed offsets may be numeric hours east of UTC (for example 5.5), or
#' strings such as \code{"UTC+05:30"}, \code{"UTC-02"} or \code{"GMT+01"}.
#' Offsets must lie between -14 and +14 hours. UTC, GMT and the empty string
#' mean zero offset; an empty string does not select the system timezone.
#'
#' Named zones listed by \code{OlsonNames()}, such as \code{"Europe/London"},
#' use the installed timezone database for the actual Gregorian date, not a
#' surrogate year. They are supported only for Gregorian years 1900-9999.
#' Earlier dates require UTC or a fixed offset. For ancient-date workflows,
#' set \code{skyscapeR.vars(timezone="UTC")} or provide timezone explicitly.
#' The package's initial timezone remains Europe/London. Database versions
#' may differ across systems; fixed offsets provide reproducible conversions.
#'
#' Gregorian and Julian calendars are proleptic: no automatic calendar reform
#' is applied. Astronomical year numbering is used: 0 is 1 BCE, -1 is 2 BCE.
#' Fixed-offset input dates support years from -100000 to 100000.
#' A timezone offset represents civil clock time, not local solar time.
#' These functions apply clock offsets to a UT Julian day; they do not apply
#' leap-second, UTC-to-UT1 or TT corrections. Leap-second input is rejected.
#'
#' Output is rounded to the nearest whole second, including date rollover.
#' The result does not contain a timezone or UTC offset. During a repeated
#' local hour, two distinct instants can have identical output strings; retain
#' the original Julian date or use a fixed offset for an unambiguous round trip.
#' @param jd Finite numeric vector of UT Julian dates.
#' @param timezone Timezone name or fixed offset; defaults to the value in \code{skyscapeR.vars()}. See Details.
#' @param calendar Output calendar: G or Gregorian, J or Julian. Defaults to the package calendar setting, initially Gregorian.
#' @param verbose Logical; display a conversion message. Default is FALSE.
#' @return A character vector in YYYY/MM/DD HH:MM:SS format. Years are signed as needed and are not padded to four digits.
#' @import swephR
#' @section Validation:
#' Build 2.0.9005 passed all 499 calendar, clock, era, angle and setting checks in A01. The largest discrepancy across 240 absolute Julian dates was 0.000021 seconds. See \href{https://stonexsky.github.io/skyscapeR/manual/#a01-calendar-checks}{A01 validation in the technical manual}.
#' @export
#' @seealso \code{\link{time2jd}}, \code{\link{timestring}}, \code{\link{skyscapeR.vars}}
#' @examples
#' jd <- time2jd("2019/12/31 23:30:00", "UTC")
#' jd2time(jd, "UTC+02") # 2020/01/01 01:30:00
#' jd2time(jd, "Europe/Paris")
#' jd <- time2jd("-1/12/31 23:30:00", "UTC", calendar="J")
#' jd2time(jd, "UTC+02", calendar="J") # 0/01/01 01:30:00
jd2time <- function(jd, timezone, calendar, verbose=F) {
  if (missing(timezone)) timezone <- skyscapeR.env$timezone
  if (missing(calendar)) calendar <- skyscapeR.env$calendar
  cal <- .calendar_flag(calendar); spec <- .timezone_spec(timezone)
  if (!is.numeric(jd) || any(!is.finite(jd))) stop('jd must be a finite numeric vector.',call.=FALSE)
  out <- character(length(jd))
  for (i in seq_along(jd)) {
    # Round the instant before selecting its offset, so DST boundaries agree.
    instant <- floor(jd[i]-.5)+.5+round((jd[i]-.5-floor(jd[i]-.5))*86400)/86400
    offset <- if(!is.null(spec$offset)) spec$offset else .named_offset(instant,spec$zone)
    local_jd <- instant+offset/86400
    day_jd <- floor(local_jd-.5)+.5
    seconds <- round((local_jd-day_jd)*86400)
    if (seconds>=86400) { day_jd <- day_jd+1; seconds <- 0 }
    date <- swephR::swe_revjul(day_jd,cal)
    out[i] <- sprintf('%s/%02d/%02d %02d:%02d:%02d',as.character(date$year),
      date$month,date$day,seconds %/% 3600,(seconds %% 3600) %/% 60,seconds %% 60)
  }
  if (verbose) .sky_message('Converted UT Julian date to local calendar time.')
  out
}


#' Assemble a timestamp from calendar and clock fields
#' @description Assemble separately recorded date and clock fields into the timestamp used by the
#' package's ephemeris routines.
#' @details Survey records often store year, month, day and clock time separately. This formatter
#' pads the month, day and clock components to two characters and joins them as
#' year/month/day and hours:minutes:seconds. Omitted time is noon. Formatting alone does
#' not supply a timezone or check that the date exists; use time2jd to validate and convert
#' it.
#' @seealso \code{\link{time2jd}}, \code{\link{jd2time}}, \code{\link{BC.AD}}, \code{\link{long.date}}
#' @return A character date-time string in the slash-separated format accepted by \code{\link{time2jd}}.
#' @param year Year. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param month Month
#' @param day Day
#' @param hour Hour
#' @param minute Minute
#' @param second Second
#' @section Validation:
#' Build 2.0.9005 passed all 499 calendar, clock, era, angle and setting checks in A01. The largest discrepancy across 240 absolute Julian dates was 0.000021 seconds. See \href{https://stonexsky.github.io/skyscapeR/manual/#a01-calendar-checks}{A01 validation in the technical manual}.
#' @export
#' @examples
#' timestring(2018, 12, 25, 2, 34)
timestring <- function(year, month, day, hour=12, minute=0, second=0) {
  if (is.character(year)) year <- astroYear(year)
  month <- as.character(month); if (nchar(month) < 2) { month <- paste0('0', month) }
  day <- as.character(day); if (nchar(day) < 2) { day <- paste0('0', day) }
  hour <- as.character(hour); if (nchar(hour) < 2) { hour <- paste0('0', hour) }
  minute <- as.character(minute); if (nchar(minute) < 2) { minute <- paste0('0', minute) }
  second <- as.character(second); if (nchar(second) < 2) { second <- paste0('0', second) }

  ts <- paste(paste(year, month, day, sep='/'), paste(hour, minute, second, sep=':'))
  return(ts)
}



#' Format a month and day as a calendar label
#' @description Format a month and day for reporting recurring solar or stellar events within a year.
#' @details A season or event table may need a day/month label without an epoch. The first two
#' characters are the month number and characters four and five the day number; the
#' separator occupies position three. The month is looked up in R's English month
#' abbreviations. This formatting step does not validate a calendar date or convert between
#' calendars.
#' @seealso \code{\link{time2jd}}, \code{\link{jd2time}}, \code{\link{timestring}}, \code{\link{BC.AD}}
#' @return Character month-and-day labels such as "21 Jun", derived from the supplied month and day.
#' @param date date in 'MM-DD' format
#' @section Validation:
#' Build 2.0.9005 passed all 499 calendar, clock, era, angle and setting checks in A01. The largest discrepancy across 240 absolute Julian dates was 0.000021 seconds. See \href{https://stonexsky.github.io/skyscapeR/manual/#a01-calendar-checks}{A01 validation in the technical manual}.
#' @export
#' @examples
#' long.date('01-01')
#' long.date('08-23')
long.date <- function(date){
  day <- as.numeric(substr(date,4,5))
  month <- month.abb[as.numeric(substr(date,1,2))]

  return( paste(day, month) )
}



#' @noRd
findWS <- function(year, calendar) {
  if (is.character(year)) year <- astroYear(year)
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }

  jd0 <- swephR::swe_julday(year, 1, 1, 12, 1)
  jd <- seq(jd0, jd0+365, length.out = 365*24)
  dec <- .file_position(jd, 0, 2048)$xx[,2]
  ind <- which.min(dec)

  aux <- jd2time(jd[ind], calendar=calendar)
  ind2 <- which(strsplit(aux, "")[[1]]==" ")

  return(list(ind=round(ind/24,0)+1, date=substr(aux,1,ind2-1), jd=jd[ind]))
}



#' @noRd
calWS <- function(WS) {
  days <- seq(1,365,1)
  days <- 365-WS$ind+days
  days[days>365] <- days[days>365]-365
  return(days)
}



#' @noRd
dd.to.DD_unvec <- function(day, char=F, WS=F) {
  if (WS==T) { day <- day-11 }
  if (day<=0) { day <- day+ 365 }
  if (char==F) {
    out <- matrix(NA, nrow=length(day), ncol=2)
  } else { out <- c() }

  for (i in 1:length(day)) {
    if (day[i] >= 1 & day[i] <= 365) {
      months <- c(31,28,31,30,31,30,31,31,30,31,30,31)
      summon <- cumsum(months)

      mm <- tail(which(summon < day[i]),1) + 1
      dd <- day[i] - summon[tail(which(summon <= day[i]),1)]
      if (NROW(mm) < 1) {
        mm <- 1
        dd <- day[i]
      }
      if (dd == 0) { dd <- months[mm] }

      if (char) {
        out[i] <- paste0(month.name[mm], ' ',dd)
      } else { out[i,] <- c(mm,dd) }
    }
  }
  return(out)
}
#' @noRd
dd.to.DD <- Vectorize(dd.to.DD_unvec, 'day')


#' Convert an astronomical year to a BC or AD label
#' @description Label an astronomical year in historical BC or AD notation when reporting an epoch.
#' @details Astronomical dating includes a year zero, whereas historical BC and AD labels run
#' directly from 1 BC to AD 1. Positive years retain their number and receive AD. For zero
#' or a negative year, subtract the astronomical year from one and append BC. The
#' conversion works element by element and changes only the label, not the calendar or
#' instant.
#' @seealso \code{\link{time2jd}}, \code{\link{jd2time}}, \code{\link{timestring}}, \code{\link{long.date}}
#' @return Character labels for astronomical years. Year 0 is 1 BC and year -1 is 2 BC.
#' @param year Numeric astronomical year or vector of years. Year zero is 1 BC. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @section Validation:
#' Build 2.0.9005 passed all 499 calendar, clock, era, angle and setting checks in A01. The largest discrepancy across 240 absolute Julian dates was 0.000021 seconds. See \href{https://stonexsky.github.io/skyscapeR/manual/#a01-calendar-checks}{A01 validation in the technical manual}.
#' @export
#' @examples
#' BC.AD(100)
#' BC.AD(0)
#' BC.AD(-1)
#' BC.AD(-99)
#' BC.AD(-100)
BC.AD <- function(year) {
  if (is.character(year)) year <- astroYear(year)
  out <- year
  ind <- which(year<=0); out[ind] <- paste(abs(year[ind]-1), 'BC')
  ind <- which(year>0); out[ind] <- paste(abs(year[ind]), 'AD')
  return(out)
}


#' @noRd
day <- function(date, as.char=F) {
  ind1 <- which(strsplit(date, "")[[1]]=="/")
  ind2 <- which(strsplit(date, "")[[1]]==" "); if (length(ind2)==0) { ind2 <- length(date)+1 }
  aux <- substr(date,ind1[2]+1,ind2-1)
  if (!as.char) { aux <- as.numeric(aux) }
  return(aux)
}


#' @noRd
month <- function(date, as.char=F) {
  ind <- which(strsplit(date, "")[[1]]=="/")
  aux <- substr(date,ind[1]+1,ind[2]-1)
  if (!as.char) { aux <- as.numeric(aux) }
  return(aux)
}


#' @noRd
year <- function(date, as.char=F) {
  ind <- which(strsplit(date, "")[[1]]=="/")
  aux <- substr(date,1,ind[1]-1)
  if (!as.char) { aux <- astroYear(aux) }
  return(aux)
}


#' Convert a historical year to astronomical numbering
#' @description Specify the epoch of an observation with historical era labels.
#' @details Astronomical years include zero. Historical BC or BCE years convert
#' by subtracting the historical year from one. AD and CE years keep their value.
#' Era labels are case-insensitive and follow a positive whole year.
#' @param x Numeric astronomical years, character integer years, or historical labels.
#' @return Numeric astronomical years in input order.
#' @seealso \code{\link{BC.AD}}, \code{\link{timestring}}
#' @examples
#' astroYear(c("3700 BC", "1 BCE", "1 AD", "2026 CE"))
#' @section Validation:
#' Build 2.0.9005 passed all 499 calendar, clock, era, angle and setting checks in A01. The largest discrepancy across 240 absolute Julian dates was 0.000021 seconds. See \href{https://stonexsky.github.io/skyscapeR/manual/#a01-calendar-checks}{A01 validation in the technical manual}.
#' @export
astroYear <- function(x) {
  if (is.numeric(x) && all(is.finite(x))) return(x)
  if (!is.character(x) || anyNA(x)) stop("Supply finite years or historical year labels.",call.=FALSE)
  x <- toupper(trimws(x))
  plain <- grepl("^[+-]?[0-9]+$",x)
  historical <- grepl("^[0-9]+[[:space:]]*(BC|BCE|AD|CE)$",x)
  if (any(!plain & !historical)) stop("Use astronomical years or labels such as '3700 BCE'.",call.=FALSE)
  y <- as.numeric(sub("[[:space:]]*[A-Z]+$","",x))
  if (any(!is.finite(y)) || any(historical & y==0)) stop("Era-labelled years must be positive and finite.",call.=FALSE)
  before <- grepl("BCE?$",x)
  y[before] <- 1-y[before]
  y
}

.era_date <- function(time) {
  era <- grepl("^[[:space:]]*[0-9]+[[:space:]]*(BCE?|AD|CE)/",time,ignore.case=TRUE)
  if (any(era)) time[era] <- paste0(astroYear(sub("/.*","",time[era])),sub("^[^/]*","",time[era]))
  time
}
