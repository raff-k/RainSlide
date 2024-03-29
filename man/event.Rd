% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/event.R
\name{event}
\alias{event}
\title{Get rain event of a landslide}
\usage{
event(
  x,
  dates = NULL,
  timesteps = NULL,
  date.of.failure = NULL,
  sub.RainEvent = TRUE,
  all.RainEvent = FALSE,
  cumu.RainFall = NULL,
  return.DataFrame = TRUE,
  S1.rainThresh = 0.2,
  S3.rainThresh = 1,
  S1.rainOffLength = c(3, 6),
  S2.rainOffLength = c(6, 12),
  S4.rainOffLength = c(48, 96),
  RD = NULL,
  MAP = NULL,
  RDN = MAP/RD,
  index.month.warm.season = c(4, 10),
  force.limit = NULL
)
}
\arguments{
\item{x}{vector containing precipitation}

\item{dates}{vector containing dates. The length of dates must be similar to the length of x. Default: NULL}

\item{timesteps}{time period or length of observation. The rev(x) and rev(dates) are subsetted to this length according to date.of.failure if set. I.e. 24 for hourly or 1 for daily data. Default: NULL}

\item{date.of.failure}{date of failure. If set data is subsetted to this date. Must be of class "POSIXct" or "POSIXt". Default: NULL}

\item{sub.RainEvent}{examine potential sub-rain-events of critical rainfall event. Default: TRUE}

\item{all.RainEvent}{if TRUE, all rain events in data are extracted. By setting this option, no critical rain event metrics are computed. Default: FALSE}

\item{cumu.RainFall}{vector containing time intervals for cumulative rainfall. I.e. c(24, 48, 96) for 1, 2 and 4 days aggregation. Default: NULL}

\item{return.DataFrame}{only the rain events are returned as a data.frame. Default: TRUE}

\item{S1.rainThresh}{isolated rainfall measurements below this thresholds are removed in the first step. Default: 0.2}

\item{S3.rainThresh}{exclusion of irrelevant rainfall sub-events under and equal to this threshold (third step). Default: 1 (in mm)}

\item{S1.rainOffLength}{dry periods between isolated rain events in the first step. Default: c(3, 6) (hours). When dates is NULL, then the smallest values is used for separation.}

\item{S2.rainOffLength}{dry periods between rainfall sub-events in the second step. Default: c(6, 12) (hours). When dates is NULL, then the smallest values is used for separation.}

\item{S4.rainOffLength}{dry periods between rainfall sub-events in the second step. Default: c(48, 96) (hours). When dates is NULL, then the smallest values is used for separation.}

\item{RD}{average number of rainy days in a year, proxy for locate climate conditions. Default: NULL}

\item{MAP}{mean annual precipitation, the long-term yearly average precipitation, see CRU - climate research units for number. Default: NULL}

\item{RDN}{a climatic index that provides better description (or proxy) than the MAP for the occurrence of extreme storm events (Guzzetti et al. 2006: 247). Default: MAP/RD}

\item{index.month.warm.season}{month indices of the warm season. First element is start, and second element represents the end (all including). Only relevant when dates are set. Default: c(4, 10) (including April, including October)}

\item{force.limit}{Useful for standard output, if result is type list. Must be integer number of specific size. If result is smaller, than results gets filled with NA, otherwise cut to this size. Default: NULL}
}
\value{
vector containing rainfall metrics (see description). If return.DataFrame is TRUE a data.frame is returned containing similar
rain metrics for all rain events.
}
\description{
This function calculates different precipitation characteristics for a specific time-series:
total precipitation, number of rainfall events, weighted mean intensity of rainfall events (normalized by MAP, RD or RDN),
cumulative critical event rainfall (normalized by MAP, RD or RDN), maximum rainfall during critical rainfall event,
duration of critical rainfall event, critical rainfall intensity (normalized by MAP, RD or RDN), rainfall at day of failure (start date),
rainfall intensity at day of failure (start date), maximum rainfall at day of failure (start date).
}
\note{
\itemize{
\item thresholds are oriented at hourly data
\item timesteps with precipitation equal 0 are included (see Melillo et al. 2015: 314)
\item Guzzetti, F., Peruccacci, S., Rossi, M., & Stark, C. P. (2007). Rainfall thresholds for the initiation of landslides in central and southern Europe. Meteorology and atmospheric physics, 98(3-4), 239-267.
\item Rossi, M., Luciani, S., Valigi, D., Kirschbaum, D., Brunetti, M. T., Peruccacci, S., & Guzzetti, F. (2017). Statistical approaches for the definition of landslide rainfall thresholds and their uncertainty using rain gauge and satellite data. Geomorphology, 285, 16-27.
\item Melillo, M., Brunetti, M. T., Peruccacci, S., Gariano, S. L., & Guzzetti, F. (2015). An algorithm for the objective reconstruction of rainfall events responsible for landslides. Landslides, 12(2), 311-320.
}
}
\keyword{event,}
\keyword{landslide,}
\keyword{metrics}
\keyword{rainfall}
\keyword{tresholds,}
