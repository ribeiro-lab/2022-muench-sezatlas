% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_ts.R
\name{normalize_ts}
\alias{normalize_ts}
\title{normalize time-series data}
\usage{
normalize_ts(
  data,
  bg_frames = 4:9,
  value_column = "value",
  method = "dff",
  regionx = NULL
)
}
\arguments{
\item{data}{input data (format as in data_t)}

\item{bg_frames}{background frames for f0 calculation.}

\item{method}{for now only 'dff'}

\item{region}{region id for "one_region" option}
}
\value{

}
\description{
Normalizes TS data e.g. by calculating DF/F.
}
