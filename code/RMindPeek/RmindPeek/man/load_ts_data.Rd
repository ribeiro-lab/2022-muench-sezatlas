% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_ts_data.R
\name{load_ts_data}
\alias{load_ts_data}
\title{Load time series data}
\usage{
load_ts_data(path = "", ds = NA, stim_exclude = NA, stim_include = NA)
}
\arguments{
\item{path}{path containing ts data}

\item{ds}{names of the data sets to load (defaults to "all")}

\item{stim_exclude}{stimuli to exclude (character vector of format "stimulus concentration")}

\item{stim_include}{stimuli to exclude (character vector of format "stimulus concentration")}
}
\value{
data.frame with all ts data cleaned up
}
\description{
Load time series data
}
