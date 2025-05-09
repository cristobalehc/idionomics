% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/i_standarbot_300.R
\name{i_standarbot_300}
\alias{i_standarbot_300}
\title{Person mean standardize time series.}
\usage{
i_standarbot_300(df, cols, idvar, explanation = TRUE, append = TRUE)
}
\arguments{
\item{df}{A dataframe.}

\item{cols}{A vector of strings with the variable names you would like to standardize.}

\item{idvar}{Your ID variable. It currently accepts only one.}

\item{explanation}{Boolean to print explanation of the function (default is TRUE).}

\item{append}{Boolean to append new standardized columns to your original dataframe or not (default is TRUE)}
}
\value{
Will return a dataframe with new person-mean standardized columns appended at the end.
}
\description{
Person mean standardize time series.
}
