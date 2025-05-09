% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/looping_machine.R
\name{looping_machine}
\alias{looping_machine}
\title{Run dynamic modelling I-ARIMAX algorithm.}
\usage{
looping_machine(
  dataframe,
  a_series,
  b_series,
  c_series,
  id_var,
  min_n_subject = 20,
  minvar = 0.01,
  hlm_compare = TRUE,
  metaanalysis = TRUE,
  timevar = NULL,
  weight_rma = FALSE,
  weight_rma_var = NULL,
  correlation_method = "pearson"
)
}
\arguments{
\item{dataframe}{Your dataframe.}

\item{a_series}{The first part of the loop.}

\item{b_series}{The second part of the loop.}

\item{c_series}{The third part of the loop.}

\item{id_var}{A string containing your id variable.}

\item{min_n_subject}{The minimum number of non NA cases to run the analyses. It will filter cases with more NA's than the threshold. Defaults to 20.}

\item{minvar}{The minimum variance for both series (&) to include a case. Defaults to 0.01.}

\item{hlm_compare}{Optional, to create a comparison with an HLM model, default is FALSE.}

\item{metaanalysis}{Bool to run a random effects meta-analysis or not.}

\item{timevar}{If hlm_compare is TRUE, then a time variable is needed, default is NULL.}

\item{weight_rma}{If adding an exogenous weight variable to the RMA model.}

\item{weight_rma_var}{Select weight RMA variable. Defaults as NULL, if NULL then the number of valid observations for Y AND X will be used (!is.na).}

\item{correlation_method}{Select method for raw semi-partial correlations. Options are: 'spearman', 'pearson' or 'kendall'. Defaults to 'pearson'.}
}
\value{
a list with dataframes from original iarimax algorithm outuputs + looping machine dataframe.
}
\description{
Run dynamic modelling I-ARIMAX algorithm.
}
