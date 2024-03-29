% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/routines.R
\name{log_likCInd2_PGW}
\alias{log_likCInd2_PGW}
\title{MEGH II}
\usage{
log_likCInd2_PGW(par, u, times, status, des, des_t, index)
}
\arguments{
\item{u}{: random effects}

\item{times}{: survival times}

\item{status}{: vital status}

\item{des}{: design matrix (p x n)}

\item{des_t}{: design matrix for time-dependent effects (q x n)}

\item{index}{: individual ID}

\item{par:}{all model parameters in the original parameterisation}
}
\value{
Conditional individual log-likelihood with PGW baseline hazard
}
\description{
MEGH II
}
