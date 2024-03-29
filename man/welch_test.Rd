% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WWER.R
\name{welch_test}
\alias{welch_test}
\title{Simple implementation of a welch test.}
\usage{
welch_test(beta1, se1, beta2, se2, welch_thresh = 0.05)
}
\arguments{
\item{beta1}{Float or vector of floats, mean of the first sample.}

\item{se1}{Float or vector of floats, SD of estimate of mu1.}

\item{beta2}{Float or vector of floats, mean of the second sample.}

\item{se2}{Float or vector of floats, SD of estimate of mu2.}

\item{welch_thresh}{Float, p_value threshold for significance.}
}
\description{
This tests the null hypothesis abs(beta1) = abs(beta2) against the two
alternatives abs(beta1) > abs(beta2) and abs(beta1) < abs(beta2).
}
