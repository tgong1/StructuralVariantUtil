% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SV_tertiary_analysis.R
\name{Spectrum_SV_type}
\alias{Spectrum_SV_type}
\title{Spectrum of SV count}
\usage{
Spectrum_SV_type(
  All_sampleID,
  SVdf_list,
  identify_hyperSV_tumour = FALSE,
  threshold_total = NULL,
  threshold_relative_freq = NULL
)
}
\arguments{
\item{All_sampleID}{sample ID}

\item{SVdf_list}{list of SV data.frame}

\item{threshold_total}{threshold of minimum total count of SVs per sample}

\item{threshold_relative_freq}{threshold of minimum relative frequency of one SV type}

\item{identify_hyperSV_tumour:}{TRUE or FALSE. Whether to identify hyper-SV mutated tumour samples for large cancer cohort. Default as FALSE.}
}
\value{
data frame of hyper SV
}
\description{
This function read bed format
}
