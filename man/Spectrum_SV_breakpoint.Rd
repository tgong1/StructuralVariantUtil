% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SV_tertiary_analysis.R
\name{Spectrum_SV_breakpoint}
\alias{Spectrum_SV_breakpoint}
\title{Genomic bins hotspots}
\usage{
Spectrum_SV_breakpoint(
  All_sampleID,
  SVdf_list,
  threshold_count_breakpoint = NULL,
  threshold_count_sample = NULL
)
}
\arguments{
\item{All_sampleID}{sample ID for all samples}

\item{SVdf_list}{list of SVs in df}

\item{threshold_count_breakpoint}{threshold of number of SD}

\item{threshold_count_sample}{threshold of number of samples}
}
\value{
data frame of genomic bins with hotspots defined
}
\description{
Define hotspots
}
