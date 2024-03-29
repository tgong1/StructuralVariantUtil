% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SV_tertiary_analysis.R
\name{SV_CNV_integration}
\alias{SV_CNV_integration}
\title{Integrate SV and CNV}
\usage{
SV_CNV_integration(
  sampleID,
  SV_data,
  CNV_data,
  overlap_f = NULL,
  bedtools_dir = NULL
)
}
\arguments{
\item{sampleID}{unique identifier of sample}

\item{SV_data}{SV data in bedpe format}

\item{CNV_data}{data frame of CNV segment}

\item{overlap_f}{the fraction of minimum overlap required of CNV segment as a fraction of SV}

\item{bedtools_dir}{bedtools for use}
}
\value{
SV data in bed format
}
\description{
This function integrate deletion and duplication with copy number segments
}
