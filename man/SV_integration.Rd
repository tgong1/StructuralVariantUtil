% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SV_integration_filtration.R
\name{SV_integration}
\alias{SV_integration}
\title{Integrate SV call sets and write output}
\usage{
SV_integration(
  SVCaller_name,
  vcf_list,
  sampleID = "sample",
  bkpt_T_callers = 200,
  SVTYPE_ignore = FALSE,
  bedtools_dir = NULL
)
}
\arguments{
\item{SVCaller_name}{name of callers}

\item{vcf_list}{list of VCF files from different callers}

\item{sampleID}{sample ID}

\item{bkpt_T_callers}{threshold of breakpoint difference}

\item{SVTYPE_ignore}{whether ignore SV type for integration}

\item{bedtools_dir}{directory of bedtools}
}
\value{
data frame
}
\description{
This function read bed format
}
