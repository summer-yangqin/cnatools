% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zygosity_functions.R
\name{restrict_to_het}
\alias{restrict_to_het}
\title{restrict_to_het}
\usage{
restrict_to_het(zygosity_data, mincov = 30, filtertype = "normal",
  numproc = 1, gender_inference_file = NULL)
}
\arguments{
\item{zygosity_data}{data.frame containing variants}

\item{mincov}{Minimum total coverage to allow in order to keep a variant (default 30)}

\item{filtertype}{"normal" or "tumor" - which sample to use for filtering, default "normal"}

\item{numproc}{Number of processors to use, default 1}
}
\value{
data.frame containing zygosity data for variants that are inferred to be heterozygous
}
\description{
Restrict a zygosity data.frame to SNPs that are inferred to be heterozygous
}

