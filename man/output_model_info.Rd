% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cn-fit-functions.R
\name{output_model_info}
\alias{output_model_info}
\title{output_model_info}
\usage{
output_model_info(coverage_data, coverage_data_col, zygosity_data, adjustment,
  chrinfo, classification_df, model_outdir, prefix, tc, locus_df = NULL,
  exceptional_regions = NULL, gender = "female")
}
\arguments{
\item{coverage_data}{data.frame containing coverage}

\item{coverage_data_col}{column name to use for coverage}

\item{zygosity_data}{data.frame containing zygosity}

\item{adjustment}{ploidy adjustment}

\item{chrinfo}{chromosome information data.frame}

\item{classification_df}{data.frame of classifications}

\item{model_outdir}{directory to which to write output files for model}

\item{prefix}{sample prefix to use for all result files}

\item{tc}{tumor content value for the model}

\item{locus_df}{data.frame of gene annotations to use for annotating the results (default NULL - no annotations)}

\item{exceptional_regions}{data.frame of exceptional regions}

\item{gender}{gender of sample (default female)}
}
\description{
Output model information
}

