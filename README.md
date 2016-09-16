# cnatools
Tools for detection of somatic and germline copy number aberrations from targeted next-generation sequencing data

1. Summary

   This R package jointly estimates tumor purity and somatic copy number across the genome, given a pair of tumor/normal samples sequenced using a targeted capture approach.  It has been successfully used to analyze whole exome data as well as data from relatively large targeted panels (1500-1700 genes).

   Features include the following:
   * Preprocessing of coverage and zygosity data, using both a direct tumor-normal comparison as well as a tumor vs. pooled normal approach
   * Joint purity and copy number estimation using a custom EM Algorithm
   * Ploidy correction by minimization of a lack-of-fit statistic, incorporating constraints to penalize biologically unlikely solutions
   * Germline copy number estimation
   * Coverage-only copy number estimation, for cases in which the sample is a mixture from different individuals (e.g. bone marrow transplants)

   (Temporary Note: tumor vs pooled, germline, and coverage-only approaches are in the process of being finalized and will be documented shortly.)

2. Installation

   The included Dockerfile contains the steps necessary to install the package starting from a clean ubuntu 14.04 installation.

   The code has been tested on R, version 3.3.1, using the following versions of packages (current as of 9/15/2016):
   * data.table (1.9.6)
   * gtools (3.5.0)
   * VGAM (1.0-2)
   * devtools (1.12.0)
   * IRanges (2.6.1)
   * Rsamtools (1.24.0)
   * DNAcopy (1.46.0)
   * Rsubread (1.22.3)
   * BSgenome.Hsapiens.UCSC.hg19 (1.4.0)

   Other than R and these packages (plus their dependencies), the only external dependency is samtools, and it is only required when interfacing with BAM files.

   To test the installation, I have provided artificial data that may be used to test whether the algorithm runs properly.

   Here are the steps necessary to run this test, using the included Dockerfile:
   * download the Dockerfile, and from the directory containing this file, build the docker image:
   ```bash
   docker build --rm=true -t mctp/cnatools .
   ```
   * Run R from within a container:
   ```bash
   docker run --rm=true -it mctp/cnatools R
   ```
   * From within R, load the package and run the cna_pipeline function:
   ```r
   library(cnatools)
   cna_pipeline(COVERAGE_FILE = system.file("data/sample_coverage_data.txt",package="cnatools"),
   			      ZYGOSITY_FILE = system.file("data/sample_zygosity_data.txt",package="cnatools"),
			      COVERAGE_INPUT = "txt", ZYGOSITY_INPUT = "txt",
			      BEDFILE = system.file("data/sample_bedfile.txt",package="cnatools"),
			      NP = 6, OUTDIR = "/out", PREFIX = "testrun")
   ```
   * This should produce output (as described below) in /out, using 'testrun' as a prefix for output filenames.
   * For regular use, one would want to run as a non-root user and write to a filesystem outside of the container; this can be accomplished by adding a 'useradd' line to the end of the Dockerfile, and by mounting directories inside the container using the -v option when running docker.  See the docker documentation for details.

   Non-docker users should follow analogous installation steps, namely installing the necessary prerequisites as described above and running devtools::install_github("mctp/cnatools") to install the package.

3. Running the algorithm on real data

   The main function in the package is the cna_pipeline function, which requires a BED file as well as coverage and zygosity data for the tumor and normal samples.  Specifically:

   * The BED file should contain the coordinates of the regions targeted by the assay. (argument BEDFILE to cna_pipeline)
   * The coverage data may be specified in one of two ways:
     1. A tab-delimited text file containing mean coverage per target region for each line of the BED file (argument COVERAGE_FILE to cna_pipeline), or
     2. Tumor and normal BAM files (arguments BAMT and BAMN to cna_pipeline), in which case the mean coverage statistics are computed.
     (The COVERAGE_INPUT argument, which takes values 'bam' or 'txt', controls which input format is assumed.)
   * The zygosity data may also be specified in one of two ways:
     1. A tab-delimited text file containing variant read count data, and optionally population frequencies, at SNP positions throughout the genome (argument ZYGOSITY_FILE to cna_pipeline)
     2. A pair of VCF files (arguments VCF and POPFREQ_VCF to cna_pipeline, the first file contains read count information and the second contains population frequencies, they could be the same file).
     The VCF_FIELDS argument controls how these files are parsed and contains five values separated by colons: (1) field encoding variant read counts, (2) field encoding read depth, 
     (3) sample name for tumor sample, (4) sample name for normal sample, and (5) field corresponding to population frequencies.  
     The default value is "AO:DP:Tumor:Normal:MAF" which corresponds to the output obtained by running freebayes on a tumor/normal pair,
     each with a single read group (called 'Tumor' and 'Normal', respectively), and annotating against a 1000 genomes VCF using SnpSift.
     (The ZYGOSITY_INPUT argument, which takes values 'txt' or 'vcf', controls which input format is assumed.)
     In either case, this file should contain as many heterozygous SNPs from the normal sample as can be detected by the assay, but does not need to be filtered (the package does its own filtering).

   After installing the package, the following commands should clarify the required format of these three files:
   ```r
   head(read.delim(system.file("data/sample_bedfile.txt",package="cnatools"),stringsAsFactors=FALSE,header=F))
   head(read.delim(system.file("data/sample_coverage_data.txt",package="cnatools"),stringsAsFactors=FALSE,header=F))
   head(read.delim(system.file("data/sample_zygosity_data.txt",package="cnatools"),stringsAsFactors=FALSE,header=F))
   ```

   The cna_pipeline function requires three additional arguments:
   
   * NP (number of processors to use; many steps are multi-threaded using the parallel package in R),
   * OUTDIR (output directory to which to write results, the user is assumed to have the necessary permissions to create and write to this directory)
   * PREFIX (a sample prefix to affix to the beginning of each output file).

   It is possible to submit arguments via a configuration file (a two-column tab-delimited text file, no header, with the first column containing parameter names, and the second column containing parameter values).
   In this case, the call to cna_pipeline would simply be cna_pipeline(config = "/path/to/config.file").  All arguments not appearing in this file would take on their default values.

4. Output Format

   The output directory is structured as follows:

   * classification_models_TN
     * model_1
       * PREFIX.adjustment.txt
       * PREFIX.classified-segments.txt
       * PREFIX.copynumber_by_locus.txt
       * PREFIX.coverage_data_classified.txt
       * PREFIX.coverage_zygosity_plot.png
       * PREFIX.tumorcontent.txt
       * PREFIX.zygosity_data_classified.txt	
     * model-summary-stats.txt
     * PREFIX.QC.TN.txt
   * params.txt
   * processing
     * PREFIX.coverage-data-exclude-poly.rds
     * PREFIX.coverage-data.rds
     * PREFIX.gender-inferred.txt
     * PREFIX.segments.TvN.rds
     * PREFIX.zygosity-data-full.rds
     * PREFIX.zygosity-data-het-N.rds
   * run_finished.txt

   At the beginning of the run, a params.txt file is produced that contains the set of input parameters passed to the algorithm.
   The processing subdirectory contains the processed coverage, zygosity, and segmentation data, stored as binary R files (they can be read via the readRDS function in R).
   The classification_models_TN directory (TN here stands for Tumor vs. Normal) contains the output from a set of candidate ploidy adjustments (each contained in its own folder: model_1, model_2, ..., model_N).
   These models are ranked by plausibility; model_1 has been ranked as the best, but output from the other models is retained in case the algorithm makes an incorrect choice.
   For convenience, we have always appended a final subdirectory (model_N) that is not really a model at all, it simply produces output files corresponding to the case where the sample has extremely low purity and all copy number classifications are copy neutral.

   An analogous classification_models_TP folder will be produced once the tumor vs. pool code is complete.

   The relevant output files in each model_i folder are as follows:
     * adjustment.txt: the adjustment made to the coverage data to compensate for the sample potentially not being diploid
     * classified-segments.txt: final classifications in a tab-delimited text file, one line per segment.
     * copynumber_by_locus.txt: copy number summarized by locus, if an input locus file is given to the cna_pipeline function
     * coverage_data_classified.txt: coverage data per targeted region, along with each region's copy number classification
     * zygosity_data_classified.txt: zygosity data per SNP, along with each SNP's final copy number classification
     * tumorcontent.txt: final estimated tumor content (between 0 and 1)
     * coverage_zygosity_plot.png: plot of coverage and zygosity data, colored by status (dark red = amplification, red = gain, black = neutral, blue = loss, dark blue = homozogyous loss, orange = region excluded from classification, e.g. T-cell receptor loci)

   The QC.TN.txt file contains some useful statistics relevant to the lower limit of tumor content estimation, and the model-summary-stats.txt file contains some summary statistics used to rank the models.  These will be described in more detail later on.

   At the end of a run, the run_finished.txt touch file is produced.

5. Notes and Limitations

   *  Throughout, chromosomes are assumed to be named {1,2,...,22,X,Y}.  If this is incompatible with the user's BAM and VCF files, by specifying COVERAGE_INPUT = 'txt' and ZYGOSITY_INPUT = 'txt' the user can produce files with these chromosome names.
   *  Samples are assumed to be human - xenograft samples may work but this has not been tested.
   *  Coordinates are assumed to be with respect to GRCh37.  In particular, locations of T-cell receptors, IgH loci, and PAR regions are currently hard-coded into the code using GRCh37 coordinates.  At some point I will attempt to support GRCh38 coordinates, but this effort is not underway yet.

