#' cna_pipeline
#'
#' run the full CNA pipeline
#' 
#' @param BAMT Tumor Bam File (required if COVERAGE_INPUT = 'bam')
#' @param BAMN Normal Bam File (required if COVERAGE_INPUT = 'bam')
#' @param VCF VCF file containing variants (required if ZYGOSITY_INPUT = 'vcf')
#' @param POPFREQ_VCF VCF file containing population frequencies (optional but recommended; used only if ZYGOSITY_INPUT = 'vcf')
#' @param VCF_FIELDS A colon-separated text string encoding relevant fields of the two VCF files, default "AO:DP:Tumor:Normal:MAF".
#' In order, these reflect the fields giving variant read counts, total depth of coverage, name of the tumor sample, name of the normal sample,
#' and the field giving population frequencies.  The first four come from VCF, the fifth comes from POPFREQ_VCF.
#' @param PREFIX Sample prefix used to label all output files (required)
#' @param OUTDIR Output directory (/out by default)
#' @param BEDFILE Path to bedfile containing target regions.
#' @param NP Number of processors to use.
#' @param COVERAGE_INPUT "bam" or "txt", if "bam" then specify BAMT and BAMN, if "txt" then specify COVERAGE_FILE
#' @param ZYGOSITY_INPUT "vcf" or "txt", if "vcf" then specify VCF, if "txt" then specify ZYGOSITY_FILE
#' @param COVERAGE_FILE path to text file (with header row) containing coverage quantifications per targeted region.
#' Required if COVERAGE_INPUT = "txt".  Required column names are Chr, Start, End, MeanCovT, and MeanCovN.
#' Chromosome names must be in {1,...,22,'X','Y'}, and the Chr/Start/End columns should exactly match the lines of BEDFILE.
#' @param ZYGOSITY_FILE path to text file (with header row) containing minimally-filtered variants, their read counts in tumor and normal,
#' and their population frequencies.  This file should include all detected heterozygous SNPs from the normal sample but can include
#' homozygous SNPS, artifacts, etc. as the pipeline does its own filtering.  Required if ZYGOSITY_INPUT = "txt".
#' Required column names are Chr, Pos, Ref, Alt, Var.T, Tot.T, Var.N, and Tot.N.
#' Ref.N and Ref.T may also be specified, otherwise they are assumed to be Tot.N-Var.N and Tot.T-Var.T, respectively.
#' An additional column Pop_Freq containing 1000 genomes or ExAC population frequencies (values between 0 and 1) is strongly recommended,
#' otherwise every variant is assumed to be a common SNP.  Chromosome names must be in {1,...,22,'X','Y'}.
#' @param NORMALPOOLFILE RDS file containing normal pool quantifications for the capture platform of interest
#' @param POLYMORPHIC_REGIONS text file containing a set of regions to mask because they are highly polymorphic
#' @param LOCUS_ANNOTATION_FILE text file with columns Locus,Chr,Start,End against which to annotate results
#' @param EXCLUDE_LOCI whether to exclude certain loci (exclude_loci_hg19) from classification, default TRUE
#' @param GENDER either 'male', 'female', or 'infer' (default = 'infer' which infers the gender from the data)
#' @param config path to text file containing non-default values of parameters
#' @export
cna_pipeline <- function(BAMT="", BAMN="", VCF = "", POPFREQ_VCF = "", VCF_FIELDS = "AO:DP:Tumor:Normal:MAF",
                         PREFIX="", OUTDIR = "/out", BEDFILE="", NP = 2,
                         COVERAGE_INPUT = "bam", ZYGOSITY_INPUT = "vcf", COVERAGE_FILE = "", ZYGOSITY_FILE = "",
                         NORMALPOOLFILE = NULL, POLYMORPHIC_REGIONS = NULL,
                         LOCUS_ANNOTATION_FILE = NULL, EXCLUDE_LOCI = TRUE, GENDER = "infer", config = "") {
    ## config file overrides any arguments given
    if (config != "") {
        if (file.exists(config)) {
            cfg <- read.delim(config, header = FALSE, stringsAsFactors = FALSE); colnames(cfg) <- c("Parameter", "Value")
            for (i in 1:nrow(cfg)) {
                val1 <- cfg$Value[i]
                ## convert to logical or numeric values if appropriate
                if (cfg$Parameter[i] %in% c("NP")) val1 <- as.numeric(cfg$Value[i])
                assign(cfg$Parameter[i], val1)
            }
        }
        rm(i); rm(cfg); rm(val1)
    }
    ## write parameters to output file
    params.list <- as.list(environment())
    params.list <- lapply(params.list,function(x) ifelse(is.null(x),"NULL",x))
    params.df <- data.frame(Parameter=names(params.list), Value = unlist(params.list), stringsAsFactors=F)
    system(paste0("mkdir -p ",OUTDIR))
    write.table(params.df, file = file.path(OUTDIR,"params.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    ## exclude_loci and PARs should have their own segments, exclude_loci should be excluded from classifications
    if (EXCLUDE_LOCI) {
        data(exclude_loci_hg19); data(PAR_coords_hg19)
        EXCEPTIONAL_REGIONS_1 <- rbind(exclude_loci_hg19[,c("Chr","Start","End")],PAR_coords_hg19[,c("Chr","Start","End")])
        EXCEPTIONAL_REGIONS_2 <- exclude_loci_hg19[,c("Chr","Start","End")]
    } else {
        EXCEPTIONAL_REGIONS_1 <- EXCEPTIONAL_REGIONS_2 <- NULL
    }
    ## files produced
    system(paste0("mkdir -p ",OUTDIR,"/processing"))
    coverage_file <- paste0(OUTDIR,"/processing/",PREFIX,".coverage-data.rds")
    coverage_file_exclude_poly <- paste0(OUTDIR,"/processing/",PREFIX,".coverage-data-exclude-poly.rds")
    genderfile_TN <- paste0(OUTDIR,"/processing/",PREFIX,".gender-inferred.txt")
    zygosity_file_full <- paste0(OUTDIR,"/processing/",PREFIX,".zygosity-data-full.rds")
    zygosity_file_het_N <- paste0(OUTDIR,"/processing/",PREFIX,".zygosity-data-het-N.rds")
    zygosity_file_het_T <- paste0(OUTDIR,"/processing/",PREFIX,".zygosity-data-het-T.rds")
    segfile_TN <- paste0(OUTDIR,"/processing/",PREFIX,".segments.TvN.rds")
    segfile_TP <- paste0(OUTDIR,"/processing/",PREFIX,".segments.TvP.rds")
    segfile_NP <- paste0(OUTDIR,"/processing/",PREFIX,".segments.NvP.rds")
    ChrInfo <- get_chrinfo(c(1:22,"X","Y"),numproc = NP)
    ## 1) prepare coverage data
    if (COVERAGE_INPUT == "bam") {
        COV_DF <- multiComputeCoverage(bamlist = as.list(c(BAMT,BAMN)),
                                       bedfile = BEDFILE, padding = 0, numproc = NP)
    }
    if (COVERAGE_INPUT == "txt") {
        COV_DF <- read.delim(COVERAGE_FILE, stringsAsFactors = FALSE)
        if (!(all(c("Chr","Start","End","MeanCovT") %in% colnames(COV_DF)))) {
            stop("COVERAGE_FILE must have column names Chr, Start, End, and MeanCovT.")
        }
        r1 <- function(x) rep(x,nrow(COV_DF))
        if (!("Midpoint" %in% colnames(COV_DF))) COV_DF$Midpoint <- round((COV_DF$Start + COV_DF$End)/2, 0)
        if (!("Region" %in% colnames(COV_DF))) COV_DF$Region <- paste0(COV_DF$Chr,r1(":"),
                                                                       as.character(as.numeric(COV_DF$Start)+1),r1("-"),COV_DF$End)
        if (!("Length" %in% colnames(COV_DF))) COV_DF$Length <- COV_DF$End - COV_DF$Start + 1
        cnames <- c("Region","Chr","Start","End","Midpoint","Length","MeanCovT","MeanCovN")
        cnames <- cnames[which(cnames %in% colnames(COV_DF))]
        COV_DF <- COV_DF[,cnames]
    }
    covdata_list <- prepare_coverage_data(coverage_df = COV_DF,
                                          normalpool_file = NORMALPOOLFILE,
                                          polymorphic_regionsfile = POLYMORPHIC_REGIONS,
                                          numproc = NP, mincov = 30,
                                          gender = GENDER,
                                          gender_inference_file = genderfile_TN)
    Cov_Data_Full <- covdata_list$coverage_df
    Cov_Data_Poly <- covdata_list$coverage_df_poly
    saveRDS(Cov_Data_Full, file = coverage_file)
    saveRDS(Cov_Data_Poly, file = coverage_file_exclude_poly)
    ## 2) prepare zygosity data
    if (ZYGOSITY_INPUT == "vcf") {
        vcf_fields <- unlist(strsplit(VCF_FIELDS,":"))
        V1 <- vcf2zygosity(VCF, var_field = vcf_fields[1], ref_field = NULL, tot_field = vcf_fields[2],
                           samplename_t = vcf_fields[3], samplename_n = vcf_fields[4])
        if (POPFREQ_VCF != "") {
            V2 <- vcf_to_popfreq_df(POPFREQ_VCF, popfreq_field = "MAF")
            if (nrow(V1) > 0) V1$Ind <- 1:nrow(V1) else V1$Ind <- numeric(0)
            Zyg_Data_Full <- merge(V1,V2,by=c("Chr","Pos","Ref","Alt"),all.x=TRUE,all.y=FALSE)
            Zyg_Data_Full <- Zyg_Data_Full[order(Zyg_Data_Full$Ind),]
            Zyg_Data_Full <- Zyg_Data_Full[,c(colnames(V1),"Pop_Freq")]
            Zyg_Data_Full$Ind <- NULL; rownames(Zyg_Data_Full) <- NULL
        } else {
            Zyg_Data_Full$Pop_Freq <- rep(0.5,nrow(Zyg_Data_Full))
        }
    }
    if (ZYGOSITY_INPUT == "txt") {
        Zyg_Data_Full <- read.delim(ZYGOSITY_FILE, stringsAsFactors = FALSE)
        if (!(all(c("Chr","Pos","Ref","Alt","Var.T","Tot.T","Var.N","Tot.N") %in% colnames(Zyg_Data_Full)))) {
            stop("ZYGOSITY_FILE must have column names Chr, Pos, Ref, Alt, Var.T, Tot.T, Var.N, Tot.N.")
        }
        r1 <- function(x) rep(x,nrow(Zyg_Data_Full))
        if (!("Ref.T" %in% colnames(Zyg_Data_Full))) Zyg_Data_Full$Ref.T <- Zyg_Data_Full$Tot.T - Zyg_Data_Full$Var.T
        if (!("Ref.N" %in% colnames(Zyg_Data_Full))) Zyg_Data_Full$Ref.N <- Zyg_Data_Full$Tot.N - Zyg_Data_Full$Var.N
        if (!("Pop_Freq" %in% colnames(Zyg_Data_Full))) Zyg_Data_Full$Pop_Freq <- r1(0.5)
        cnames <- c("Chr","Pos","Ref","Alt","Ref.T","Var.T","Tot.T","Ref.N","Var.N","Tot.N","Pop_Freq")
        cnames <- cnames[which(cnames %in% colnames(Zyg_Data_Full))]
        Zyg_Data_Full <- Zyg_Data_Full[,cnames]
        Zyg_Data_Full$Plot_Coord <- get_plotting_coords(position_vec = Zyg_Data_Full$Pos,
                                                        chr_vec = Zyg_Data_Full$Chr,
                                                        numproc = NP)
    }
    Zyg_Data_Full$Plot_Coord <- get_plotting_coords(position_vec = Zyg_Data_Full$Pos,
                                                    chr_vec = Zyg_Data_Full$Chr,
                                                    numproc = NP)
    saveRDS(Zyg_Data_Full, file = zygosity_file_full)
    ## 3) restrict to heterozygous SNPs
    Zyg_Data_TN <- restrict_to_het(Zyg_Data_Full, mincov = 30, filtertype = "normal",
                                   gender_inference_file = genderfile_TN)
    saveRDS(Zyg_Data_TN, file = zygosity_file_het_N)
    if (!is.null(NORMALPOOLFILE)) {
        Zyg_Data_TP <- restrict_to_het(Zyg_Data_Full, mincov = 30, filtertype = "tumor",
                                       gender_inference_file = genderfile_TN)
        saveRDS(Zyg_Data_TP, file = zygosity_file_het_T)
    }
    ## 4) segmentation: tumor vs normal (Seg_Data_TN) and tumor vs pool (Seg_Data_TP)
    Seg_Data_TN <- integrated_segmentation(coverage_data = Cov_Data_Full, zygosity_data = Zyg_Data_TN,
                                           coverage_colname = "LR.Corrected.TN",
                                           use_normal_sample = TRUE,
                                           exceptional_regions = EXCEPTIONAL_REGIONS_1)
    saveRDS(Seg_Data_TN, file = segfile_TN)
    if (!is.null(NORMALPOOLFILE)) {
        Seg_Data_TP <- integrated_segmentation(coverage_data = Cov_Data_Full, zygosity_data = Zyg_Data_TP,
                                               coverage_colname = "LR.Corrected.TP",
                                               use_normal_sample = FALSE,
                                               exceptional_regions = EXCEPTIONAL_REGIONS_1)
        Seg_Data_NP <- coverage_only_segmentation(coverage_data = Cov_Data_Full,
                                                  coverage_colname = "LR.Corrected.NP", chrinfo = ChrInfo)
        saveRDS(Seg_Data_TP, file = segfile_TP)
        saveRDS(Seg_Data_NP, file = segfile_NP)
    }
    ## 5) classification
    gender_N <- readLines(genderfile_TN)[2]
    integrated_classification(covdata = Cov_Data_Full, zygdata = Zyg_Data_TN, segdata = Seg_Data_TN,
                              outdir = OUTDIR, refsample = "normal", prefix = PREFIX, numproc = NP,
                              locus_annotation_file = LOCUS_ANNOTATION_FILE,
                              exceptional_regions = EXCEPTIONAL_REGIONS_2, gender = gender_N)
    if (!is.null(NORMALPOOLFILE)) {
        gender_T <- readLines(genderfile_TN)[1]
        integrated_classification(covdata = Cov_Data_Full, zygdata = Zyg_Data_TP, segdata = Seg_Data_TP,
                                  outdir = OUTDIR, refsample = "pool", prefix = PREFIX, numproc = NP,
                                  locus_annotation_file = LOCUS_ANNOTATION_FILE,
                                  exceptional_regions = EXCEPTIONAL_REGIONS_2, gender = gender_T)
        germline_classification(covdata = Cov_Data_Full, segdata = Seg_Data_NP, outdir = OUTDIR,
                                prefix = PREFIX, chrinfo = ChrInfo, gene_df = GENE_STRUCTURE)
        ## 6) make germline zygosity plot for QC purposes
        make_germline_zygosity_plot(zygdata = Zyg_Data_Full,
                                    polymorphic_regionsfile = POLYMORPHIC_REGIONS,
                                    outdir = OUTDIR, prefix = PREFIX, chrinfo = ChrInfo)
    }
    system(paste0("touch ",OUTDIR,"/run_finished.txt"))
}
