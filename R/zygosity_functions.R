#' vcf_to_df
#'
#' Convert a VCF file to a data.frame
#'
#' @param vcf Input VCF file (required)
#' @return data.frame
#'
#' @export
vcf_to_df <- function(vcf) {
    filetype <- ifelse(substr(vcf, nchar(vcf) - 2, nchar(vcf)) == ".gz","vcf.gz","vcf")
    if (filetype == "vcf.gz") { con <- gzfile(vcf); A <- readLines(con); close(con) }
    if (filetype == "vcf") A <- readLines(vcf)
    header_cols <- unlist(strsplit(A[which(substr(A,1,1) == "#" & substr(A,1,2) != "##")],"\t"))
    A <- A[which(substr(A,1,1) != "#")]
    A <- as.data.frame(do.call(rbind,strsplit(A,"\t")),stringsAsFactors = F)
    colnames(A) <- header_cols
    return(A)
}

#' vcf2zygosity
#'
#' Convert a VCF file to a data.frame containing zygosity data for tumor/normal files
#'
#' @param input_vcf VCF file (may be gzipped) containing variants, variant IDs, and read counts
#' @param var_field Field containing per-variant variant read counts (e.g. "AO")
#' @param ref_field Field containing per-variant reference read counts (e.g. "RO")
#' @param tot_field Field containing per-variant total read counts (e.g. "DP")
#' @param samplename_t Sample name under which tumor sample appears (NULL if no tumor sample)
#' @param samplename_n Sample name under which normal sample appears (NULL if no normal sample)
#' @return data.frame with columns ID, Chr, Pos, Ref, Alt, Ref.T, Var.T, Tot.T, Ref.N, Var.N, Tot.N
#'
#' @export
vcf2zygosity <- function(input_vcf, var_field = NULL, ref_field = NULL, tot_field = NULL,
                           samplename_t = "Tumor", samplename_n = "Normal") {
    samples_to_process <- c(samplename_t, samplename_n)
    if (is.null(samples_to_process)) stop("No samples to process, exiting.\n")
    ## read VCF files
    if (!file.exists(input_vcf)) stop("VCF file does not exist, exiting.\n")
    fields_to_process <- c(var_field, ref_field, tot_field)
    if (length(fields_to_process) < 2) stop("At least 2 of var_field, ref_field, and tot_field must be provided.\n")
    V1_df <- vcf_to_df(input_vcf)
    ## retrieve variant / reference / total counts
    if (!is.null(samplename_t)) {
        if (!(samplename_t %in% colnames(V1_df))) {
            stop("provided tumor sample name does not exist in VCF file, exiting.\n")
        }
        sT <- strsplit(V1_df[,samplename_t],":")
        stmp <- sT
    }
    if (!is.null(samplename_n)) {
        if (!(samplename_n %in% colnames(V1_df))) {
            stop("provided normal sample name does not exist in VCF file, exiting.\n")
        }
        sN <- strsplit(V1_df[,samplename_n],":")
        stmp <- sN
    }
    for (i in 1:length(stmp)) {
        u <- unlist(strsplit(V1_df$FORMAT[i],":"))
        if (!is.null(samplename_t)) {
            if (length(sT[[i]])==1) {
                if (sT[[i]][1]==".") {
                    sT[[i]] <- rep(NA,length(u))
                }
            }
            names(sT[[i]]) <- u
        }
        if (!is.null(samplename_n)) {
            if (length(sN[[i]])==1) {
                if (sN[[i]][1]==".") {
                    sN[[i]] <- rep(NA,length(u))
                }
            }
            names(sN[[i]]) <- u
        }
    }
    for (cname in fields_to_process) {
        if (!is.null(samplename_t)) {
            x <- unlist(lapply(sT,function(x) x[cname]))
            x <- ifelse(x == ".","0", x)
            V1_df$tmp <- as.numeric(x)
            if (!is.null(var_field)) {
                if (cname == var_field) colnames(V1_df)[which(colnames(V1_df) == "tmp")] <- "Var.T"
            }
            if (!is.null(ref_field)) {
                if (cname == ref_field) colnames(V1_df)[which(colnames(V1_df) == "tmp")] <- "Ref.T"
            }
            if (!is.null(tot_field)) {
                if (cname == tot_field) colnames(V1_df)[which(colnames(V1_df) == "tmp")] <- "Tot.T"
            }
        }
        if (!is.null(samplename_n)) {
            x <- unlist(lapply(sN,function(x) x[cname]))
            x <- ifelse(x == ".","0", x)
            V1_df$tmp <- as.numeric(x)
            if (!is.null(var_field)) {
                if (cname == var_field) colnames(V1_df)[which(colnames(V1_df) == "tmp")] <- "Var.N"
            }
            if (!is.null(ref_field)) {
                if (cname == ref_field) colnames(V1_df)[which(colnames(V1_df) == "tmp")] <- "Ref.N"
            }
            if (!is.null(tot_field)) {
                if (cname == tot_field) colnames(V1_df)[which(colnames(V1_df) == "tmp")] <- "Tot.N"
            }
        }
    }
    if (is.null(var_field)) {
        if (!is.null(samplename_t)) V1_df$Var.T <- V1_df$Tot.T - V1_df$Ref.T
        if (!is.null(samplename_n)) V1_df$Var.N <- V1_df$Tot.N - V1_df$Ref.N
    }
    if (is.null(ref_field)) {
        if (!is.null(samplename_t)) V1_df$Ref.T <- V1_df$Tot.T - V1_df$Var.T
        if (!is.null(samplename_n)) V1_df$Ref.N <- V1_df$Tot.N - V1_df$Var.N
    }
    if (is.null(tot_field)) {
        if (!is.null(samplename_t)) V1_df$Tot.T <- V1_df$Var.T + V1_df$Ref.T
        if (!is.null(samplename_n)) V1_df$Tot.N <- V1_df$Var.N + V1_df$Ref.N
    }
    ## only keep positions at which ref + alt together account for at least 95% of the total coverage
    if (!is.null(samplename_t)) V1_df <- subset(V1_df, Var.T + Ref.T >= 0.95 * Tot.T)
    if (!is.null(samplename_n)) V1_df <- subset(V1_df, Var.N + Ref.N >= 0.95 * Tot.N)
    V1_df$Chr <- V1_df$"#CHROM"
    V1_df$Pos <- as.numeric(V1_df$POS)
    V1_df$Ref <- V1_df$REF
    V1_df$Alt <- V1_df$ALT
    cnames <- c("ID","Chr","Pos","Ref","Alt","Ref.T","Var.T","Tot.T","Ref.N","Var.N","Tot.N")
    V1_df <- V1_df[,intersect(cnames, colnames(V1_df))]
    V1_df <- subset(V1_df,nchar(Ref) == nchar(Alt) & nchar(Ref) <= 2) ## only single or dinucleotide variants
    V1_df <- V1_df[order(factor(V1_df$Chr, levels = c(1:22,"X","Y")), V1_df$Pos), ]
    rownames(V1_df) <- NULL
    return(V1_df)
}

#' vcf_to_popfreq_df
#'
#' Convert a VCF containing population frequencies into a data.frame
#' @param input_vcf input VCF file to parse
#' @param popfreq_field (default "MAF") INFO field containing the population frequency of interest
#' @return data.frame with columns Chr, Pos, Ref, Alt, and Pop_Freq.  If popfreq_field is missing for a
#' variant, its population frequency is assumed to be zero.
#' @export
vcf_to_popfreq_df <- function(input_vcf, popfreq_field = "MAF") {
    X <- vcf_to_df(input_vcf)
    if (!("INFO" %in% colnames(X))) {
        w <- which(is.na(colnames(X)))
        if (length(w)==1) colnames(X)[w] <- "INFO"
    }
    popfreq_field_equals <- paste0(popfreq_field, "=")
    s1 <- lapply(strsplit(X$INFO,";"), function(x) {
        gsub(popfreq_field_equals, "", x[grep(popfreq_field_equals,x)])
    })
    s1 <- lapply(s1,function(x) ifelse(length(x) == 0,"0",x))
    X$Pop_Freq <- as.numeric(unlist(s1))
    X <- X[,c("#CHROM","POS","REF","ALT","Pop_Freq")]
    colnames(X)[1:4] <- c("Chr","Pos","Ref","Alt")
    rownames(X) <- NULL
    return(X)
}

#' restrict_to_het
#'
#' Restrict a zygosity data.frame to SNPs that are inferred to be heterozygous
#'
#' @param zygosity_data data.frame containing variants
#' @param mincov Minimum total coverage to allow in order to keep a variant (default 30)
#' @param filtertype "normal" or "tumor" - which sample to use for filtering, default "normal"
#' @param numproc Number of processors to use, default 1
#' @return data.frame containing zygosity data for variants that are inferred to be heterozygous
#'
#' @export
restrict_to_het <- function(zygosity_data, mincov = 30, filtertype = "normal",
                            numproc = 1, gender_inference_file = NULL) {
    suppressPackageStartupMessages(require(DNAcopy))
    suppressPackageStartupMessages(require(IRanges))
    suppressPackageStartupMessages(require(parallel))
    suppressPackageStartupMessages(require(gtools))
    ## remove variants based on read counts
    if (filtertype == "normal") {
        zygosity_data <- subset(zygosity_data, Tot.N >= mincov)
        zygosity_data <- subset(zygosity_data, Var.N >= 10 & Ref.N >= 10)
    }
    if (filtertype == "tumor") {
        zygosity_data <- subset(zygosity_data, Tot.T >= mincov & Var.T >= 5 & Ref.T >= 2 & Var.T/Tot.T <= 0.99)
        zygosity_data <- subset(zygosity_data, Pop_Freq > 0)
    }
    if (nrow(zygosity_data) > 0) {
        ## multiple variants at the same position
        chrpos <- paste(zygosity_data$Chr,zygosity_data$Pos)
        zygosity_data <- subset(zygosity_data,!(chrpos %in% chrpos[which(duplicated(chrpos))]))
        ## eliminate points with vf = p if no corresponding point with vf = (2*m2-p) within +/- 10 points
        if (filtertype == "tumor") {
            m2 <- mean(zygosity_data$Var.T/zygosity_data$Tot.T, na.rm = TRUE)
            m0 <- mclapply(1:nrow(zygosity_data),function(i) {
                inds <- setdiff(intersect((i-10):(i+10),1:nrow(zygosity_data)),i)
                b1 <- zygosity_data$Var.T[i]/zygosity_data$Tot.T[i]
                b2 <- 2*m2-b1
                mm <- min(abs(b2-zygosity_data$Var.T[inds]/zygosity_data$Tot.T[inds]))
                if (mm > 0.05) {
                    w1 <- which(zygosity_data$Chr==zygosity_data$Chr[i] &
                                zygosity_data$Pos==zygosity_data$Pos[i])
                } else {
                    w1 <- NA
                }
                return(w1)
            }, mc.cores = numproc)
            m0 <- setdiff(unlist(m0),NA)
            if (nrow(zygosity_data) > 0) zygosity_data <- zygosity_data[setdiff(1:nrow(zygosity_data), m0),]
        }
    }
    ## segment normal variant fractions and remove those deviating too much from from 50%
    if (filtertype == "normal") {
        regions.to.remove <- 5000
        if (nrow(zygosity_data) > 0) {
            while (regions.to.remove > 0) {
                C1 <- CNA(genomdat = abs(zygosity_data$Var.N/zygosity_data$Tot.N - 0.5),
                          chrom = zygosity_data$Chr, maploc = zygosity_data$Pos, data.type = "logratio")
                set.seed(1000)
                S1 <- segment(C1)
                to.remove <- subset(S1$out, seg.mean >= 0.1)
                regions.to.remove <- nrow(to.remove)
                for (u in unique(to.remove$chrom)) {
                    wZ <- which(zygosity_data$Chr==u)
                    wT <- which(to.remove$chrom==u)
                    targets.IRanges <- IRanges(start=to.remove$loc.start[wT],end=to.remove$loc.end[wT])
                    snps.IRanges <- IRanges(start=zygosity_data$Pos[wZ],end=zygosity_data$Pos[wZ])
                    overlaps <- findOverlaps(query=snps.IRanges,subject=targets.IRanges)
                    w.new <- wZ[queryHits(overlaps)]
                    if (nrow(zygosity_data) > 0) {
                        zygosity_data <- zygosity_data[setdiff(1:nrow(zygosity_data), w.new), ]
                    }
                }
            }
        }
    }
    ## remove outliers
    if (nrow(zygosity_data) > 0) {
        if (filtertype == "normal") global_mean <- mean(zygosity_data$Var.N/zygosity_data$Tot.N, na.rm = TRUE)
        if (filtertype == "tumor") global_mean <- mean(zygosity_data$Var.T/zygosity_data$Tot.T, na.rm = TRUE)
        if (filtertype == "normal") {
            zygosity_data <- subset(zygosity_data, abs(zygosity_data$Var.N/zygosity_data$Tot.N - global_mean) < 0.1)
        }
        if (filtertype == "tumor") {
            r1 <- running(abs(zygosity_data$Var.T/zygosity_data$Tot.T - global_mean),
                          width = 5, fun = mean, align = "center", allow.fewer = TRUE)
            r1 <- r1[setdiff(1:length(r1),c(1:2,(length(r1)-1):length(r1)))]
            zygosity_data <- subset(zygosity_data,
                                    abs(r1 - abs(zygosity_data$Var.T/zygosity_data$Tot.T - global_mean)) < 0.1)
        }
    }
    ## remove variants that are too close to each other
    if (nrow(zygosity_data) > 0) {
        d1 <- diff(zygosity_data$Pos)
        w <- which(d1 > 0 & d1 < 100)
        w2 <- unique(c(w,w+1))
        zygosity_data <- zygosity_data[setdiff(1:nrow(zygosity_data), w2),]
    }
    ## no heterozygous snps on Y
    zygosity_data <- subset(zygosity_data, Chr != "Y")
    ## no heterozygous snps on X if male, except for pseudo-autosomal regions
    if (!is.null(gender_inference_file)) {
        g1 <- readLines(gender_inference_file)
        if (filtertype == "tumor") {
            if (g1[1]=="male") {
                data(PAR_coords_hg19)
                w.ok <- which(zygosity_data$Chr != "X")
                for (i in 1:nrow(PAR_coords_hg19)) {
                    w.ok <- c(w.ok,
                              which(zygosity_data$Chr == "X" &
                                    zygosity_data$Pos >= PAR_coords_hg19$Start[i] &
                                    zygosity_data$Pos <= PAR_coords_hg19$End[i]))
                }
                zygosity_data <- zygosity_data[w.ok,]
            }
        }
        if (filtertype == "normal") {
            if (g1[2]=="male") {
                data(PAR_coords_hg19)
                w.ok <- which(zygosity_data$Chr != "X")
                for (i in 1:nrow(PAR_coords_hg19)) {
                    w.ok <- c(w.ok,
                              which(zygosity_data$Chr == "X" &
                                    zygosity_data$Pos >= PAR_coords_hg19$Start[i] &
                                    zygosity_data$Pos <= PAR_coords_hg19$End[i]))
                }
                zygosity_data <- zygosity_data[w.ok,]
            }
        }
    }
    return(zygosity_data)
}

