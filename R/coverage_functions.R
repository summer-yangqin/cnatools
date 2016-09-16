#' get_chrinfo
#'
#' @param chrnames Vector of chromosome names to use (default 1:22,X,Y)
#' @param numproc Number of processors to use
#' @return data.frame containing ChrName, Length_Orig (in base pairs), Length, and PrevLengths
#'
#' @export
get_chrinfo <- function(chrnames, numproc=2) {
    suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg19))
    options(scipen=99)
    ChrInfo <- data.frame(ChrName=chrnames,stringsAsFactors=F)
    ChrInfo$Length_Orig <- unlist(mclapply(as.list(ChrInfo$ChrName),function(chrname) {
        return(length(Hsapiens[[paste0("chr",chrname)]]))
    },mc.cores=numproc))
    ChrInfo$Length <- ChrInfo$Length_Orig
    ChrInfo$Length <- ChrInfo$Length+20000000
    ChrInfo$Length <- ChrInfo$Length/10000
    ChrInfo$PrevLengths <- c(0,cumsum(ChrInfo$Length[-nrow(ChrInfo)]))
    return(ChrInfo)
}

#' get_max_readlength
#' @param bamfile path to bam file
#' @param N number of reads to use
#' @return numeric value representing maximum read length (NULL if bamfile doesn't exist)
#'
#' @export
get_max_readlength <- function(bamfile, N = 10000) {
    y <- NULL
    if (file.exists(bamfile)) {
        x <- system(paste0("samtools view ",bamfile," | head -n ",N),intern=T)
        y <- max(unlist(lapply(strsplit(x,"\t"),function(a) nchar(a[10]))))
    }
    return(y)
}

#' bed2df
#' @param bedfile path to bedfile (can contain optional header, ignored if it exists.)
#' @param add_region_names if TRUE, adds 1-based region names in the form chr:start-end
#' @return data.frame containing columns Chr, Start, End, and optionally Region for targets on 1-22,X,Y
#'
#' @export
bed2df <- function(bedfile, add_region_names = TRUE) {
    b1 <- NULL
    if (file.exists(bedfile)) {
        b1 <- readLines(bedfile)
        b1 <- b1[which(substr(b1,1,1) != "#")]
        b1 <- do.call(rbind,lapply(strsplit(b1,"\t"),function(a) a[1:3]))
        b1 <- as.data.frame(b1,stringsAsFactors=F)
        colnames(b1) <- c("Chr","Start","End")
        b1$Start <- as.numeric(b1$Start)
        b1$End <- as.numeric(b1$End)
        b1$Chr <- gsub("chr","",b1$Chr)
        b1 <- subset(b1,Chr %in% c(1:22,"X","Y"))
        o1 <- order(factor(b1$Chr,levels=c(1:22,"X","Y")),
                    b1$Start,b1$End)
        b1 <- b1[o1,]
        if (add_region_names) b1$Region <- paste0(b1$Chr,":",b1$Start+1,"-",b1$End)
        rownames(b1) <- NULL
    }
    return(b1)
}

#' prepare_subread_annotation
#' @param bedfile path to bedfile (can contain optional header, ignored if it exists.)
#' @param padding number of bases (default 0) on each side of target on which to compute coverage
#' @return data.frame compatible with featureCounts function
#' 
#' @export
prepare_subread_annotation <- function(bedfile, padding = 0) {
    R1 <- bed2df(bedfile, add_region_names = FALSE)
    R1$GeneID <- 1:nrow(R1)
    R1$Strand <- "+"
    R1 <- R1[,c("GeneID","Chr","Start","End","Strand")]
    R1$Start <- R1$Start-padding
    R1$End <- R1$End+padding
    return(R1)
}

#' getFeatureCounts
#'
#' @param bamfile full path to bam file
#' @param bedfile full path to bed file
#' @param padding number of bases (default 0) on each side of target on which to quantify reads
#' @param threads number of threads to use (default 1)
#' @return data.frame containing per-target counts from the bam file for each target in the bed file
#' 
#' @export
getFeatureCounts <- function(bamfile, bedfile, padding=0, threads=1, writedir=NULL) {
    suppressPackageStartupMessages(require(Rsubread))
    d1 <- getwd()
    t1 <- tempfile()
    setwd(dirname(t1))
    R1 <- prepare_subread_annotation(bedfile, padding)
    F1 <- featureCounts(files=bamfile,
                        annot.ext=R1,
                        isPairedEnd=FALSE,
                        allowMultiOverlap=TRUE,
                        ignoreDup=TRUE,
                        nthreads=threads,
                        reportReads=FALSE)
    F2 <- F1$annotation; F2$Counts <- F1$counts[,1]
    setwd(d1)
    return(F2)
}

#' multiComputeCoverage
#'
#' @param bamlist list containing paths to bam files.  may have length 1 or 2.
#' @param bedfile bed file on which to compute coverage
#' @param padding number of bases (default 0) by which to pad targets in bed file
#' @param numproc number of processors to use (default 2)
#' @return data.frame with columns Region, Chr, Start, End, Midpoint, Length, and either {Coverage} or {Coverage.T,Coverage.N} depending on length of bamlist
#'
#' @export
multiComputeCoverage <- function(bamlist, bedfile, padding = 0, numproc = 2) {
    suppressPackageStartupMessages(require(parallel))
    number_bams <- length(bamlist)
    readLengths <- unlist(lapply(bamlist, get_max_readlength))
    M1 <- mclapply(bamlist, getFeatureCounts, bedfile = bedfile,
                   padding = padding, threads = floor(numproc/number_bams),
                   mc.cores = numproc)
    if (number_bams == 2) {
        colnames(M1[[1]])[7] <- "Counts_T"; colnames(M1[[2]])[7] <- "Counts_N"
        X <- merge(M1[[1]],M1[[2]],by=c("GeneID","Chr","Start","End","Strand","Length"),sort=FALSE)
        X$MeanCovT <- round(X$Counts_T*readLengths[1]/X$Length,2)
        X$MeanCovN <- round(X$Counts_N*readLengths[2]/X$Length,2)
        covcolnames <- c("MeanCovT","MeanCovN")
    }
    if (number_bams == 1) {
        colnames(M1[[1]])[7] <- "Counts_T"
        X <- M1[[1]]
        X$MeanCovT <- round(X$Counts_T*readLengths[1]/X$Length,2)
        covcolnames <- "MeanCovT"
    }
    ## bedfiles are 0-based, this should be 1-based
    X$Region <- paste0(X$Chr,":",X$Start+1,"-",X$End)
    X$Midpoint <- round((X$Start+X$End)/2,0)
    X$ChrNum <- as.numeric(factor(X$Chr,levels=c(1:22,"X","Y")))
    X <- X[order(X$ChrNum,X$Start,X$End),]
    X <- X[,c("Region","Chr","Start","End","Midpoint","Length",covcolnames)]
    return(X)
}

#' get_gc
#'
#' @param startvec vector of start positions
#' @param lengthvec vector of exon lengths
#' @param chrom chromosome name (from 1-22,X,Y)
#' @param numproc number of processors (default 2) to use
#' @return numeric vector containing mean GC contents for the regions specified
#'
#' @export
get_gc <- function(startvec, lengthvec, chrom, numproc = 2) {
	suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg19))
        suppressPackageStartupMessages(require(parallel))
	chrseq1 <- Hsapiens[[paste0("chr",chrom)]]
	get_gc_1 <- function(startpos, numbases) {
		a <- alphabetFrequency(DNAString(chrseq1, start = startpos, nchar = numbases))
		return((a["C"]+a["G"])/sum(a))
	}
	m1 <- mcmapply(get_gc_1, startvec, lengthvec, mc.cores = numproc)
	return(m1)
}

#' add_gc_df
#'
#' @param data_frame_input input data.frame, with column names including Start and Length
#' @param numproc number of processors to use (default 2)
#' @return data.frame containing the input columns plus an extra column 'GC' representing mean GC content
#'
#' @export
add_gc_df <- function(data_frame_input, numproc = 2) {
    data_frame_output <- data_frame_input
    data_frame_output$GC <- NA
    for (chr in unique(data_frame_output$Chr)) {
        w <- which(data_frame_output$Chr == chr)
        data_frame_output$GC[w] <- get_gc(startvec = data_frame_output$Start[w],
                                          lengthvec = data_frame_output$Length[w],
                                          chrom = chr, numproc = numproc)
    }
    data_frame_output$GC <- round(data_frame_output$GC, 4)
    return(data_frame_output)
}

#' prepare_normal_pool
#'
#' Prepare a pool containing normal coverage quantifications
#' @param pool_config_file input tab-delimited text file (no header) with two columns containing paths to bam files and corresponding sample IDs, respectively
#' @param bedfile path to bed file containing target regions on which to compute coverage
#' @param outfile output file in which to save quantifications (.rds format)
#' @param numproc number of processors (default 2) to use
#' @param padding number of bases by which to extend targets on each side (default 0)
#' @return data.frame containing pool data (targets x samples), which is also saved to outfile
#'
#' @export
prepare_normal_pool <- function(pool_config_file, bedfile, outfile, numproc = 2, padding = 0) {
    suppressPackageStartupMessages(require(parallel))
    options(scipen = 99)
    bam.info <- read.delim(pool_config_file, stringsAsFactors = F, header = F)
    bam.info <- subset(bam.info,file.exists(bam.info[,1]))
    readLengths <- unlist(lapply(as.list(bam.info[,1]),get_max_readlength))
    M1 <- mclapply(as.list(bam.info[,1]), getFeatureCounts, bedfile = BEDFILE,
                   padding = padding, threads = 1,
                   mc.cores = numproc)
    for (i in 1:length(M1)) {
        M1[[i]]$Read.Length <- readLengths[i]
        M1[[i]]$Sample.ID <- bam.info[i,2]
    }
    M1 <- lapply(M1, function(p) {
        p$Coverage <- round(p$Counts * p$Read.Length/p$Length, 2)
        colnames(p)[which(colnames(p)=="Coverage")] <- unique(p$Sample.ID)
        p$GeneID <- p$Strand <- p$Length <- p$Counts <- p$Sample.ID <- p$Read.Length <- NULL
        return(p)
    })
    for (i in 1:length(M1)) {
        if (i==1) X <- M1[[i]] else X <- merge(X,M1[[i]],all = TRUE)
    }
    Y <- bed2df(bedfile, add_region_names = TRUE); Y$Ind <- 1:nrow(Y)
    Z <- merge(Y,X,by = c("Chr","Start","End"),all.x = TRUE,all.y = FALSE)
    Z <- Z[order(Z$Ind),]
    Z <- Z[,c("Chr","Start","End",bam.info[,2])]
    rownames(Z) <- NULL
    pool_data <- Z
    saveRDS(pool_data, file = outfile)
    return(pool_data)
}

#' infer_gender
#'
#' Infers gender based on Y-chromosome coverage
#'
#' @param coverage_vec vector containing coverage values across all target regions
#' @param chrom_vec vector of chromosomes (from 1-22,X,Y)
#' @return character vector containing inferred gender ("female" or "male")
#'
#' @export
infer_gender <- function(coverage_vec, chrom_vec) {
    cov_by_chrom <- tapply(coverage_vec, chrom_vec, mean)
    cov_pct_y <- cov_by_chrom["Y"]/sum(cov_by_chrom)
    gender <- ifelse(cov_pct_y < 0.005,"female","male")
    names(gender) <- NULL
    return(gender)
}

#' get_quant_mx
#'
#' Quantify a sample against a matrix of reference samples
#'
#' @param input_vec vector containing coverage quantifications for the sample of interest
#' @param ref_mx matrix (rows = targets, columns = samples) containing quantifications for the reference samples
#' @param gc_vec vector of mean GC content for the targets
#' @param mincov minimum reference coverage to use for quantification
#' @return matrix of corrected values (in the same order as the input reference matrix)
#'
#' @export
get_quant_mx <- function(input_vec, ref_mx, gc_vec, mincov = 30) {
    Z <- cbind(Sample = input_vec, ref_mx)
    samplenames <- colnames(ref_mx)
    ## normalize by average coverage
    a1 <- apply(Z, 2, function(a) sum(a/1000))
    m1 <- median(a1)
    for (col in 1:ncol(Z)) {
        obs1 <- a1[col]
        Z[,col] <- round(Z[,col]*m1/obs1,1)
    }
    corrected_mx <- apply(Z, 2, function(y) {
        lr.raw <- log2((Z[,1]+1)/(y+1))
        w <- which(!is.na(lr.raw))
        L <- lowess(x = gc_vec[w], y = lr.raw[w], f = 0.05)
        lr.corrected <- rep(NA, length(lr.raw))
        lr.corrected[w] <- lr.raw[w] - L$y[match(gc_vec[w], L$x)]
        lr.corrected <- lr.corrected - median(lr.corrected,na.rm = T)
        return(lr.corrected)
    })
    ## impute NA for low-coverage values
    corrected_mx <- corrected_mx[,setdiff(colnames(corrected_mx),"Sample")]
    corrected_mx <- ifelse(ref_mx < mincov, NA, corrected_mx)
    return(corrected_mx)
}

#' prepare_coverage_data
#'
#' Prepare coverage data from bam files, a bed file, and a file containing quantifications across a normal pool
#'
#' @param coverage_df data.frame of mean coverage quantifications (required)
#' @param normalpool_file filepath containing normal pool constructed with prepare_normal_pool function (optional)
#' @param pool_matched_gender if TRUE, pool has already been matched to the sample of interest.  default FALSE.
#' @param add_plotting_coords if TRUE, adds an extra column (Plot_Coord) containing coordinates for plotting.
#' @param numproc number of processors to use (default 2)
#' @return data.frame containing GC-corrected coverage data of tumor (and optionally normal) sample against the normal sample and/or the normal pool
#'
#' @export 
prepare_coverage_data <- function(coverage_df, normalpool_file = NULL, pool_matched_gender = FALSE,
                                  add_plotting_coords = TRUE, polymorphic_regionsfile = NULL,
                                  numproc = 2, mincov = 30, gender = "infer", gender_inference_file = NULL) {
    suppressPackageStartupMessages(require(IRanges))
    if ("MeanCovN" %in% colnames(coverage_df)) normal_sample <- TRUE else normal_sample <- FALSE
    gender_pred_T <- gender_pred_N <- NA
    if (normal_sample) {
        if (gender == "infer") {
            gender_pred_N <- infer_gender(coverage_vec = coverage_df$MeanCovN,
                                          chrom_vec = coverage_df$Chr)
        } else {
            gender_pred_N <- gender
        }
    }
    ## add GC content
    cat("Adding GC content...\n")
    coverage_df <- add_gc_df(coverage_df, numproc = numproc)
    ## if normal pool, predict gender and select most similar samples.
    if (!is.null(normalpool_file)) {
        cat("Selecting most similar samples and quantifying against normal pool...\n")
        normal_pool_df <- readRDS(normalpool_file)
        ## put this in the same order as coverage_df
        m <- match(paste(coverage_df$Chr,coverage_df$Start,coverage_df$End),
                   paste(normal_pool_df$Chr,normal_pool_df$Start,normal_pool_df$End))
        normal_pool_df <- normal_pool_df[m,]
        if (gender == "infer") {
            gender_pred_T <- infer_gender(coverage_vec = coverage_df$MeanCovT,
                                          chrom_vec = coverage_df$Chr)
        } else {
            gender_pred_T <- gender
        }
        pool_samples <- setdiff(colnames(normal_pool_df), c("Chr","Start","End"))
        if (pool_matched_gender) {
            pool_samples_same_gender <- pool_samples
        } else {
            gender_preds_pool <- apply(normal_pool_df[,pool_samples], 2, infer_gender,
                                       chrom_vec = normal_pool_df$Chr)
            pool_samples_same_gender <- names(gender_preds_pool)[which(gender_preds_pool ==
                                                                       gender_pred_T)]
        }
        f1 <- get_quant_mx(input_vec = coverage_df$MeanCovT,
                           ref_mx = normal_pool_df[,pool_samples_same_gender],
                           gc_vec = coverage_df$GC,
                           mincov = mincov)
        similarity_score <- apply(f1, 2, function(x) sd(diff(x), na.rm = TRUE))
        most_similar_samples <- names(sort(similarity_score, decreasing = FALSE))[1:10]
        final_quants <- apply(f1[,most_similar_samples], 1, median, na.rm = TRUE)
        inds_low_coverage <- which(apply(is.na(f1),1,mean) >= 0.05)
        final_quants[inds_low_coverage] <- NA
        coverage_df$LR.Corrected.TP <- final_quants
        ## normal vs. pool quantification
        if (normal_sample) {
            if (pool_matched_gender) {
                pool_samples_same_gender <- pool_samples
            } else {
                gender_preds_pool <- apply(normal_pool_df[,pool_samples], 2, infer_gender,
                                           chrom_vec = normal_pool_df$Chr)
                pool_samples_same_gender <- names(gender_preds_pool)[which(gender_preds_pool == gender_pred_N)]
            }
            f1 <- get_quant_mx(input_vec = coverage_df$MeanCovN,
                               ref_mx = normal_pool_df[,pool_samples_same_gender],
                               gc_vec = coverage_df$GC,
                               mincov = mincov)
            similarity_score <- apply(f1, 2, function(x) sd(diff(x), na.rm = TRUE))
            most_similar_samples <- names(sort(similarity_score, decreasing = FALSE))[1:10]
            final_quants <- apply(f1[,most_similar_samples], 1, median, na.rm = TRUE)
            inds_low_coverage <- which(apply(is.na(f1),1,mean) >= 0.05)
            final_quants[inds_low_coverage] <- NA
            coverage_df$LR.Corrected.NP <- final_quants
        }
    } else {
        cat("Normal pool unspecified, skipping normal pool quantification...\n")
    }
    ## tumor vs. normal quantification
    if (normal_sample) {
        mx1 <- as.matrix(coverage_df$MeanCovN, ncol = 1); colnames(mx1) <- "Normal"
        final_quants <- get_quant_mx(input_vec = coverage_df$MeanCovT,
                                     ref_mx = mx1,
                                     gc_vec = coverage_df$GC,
                                     mincov = mincov)[,1]
        coverage_df$LR.Corrected.TN <- final_quants
    }
    ## add plotting coordinates
    if (add_plotting_coords) {
       coverage_df$Plot_Coord  <- get_plotting_coords(position_vec = coverage_df$Midpoint,
                                                      chr_vec = coverage_df$Chr, numproc = numproc)
    }
    if (!is.null(polymorphic_regionsfile)) {
        poly_regions <- read.delim(polymorphic_regionsfile, stringsAsFactors=F)
        inds_remove <- NULL
        for (chrname in c(1:22,"X","Y")) {
            wC <- which(coverage_df$Chr==chrname)
            wG <- which(poly_regions$Chr==chrname)
            IC <- IRanges(start=coverage_df$Start[wC],end=coverage_df$End[wC])
            IG <- IRanges(start=poly_regions$Start[wG],poly_regions$End[wG])
            F1 <- findOverlaps(query=IG,subject=IC)
            inds_remove <- c(inds_remove,wC[queryHits(F1)])
        }
        coverage_df$Polymorphic <- rep(0,nrow(coverage_df))
        coverage_df$Polymorphic[inds_remove] <- 1
        coverage_df_poly <- coverage_df[inds_remove,]
        for (cname in c("LR.Corrected.TP","LR.Corrected.TN","LR.Corrected.NP")) {
            if (cname %in% colnames(coverage_df)) {
                coverage_df[which(coverage_df$Polymorphic==1),cname] <- NA
            }
        }
    } else {
        coverage_df_poly <- NULL
    }
    if (!is.null(gender_inference_file)) cat(c(gender_pred_T, gender_pred_N),
                                              sep = "\n", file = gender_inference_file)
    return(list(coverage_df=coverage_df, coverage_df_poly=coverage_df_poly))
}
