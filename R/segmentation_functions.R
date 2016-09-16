#' merge_short_segments
#'
#' Merge short segments of a segmentation data.frame: if a region has one non-missing exon and no snps, then merge into closest region
#' @param segdata input segmentation data.frame
#' @return data.frame containing merged segments
#' @export
merge_short_segments <- function(segdata) {
    w <- which(segdata$Num.Exons.NM==1 & segdata$Num.SNPs==0)
    for (i in w) {
        candidate.segdata <- subset(rbind(segdata[i-1,],segdata[i+1,]),Chr==segdata$Chr[i])
        m1 <- NULL
        for (j in 1:nrow(candidate.segdata)) {
            m1[j] <- min(min(abs(c(candidate.segdata$Start[j],candidate.segdata$End[j])-segdata$Start[i])),
                         min(abs(c(candidate.segdata$Start[j],candidate.segdata$End[j])-segdata$End[i])))
        }    
        segment.to.merge <- candidate.segdata$Segment[which.min(m1)]
        ind.merge <- which(segdata$Segment==segment.to.merge)
        segdata$Start[ind.merge] <- min(segdata$Start[i],segdata$Start[ind.merge])
        segdata$End[ind.merge] <- max(segdata$End[i],segdata$End[ind.merge])
        segdata$Num.Exons[ind.merge] <- segdata$Num.Exons[ind.merge]+segdata$Num.Exons[i]
        segdata$Num.Exons.NM[ind.merge] <- segdata$Num.Exons.NM[ind.merge]+segdata$Num.Exons.NM[i]
        segdata$Log2.Coverage.Ratio[ind.merge] <- (segdata$Log2.Coverage.Ratio[ind.merge]*segdata$Num.Exons.NM[ind.merge] + 
                                              segdata$Log2.Coverage.Ratio[i]*segdata$Num.Exons.NM[i])/
            (segdata$Num.Exons.NM[i]+segdata$Num.Exons.NM[ind.merge])
    }
    segdata.new <- segdata[setdiff(1:nrow(segdata),w),]
    return(segdata.new)
}

#' integrated_segmentation
#'
#' Joint segmentation of coverage and zygosity data.
#'
#' @param coverage_data data.frame containing coverage quantifications
#' @param zygosity_data data.frame containing variants
#' @param coverage_colname column name of coverage_data to use as log2 coverage ratios
#' @param use_normal_sample if TRUE, uses normal sample; if FALSE uses only tumor sample
#' @param exceptional_regions set of non-overlapping regions, each of which is its own segment.  data.frame with columns Chr, Start, End.  default = NULL.
#' @return data.frame of segments with columns Segment, Chr, Start, End, Num.Agg.SNPs, Dev.Zyg.Seg, Num.Exons, Log2.Coverage.Ratio, Num.SNPs, Num.Exons.NM
#'
#' @export
integrated_segmentation <- function(coverage_data, zygosity_data, coverage_colname, use_normal_sample = TRUE, exceptional_regions = NULL) {
    suppressPackageStartupMessages(require(IRanges))
    suppressPackageStartupMessages(require(DNAcopy))
    ## use L2.CovRatio as coverage column name
    coverage_data$L2.CovRatio <- coverage_data[, coverage_colname]
    ## if zero coverage in the tumor, change to 1 in order to avoid problems with division by zero
    zygosity_data$Tot.T[which(zygosity_data$Tot.T == 0)] <- 1
    ## compute Zyg.Dev (deviation from global zygosity mean)
    global_mean <- sum(zygosity_data$Var.T, na.rm = TRUE)/sum(zygosity_data$Tot.T, na.rm = TRUE)
    zygosity_data$Zyg.Dev <- abs(zygosity_data$Var.T/zygosity_data$Tot.T - global_mean)
    ## find correspondence between snps and targets
    zygosity_data$Region <- rep(NA, nrow(zygosity_data))
    for (chrom in unique(c(zygosity_data$Chr, coverage_data$Chr))) {
        wZ <- which(zygosity_data$Chr == chrom)
        wC <- which(coverage_data$Chr == chrom)
        targets.IRanges <- IRanges(start = coverage_data$Start[wC], end = coverage_data$End[wC])
        snps.IRanges <- IRanges(start = zygosity_data$Pos[wZ], end = zygosity_data$Pos[wZ])
        overlaps <- findOverlaps(query = snps.IRanges, subject = targets.IRanges)
        zygosity_data$Region[wZ][queryHits(overlaps)] <- coverage_data$Region[wC][subjectHits(overlaps)]
    }
    ## subdivide snps into on-target and off-target snps.
    ## if multiple snps lie within the same target, select the one with max coverage in appropriate sample
    ## so that breakpoints cannot occur within a target.
    Z_offtarget <- subset(zygosity_data, is.na(Region))
    Z_ontarget <- subset(zygosity_data, !is.na(Region))
    Z_ontarget.split <- lapply(split(Z_ontarget, Z_ontarget$Region), function(p) {
        if (use_normal_sample) w <- which.max(p$Tot.N) else w <- which.max(p$Tot.T)
        p[w, ]
    })
    Z_ontarget <- do.call(rbind, Z_ontarget.split)
    Z <- rbind(Z_offtarget, Z_ontarget)
    Z <- Z[order(factor(Z$Chr, levels = c(1:22, "X", "Y")), Z$Pos), ]
    rownames(Z) <- NULL
    if (!is.null(exceptional_regions)) {
        Z$Addition <- rep(0,nrow(Z))
        if (nrow(exceptional_regions) > 0) {
            for (i in 1:nrow(exceptional_regions)) {
                wi <- which(Z$Chr == exceptional_regions$Chr[i] &
                            Z$Pos <= exceptional_regions$End[i] &
                            Z$Pos >= exceptional_regions$Start[i])
                if (length(wi) > 0) {
                    Z$Addition[wi] <- 100*i
                    Z$Zyg.Dev[wi] <- Z$Zyg.Dev[wi]+Z$Addition[wi]
                }
            }
        }
    }
    ## segment the allelic fractions
    if (nrow(Z) > 0) {
        C1 <- CNA(genomdat = Z$Zyg.Dev, chrom = Z$Chr, maploc = Z$Pos, data.type = "logratio")
        set.seed(1000)
        S1 <- segment(C1, verbose = 0)
        seg.out <- S1$out
        colnames(seg.out)[which(colnames(seg.out) == "ID")] <- "Segment"
        for (cname in c("loc.start", "loc.end", "num.mark", "seg.mean")) {
            seg.out[, cname] <- as.numeric(seg.out[, cname])
        }
    } else {
        seg.out <- data.frame(Segment = numeric(0), chrom = character(0), loc.start = numeric(0),
                              loc.end = numeric(0), num.mark = numeric(0), seg.mean = numeric(0),
                              stringsAsFactors = F)
    }
    if (!is.null(exceptional_regions)) {
        seg.out$Exclude <- rep(0,nrow(seg.out))
        if (nrow(exceptional_regions) > 0) {
            for (i in 1:nrow(exceptional_regions)) {
                wi <- which(seg.out$chrom == exceptional_regions$Chr[i] &
                            seg.out$loc.start <= exceptional_regions$End[i] &
                            seg.out$loc.end >= exceptional_regions$Start[i])
                if (length(wi) > 0) {
                    seg.out$seg.mean[wi] <- seg.out$seg.mean[wi] - 100*i
                    seg.out$loc.start[wi] <- exceptional_regions$Start[i]
                    seg.out$loc.end[wi] <- exceptional_regions$End[i]
                    if (length(wi) > 1) {
                        sm_new <- sum(seg.out$num.mark[wi]*seg.out$seg.mean[wi])/sum(seg.out$num.mark[wi])
                        nm_new <- sum(seg.out$num.mark[wi])
                        seg.out$seg.mean[wi] <- sm_new
                        seg.out$num.mark[wi] <- nm_new
                        seg.out$Exclude[setdiff(wi,wi[1])] <- 1
                    }
                }
            }
        }
        seg.out <- subset(seg.out, Exclude == 0)
        seg.out$Exclude <- NULL
        Z$Zyg.Dev <- Z$Zyg.Dev - Z$Addition
        Z$Addition <- NULL
    }
    ## extend segments to the starts and ends of chromosomes
    seg.out <- seg.out[order(seg.out$chrom, seg.out$loc.start, seg.out$loc.end), ]
    seg.out <- do.call(rbind, lapply(split(seg.out, seg.out$chrom), function(p) {
        p$loc.start[1] <- 0
        p$loc.end[nrow(p)] <- 300000001
        p
    }))
    ## if chromosomes are missing, add rows corresponding to these whole chromosomes:
    missing.chrs <- setdiff(as.character(c(1:22, "X", "Y")), seg.out$chrom)
    if (length(missing.chrs) > 0) {
        for (chr in missing.chrs) {
            seg.out <- rbind(seg.out, data.frame(Segment = "NA", chrom = chr, 
                                                 loc.start = 0, loc.end = 300000001, num.mark = 0, seg.mean = NA, 
                                                 stringsAsFactors = F))	
        }
    }
    seg.out <- seg.out[order(factor(seg.out$chrom, levels = c(1:22, "X", "Y")),
                             seg.out$loc.start, seg.out$loc.end), ]
    seg.out$Segment <- 1:nrow(seg.out)
    rownames(seg.out) <- NULL
    ## add missing segments (there could be a set of exons between consecutive snps)
    seg.out.split <- split(seg.out, seg.out$chrom)
    seg.out.split <- lapply(seg.out.split, function(p) {
        if (nrow(p) > 1) {
            for(i in 1:(nrow(p)-1)) {
                p2 <- data.frame(Segment = NA, chrom = unique(p$chrom), 
                                 loc.start = p$loc.end[i]+1, loc.end = p$loc.start[i+1]-1,
                                 num.mark = 0, seg.mean = NA, stringsAsFactors = F)
                p <- rbind(p, p2)
            }
            p <- p[order(p$loc.start), ]
        }
        return(p)
    })
    seg.out <- do.call(rbind, seg.out.split)
    seg.out <- seg.out[order(factor(seg.out$chrom, levels = c(1:22, "X", "Y")),
                             seg.out$loc.start, seg.out$loc.end), ]
    seg.out$Segment <- 1:nrow(seg.out)
    rownames(seg.out) <- NULL
    ## look for large coverage ratio changes within a segment, and break into smaller segments if necessary
    seg.refined <- list()
    for (i in 1:nrow(seg.out)) {
        seg.refined[[i]] <- data.frame(ID = character(0), chrom = character(0), loc.start = integer(0), 
                                       loc.end = integer(0), num.mark = integer(0), seg.mean = numeric(0), 
                                       stringsAsFactors = F)
        Ctemp <- subset(coverage_data, Chr == seg.out$chrom[i] & Midpoint >= seg.out$loc.start[i] &
                                       Midpoint <= seg.out$loc.end[i] & !is.na(L2.CovRatio))
        if (!is.null(exceptional_regions)) {
            Ctemp$Addition <- rep(0,nrow(Ctemp))
            if (nrow(exceptional_regions) > 0) {
                for (j in 1:nrow(exceptional_regions)) {
                    wj <- which(Ctemp$Chr == exceptional_regions$Chr[j] &
                                Ctemp$Midpoint <= exceptional_regions$End[j] &
                                Ctemp$Midpoint >= exceptional_regions$Start[j])
                    if (length(wj) > 0) {
                        Ctemp$Addition[wj] <- 100*j
                        Ctemp$L2.CovRatio[wj] <- Ctemp$L2.CovRatio[wj]+Ctemp$Addition[wj]
                    }
                }
            }
        }
        if (nrow(Ctemp) > 0) {
            C2 <- CNA(genomdat = Ctemp$L2.CovRatio, chrom = Ctemp$Chr, 
                      maploc = Ctemp$Midpoint, data.type = "logratio")
            set.seed(1000); S2 <- segment(C2, undo.splits = "sdundo", undo.SD = 2.5, min.width = 5, verbose = 0)
            seg.refined[[i]] <- S2$out
            if (!is.null(exceptional_regions)) {
                if (nrow(exceptional_regions) > 0) {
                    for (j in 1:nrow(exceptional_regions)) {
                        wj <- which(seg.refined[[i]]$chrom == exceptional_regions$Chr[j] &
                                    seg.refined[[i]]$loc.start <= exceptional_regions$End[j] &
                                    seg.refined[[i]]$loc.end >= exceptional_regions$Start[j])
                        if (length(wj) > 0) {
                            seg.refined[[i]]$seg.mean[wj] <- seg.refined[[i]]$seg.mean[wj] - 100*j
                            seg.refined[[i]]$loc.start[wj] <- exceptional_regions$Start[j]
                            seg.refined[[i]]$loc.end[wj] <- exceptional_regions$End[j]
                        }
                    }
                }
            }
        }
    }
    if (length(seg.refined) > 0) {
        wsplit <- which(unlist(lapply(seg.refined, nrow)) > 1)
        for (j in wsplit) {
            Stemp <- seg.refined[[j]]; colnames(Stemp)[1] <- "Segment"
            Stemp$Segment <- j+0.001*(1:nrow(Stemp))
            chrtemp <- Stemp$chrom
            Stemp$loc.start[which.min(Stemp$loc.start)] <- seg.out$loc.start[j]
            Stemp$loc.end[which.max(Stemp$loc.end)] <- seg.out$loc.end[j]
            ## look up the zygosity segment mean and number of aggregated snps on these new segments
            for (k in 1:nrow(Stemp)) {
                w <- which(Z$Chr == Stemp$chrom[k] & Z$Pos >= Stemp$loc.start[k] & 
                           Z$Pos <= Stemp$loc.end[k])
                Stemp$num.mark[k] <- length(w); Stemp$seg.mean[k] <- mean(Z$Zyg.Dev[w])
            }
            ## there can be snps lying between these segments - assign these to the nearest segment
            for (k in 1:nrow(Stemp)) {
                w <- which(Z$Chr == Stemp$chrom[k] & Z$Pos > Stemp$loc.end[k] & 
                           Z$Pos < Stemp$loc.start[k+1])
                for (ind in w) {
                    wmin <- which.min(abs(Z$Pos[ind]-c(Stemp$loc.end[k], Stemp$loc.start[k+1])))
                    if (wmin == 1) Stemp$num.mark[k] <- 1+Stemp$num.mark[k]
                    if (wmin == 2) Stemp$num.mark[k+1] <- 1+Stemp$num.mark[k+1]
                }
            }
            seg.out <- rbind(seg.out, Stemp)
        }
        seg.out <- subset(seg.out, !(Segment %in% wsplit))
        seg.out <- seg.out[order(seg.out$Segment), ]
        seg.out$seg.mean[which(is.na(seg.out$seg.mean))] <- NA
        seg.out$Segment <- 1:nrow(seg.out)
    }
    ## map each snp and target to these segments, and add coverage means and number of targets to the segments
    #Z$Segment <- rep(1:nrow(seg.out), seg.out$num.mark)
    Z$Segment <- rep(NA,nrow(Z))
    for (chrom in unique(Z$Chr)) {
        wZ <- which(Z$Chr==chrom)
        wS <- which(seg.out$chrom==chrom)
        IZ <- IRanges(start=Z$Pos[wZ],end=Z$Pos[wZ])
        IS <- IRanges(start=seg.out$loc.start[wS],end=seg.out$loc.end[wS])
        FZS <- findOverlaps(query=IZ,subject=IS)
        Z$Segment[wZ[queryHits(FZS)]] <- wS[subjectHits(FZS)]
    }
    coverage_data$Segment <- NA
    for (chrom in unique(c(seg.out$chrom, coverage_data$Chr))) {
        wS <- which(seg.out$chrom == chrom)
        wC <- which(coverage_data$Chr == chrom)
        segments.IRanges <- IRanges(start = seg.out$loc.start[wS], end = seg.out$loc.end[wS])
        coverage.IRanges <- IRanges(start = coverage_data$Midpoint[wC], end = coverage_data$Midpoint[wC])
        overlaps <- findOverlaps(query = coverage.IRanges, subject = segments.IRanges)
        coverage_data$Segment[wC[queryHits(overlaps)]] <- wS[subjectHits(overlaps)]
    }
    colnames(seg.out) <- c("Segment", "Chr", "Start", "End", "Num.Agg.SNPs", "Dev.Zyg.Seg")
    seg.out$Num.Exons <- NA; seg.out$Log2.Coverage.Ratio <- NA
    seg.out <- seg.out[order(seg.out$Segment), ]; rownames(seg.out) <- NULL
    seg.out$Num.SNPs <- 0
    for (chrom in unique(c(seg.out$Chr, zygosity_data$Chr))) {
        wS <- which(seg.out$Chr == chrom)
        wZ <- which(zygosity_data$Chr == chrom)
        segments.IRanges <- IRanges(start = seg.out$Start[wS], end = seg.out$End[wS])
        zyg.IRanges <- IRanges(start = zygosity_data$Pos[wZ], end = zygosity_data$Pos[wZ])
        overlaps <- findOverlaps(query = zyg.IRanges, subject = segments.IRanges)
        t1 <- tapply(queryHits(overlaps), subjectHits(overlaps), length)
        seg.out$Num.SNPs[wS[as.numeric(names(t1))]] <- t1
    }
    seg.out$Num.Exons <- rep(0, nrow(seg.out))
    seg.out$Num.Exons.NM <- rep(0, nrow(seg.out))
    seg.out$Log2.Coverage.Ratio <- rep(NA, nrow(seg.out))
    for (chrom in unique(c(seg.out$Chr, coverage_data$Chr))) {
        wS <- which(seg.out$Chr == chrom)
        wC <- which(coverage_data$Chr == chrom)
        segments.IRanges <- IRanges(start = seg.out$Start[wS], end = seg.out$End[wS])
        coverage.IRanges <- IRanges(start = coverage_data$Midpoint[wC], end = coverage_data$Midpoint[wC])
        overlaps <- findOverlaps(query = coverage.IRanges, subject = segments.IRanges)
        t1 <- tapply(queryHits(overlaps), subjectHits(overlaps), length)
        seg.out$Num.Exons[wS[as.numeric(names(t1))]] <- t1
        t2 <- tapply(coverage_data$L2.CovRatio[wC[queryHits(overlaps)]], 
                     subjectHits(overlaps), function(x) length(which(!is.na(x))))
        seg.out$Num.Exons.NM[wS[as.numeric(names(t2))]] <- t2
        t3 <- tapply(coverage_data$L2.CovRatio[wC[queryHits(overlaps)]], 
                     subjectHits(overlaps), mean, na.rm = T)
        seg.out$Log2.Coverage.Ratio[wS[as.numeric(names(t3))]] <- t3
    }
    seg.out$Log2.Coverage.Ratio <- round(seg.out$Log2.Coverage.Ratio, 4)
    ## segments with zero snps and zero non-missing exon quantifications get merged to the longest
    ## neighboring segment on the same chromosome
    w.missing <- which(seg.out$Num.SNPs == 0 & seg.out$Num.Exons.NM == 0)
    for (i in w.missing) {
        i2 <- intersect(c(i-1, i+1), which(seg.out$Chr == seg.out$Chr[i]))
        if (length(i2) > 0) {
            l2 <- seg.out$End[i2]-seg.out$Start[i2]
            ind.merge <- i2[which.max(l2)]
            seg.out$Start[ind.merge] <- min(c(seg.out$Start[ind.merge], seg.out$Start[i]))
            seg.out$End[ind.merge] <- max(c(seg.out$End[ind.merge], seg.out$End[i]))
        }
    }
    ## remove segments with no non-missing exons
    seg.out <- subset(seg.out, !(Num.Exons.NM == 0 & !(Chr %in% missing.chrs)))
    seg.out$Log2.Coverage.Ratio[which(is.na(seg.out$Log2.Coverage.Ratio))] <- NA ## NaN becomes NA
    ## find each segment with no snps; if its log-ratio is within 0.02 units of a neighbor, then merge.
    seg.out$To_Delete <- seg.out$Processed <- rep(0, nrow(seg.out))
    w0 <- which(seg.out$Num.SNPs == 0 & seg.out$Processed == 0)
    while (length(w0) > 0) {
        w <- w0[1]
        wnear <- intersect(which(seg.out$Chr == seg.out$Chr[w]), c(w-1, w+1))
        wnear2 <- intersect(wnear, which(abs(seg.out$Log2.Coverage.Ratio-seg.out$Log2.Coverage.Ratio[w]) <= 0.02))
        wnear3 <- wnear2[which.min(abs(seg.out$Log2.Coverage.Ratio[wnear2]-seg.out$Log2.Coverage.Ratio[w]))]
        wnear4 <- intersect(wnear3, which(seg.out$Num.SNPs > 0))
        if (length(wnear4) == 1) {
            seg.out$End[wnear4] <- max(c(seg.out$End[w], seg.out$End[wnear4]))
            ne <- seg.out$Num.Exons[w] + seg.out$Num.Exons[wnear4]
            ne.nm <- seg.out$Num.Exons.NM[w] + seg.out$Num.Exons.NM[wnear4]
            l2cr <- (seg.out$Log2.Coverage.Ratio[w]*seg.out$Num.Exons.NM[w] +
                     seg.out$Log2.Coverage.Ratio[wnear4]*seg.out$Num.Exons.NM[wnear4]) /
                (seg.out$Num.Exons.NM[w] + seg.out$Num.Exons.NM[wnear4])
            seg.out$Num.Exons[wnear4] <- ne
            seg.out$Num.Exons.NM[wnear4] <- ne.nm
            seg.out$Log2.Coverage.Ratio[wnear4] <- l2cr
            seg.out$To_Delete[w] <- 1
        }
        seg.out$Processed[w] <- 1
        seg.out <- subset(seg.out, To_Delete == 0)
        w0 <- which(seg.out$Num.SNPs == 0 & seg.out$Processed == 0)
    }
    seg.out$Processed <- seg.out$To_Delete <- NULL
    seg.out <- merge_short_segments(seg.out)
    ## replace segment numbers by integers
    seg.out$Segment <- 1:nrow(seg.out)
    return(seg.out)
}


#' coverage_only_segmentation
#'
#' Germline segmentation
#'
#' @param coverage_data coverage data.frame
#' @param coverage_colname column name to use for coverage
#' @param chrinfo chrinfo data.frame
#' @return data.frame containing segmentation
#' @export
coverage_only_segmentation <- function(coverage_data, coverage_colname, chrinfo) {
    suppressPackageStartupMessages(require(IRanges))
    suppressPackageStartupMessages(require(DNAcopy))
    coverage_data$L2.CovRatio <- coverage_data[, coverage_colname]
    coverage_data <- subset(coverage_data, !is.na(L2.CovRatio))
    C1 <- CNA(genomdat = coverage_data$L2.CovRatio, chrom = coverage_data$Chr,
              maploc = coverage_data$Midpoint, data.type = "logratio")
    set.seed(1000)
    S1 <- segment(C1, verbose = 0)
    seg.out <- S1$out
    for (cname in c("loc.start", "loc.end", "num.mark", "seg.mean")) {
        seg.out[, cname] <- as.numeric(seg.out[, cname])
    }
    colnames(seg.out)[which(colnames(seg.out)=="ID")] <- "Segment"
    ## extend segments to the starts and ends of chromosomes
    seg.out <- seg.out[order(seg.out$chrom, seg.out$loc.start, seg.out$loc.end), ]
    seg.out <- do.call(rbind, lapply(split(seg.out, seg.out$chrom), function(p) {
        p$loc.start[1] <- 0
        p$loc.end[nrow(p)] <- 300000001
        p
    }))
    ## if chromosomes are missing, add rows corresponding to these whole chromosomes:
    missing.chrs <- setdiff(as.character(c(1:22, "X", "Y")), seg.out$chrom)
    if (length(missing.chrs) > 0) {
        for (chr in missing.chrs) {
            seg.out <- rbind(seg.out, data.frame(Segment = "NA", chrom = chr, 
                                                 loc.start = 0, loc.end = 300000001, num.mark = 0, seg.mean = NA, 
                                                 stringsAsFactors = F))	
        }
    }
    seg.out <- seg.out[order(factor(seg.out$chrom, levels = c(1:22, "X", "Y")),
                             seg.out$loc.start, seg.out$loc.end), ]
    seg.out$Segment <- 1:nrow(seg.out)
    rownames(seg.out) <- NULL
    colnames(seg.out) <- c("Segment","Chr","Start","End","Exons","Log2_Coverage_Ratio")
    seg.out$Start[which(seg.out$Start==0)] <- 1
    w.end <- which(seg.out$End==300000001)
    seg.out$End[w.end] <- chrinfo$Length_Orig[match(seg.out$Chr[w.end],chrinfo$ChrName)]
    return(seg.out)
}
