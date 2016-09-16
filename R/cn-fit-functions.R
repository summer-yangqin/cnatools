#' get_tc_lower_limit
#'
#' Compute lower limit for tumor content estimation
#' @param coverage_df data.frame of coverage values
#' @param colname_LR column name containing log-ratio values to use
#' @param colname_cov column name containing coverage value to use as a minimum cutoff
#' @return vector containing lower limit for tumor content estimation, standard deviation, and minimum sequencing coverage
#' @export
get_tc_lower_limit <- function(coverage_df, colname_LR, colname_cov) {
    ## compute the sd for various cutoffs
    mincovseq <- seq(30,400,by=5)
    sigmas <- Ns <- NULL
    for (i in 1:length(mincovseq)) {
        xtmp <- subset(coverage_df,coverage_df[,colname_cov] >= mincovseq[i])
        sigmas[i] <- sd(diff(xtmp[,colname_LR]),na.rm=T)
        Ns[i] <- length(which(!is.na(xtmp[,colname_LR])))
    }
    ## some multiple of 1 sd should correspond to log2(1-p/2)
    lower.tc <- round(2*(1-2^(-1*sigmas[1])),2)
    if (is.na(lower.tc)) lower.tc <- 1
    lower.tc <- min(lower.tc,0.3)
    return(c(lower.tc,round(sigmas[1],3),mincovseq[1]))
}

#' initialize_adjustments
#'
#' Initialize adjustments for ploidy correction
#' @param adjustment_range vector of values to consider for ploidy adjustment
#' @param full_search if TRUE, scans all values in adjustment_range (slow), if FALSE, limits to likely candidates
#' @param seg_data segmentation data.frame (default NULL); required if full_search is FALSE
#' @param tc_lower lower limit of tumor content to consider
#' @return data.frame containing Adjustment and TC_Lower
#' @export
initialize_adjustments <- function(adjustment_range = seq(-0.6, 1, by = 0.02),
                                   full_search = FALSE,
                                   seg_data = NULL,
                                   tc_lower = 0.3) {
    if (!full_search) {
        b1 <- c(0,sort(unique(round(-seg_data$Log2.Coverage.Ratio,2))))
        b1 <- c(b1,b1+0.01); b1 <- c(b1,b1-0.01)
        adjustment_range <- sort(intersect(b1,round(adjustment_range,2)))
    }
    d1 <- data.frame(Adjustment = adjustment_range, stringsAsFactors=F)
    d1$TC_Lower <- rep(tc_lower, nrow(d1))
    return(d1)
}

#' append_covzyg_predictions
#'
#' Append predicted values to coverage and zygosity data
#' @param classifs data.frame of classifications
#' @param tc Inferred tumor content
#' @param C1 coverage data
#' @param Z1 zygosity data
#' @return a list with components coverage and zygosity
#' @export
append_covzyg_predictions <- function(classifs,tc,C1,Z1) {
    suppressPackageStartupMessages(require(IRanges))
    C1.with.class <- C1; C1.with.class$LR.Predicted <- NA; C1.with.class$Copy.Number.Predicted <- NA
    Z1.with.class <- Z1; Z1.with.class$Dev.VF.Observed <- abs(Z1.with.class$Var.T/Z1.with.class$Tot.T-0.47)
    Z1.with.class$VF.Predicted <- rep(NA,nrow(Z1.with.class))
    Z1.with.class$Dev.VF.Predicted <- rep(NA,nrow(Z1.with.class))
    Z1.with.class$Genotype.Predicted <- rep(NA,nrow(Z1.with.class))
    chroms <- unique(classifs$Chr)
    for (chrom in chroms) {
        wB1 <- which(classifs$Chr==chrom)
        wC1 <- which(C1.with.class$Chr==chrom)
        wZ1 <- which(Z1.with.class$Chr==chrom)
        IB1 <- IRanges(classifs$Start[wB1],classifs$End[wB1])
        IC1 <- IRanges(C1.with.class$Start[wC1],C1.with.class$End[wC1])
        IZ1 <- IRanges(Z1.with.class$Pos[wZ1],Z1.with.class$Pos[wZ1])
        F1 <- findOverlaps(query=IB1,subject=IC1)
        C1.with.class$Copy.Number.Predicted[wC1[subjectHits(F1)]] <- classifs$Copy.Number[wB1[queryHits(F1)]]
        F2 <- findOverlaps(query=IB1,subject=IZ1)
        Z1.with.class$Genotype.Predicted[wZ1[subjectHits(F2)]] <- classifs$Classification[wB1[queryHits(F2)]]
    }
    C1.with.class$LR.Predicted <- log2(2*(1-tc)+C1.with.class$Copy.Number.Predicted*tc)-1
    compute_predicted_varfrac <- function(p,A,B) {
        num <- 1+p*(B-1)
        denom <- 2+p*(A+B-2)
        return(num/denom)
    }
    w <- which(Z1.with.class$Genotype.Predicted=="1 Copy")
    Z1.with.class$VF.Predicted[w] <- compute_predicted_varfrac(tc,1,0)
    w <- which(Z1.with.class$Genotype.Predicted=="0 Copies")
    Z1.with.class$VF.Predicted[w] <- compute_predicted_varfrac(tc,0,0)
    w <- which(Z1.with.class$Genotype.Predicted=="2 Copies, LOH")
    Z1.with.class$VF.Predicted[w] <- compute_predicted_varfrac(tc,0,2)
    w <- which(Z1.with.class$Genotype.Predicted=="2 Copies, Normal Genotype")
    Z1.with.class$VF.Predicted[w] <- compute_predicted_varfrac(tc,1,1)
    g <- grep(", Genotype ",Z1.with.class$Genotype.Predicted)
    if (length(g) > 0) {
        for (i in 1:length(g)) {
            s1 <- unlist(strsplit(Z1.with.class$Genotype.Predicted[g[i]],", Genotype "))
            A <- length(which(unlist(strsplit(s1[2],"/"))=="A"))
            B <- length(which(unlist(strsplit(s1[2],"/"))=="B"))
            Z1.with.class$VF.Predicted[g[i]] <- compute_predicted_varfrac(tc,A,B)
        }
    }
    ## get predicted deviation
    Z1.with.class$Dev.VF.Predicted <- abs(Z1.with.class$VF.Predicted-0.5)
    Z1.with.class$Dev.VF.Predicted[which(Z1.with.class$VF.Predicted==0.5)] <- 0.04
    Z1.with.class <- subset(Z1.with.class,!is.na(Dev.VF.Predicted))
    return(list(coverage=C1.with.class,zygosity=Z1.with.class))
}

#' compute_fit
#'
#' Compute the fit of estimated zygosity and coverage values
#' @param zygosity data.frame containing zygosity data
#' @param coverage data.frame containing coverage data
#' @param coverage_colname column name to use for coverage
#' @param exceptional_regions regions to exclude from computation of fit statistics
#' @return vector of length two containing ratios of zygosity and coverage data deviations scaled to their theoretical lower limits
#' @export
compute_fit <- function(zygosity, coverage, coverage_colname, exceptional_regions = NULL) {
    wZ <- wC <- integer(0)
    if (!is.null(exceptional_regions)) {
        if (nrow(exceptional_regions) > 0) {
            for (i in 1:nrow(exceptional_regions)) {
                wZ <- c(wZ,which(zygosity$Chr == exceptional_regions$Chr[i] &
                                 zygosity$Pos >= exceptional_regions$Start[i] &
                                 zygosity$Pos <= exceptional_regions$End[i]))
                wC <- c(wC,which(coverage$Chr == exceptional_regions$Chr[i] &
                                 coverage$Pos >= exceptional_regions$Start[i] &
                                 coverage$Pos <= exceptional_regions$End[i]))
            }
        }
    }
    if (nrow(zygosity) > 0) zygosity <- zygosity[setdiff(1:nrow(zygosity), wZ),]
    if (nrow(coverage) > 0) coverage <- coverage[setdiff(1:nrow(coverage), wC),]
    ## theoretical limit of fit vs. observed fit
    devZ.limit <- sd(diff(zygosity$Dev.VF.Observed),na.rm=TRUE)/sqrt(2)
    devZ.observed <- sd(zygosity$Dev.VF.Observed-zygosity$Dev.VF.Predicted,na.rm=TRUE)
    devC.limit <- sd(diff(coverage[,coverage_colname]),na.rm=TRUE)/sqrt(2)
    devC.observed <- sd(coverage[,coverage_colname]-coverage$LR.Predicted,na.rm=TRUE)
    ## the two factors (lower is better)
    return(c(devZ.observed/devZ.limit,devC.observed/devC.limit))
}

#' fit_models
#'
#' Fit classification models across a range of ploidy adjustments
#' @param adj_df data.frame of adjustments to consider
#' @param coverage_data coverage data.frame
#' @param coverage_data_colname colname to use for coverage
#' @param seg_data segmented data
#' @param zygosity_data zygosity data.frame
#' @param numproc Number of processes to use
#' @param chrinfo data.frame containing chromosome information
#' @param exceptional_regions regions to exclude when computing fit statistics
#' @param gender gender of sample (default female)
#' @return a list, one element per ploidy adjustment, containing such information as the overall classifications, ploidy, and goodness of fit
#' @export
fit_models <- function(adj_df, coverage_data, coverage_data_colname, seg_data, zygosity_data, numproc, chrinfo,
                       exceptional_regions = NULL, gender = "female") {
    suppressPackageStartupMessages(require(parallel))
    cat("Evaluating fits for",nrow(adj_df),"adjustments\n")
    fitlist <- mclapply(1:nrow(adj_df),function(ind) {
        adjustment <- adj_df$Adjustment[ind]
        tclower <- adj_df$TC_Lower[ind]
        ## add adjustment to coverage data and segment data
        coverage_data_adjusted <- coverage_data
        coverage_data_adjusted[,coverage_data_colname] <- coverage_data_adjusted[,coverage_data_colname] + adjustment
        seg_data_adjusted <- seg_data
        seg_data_adjusted$Log2.Coverage.Ratio <- seg_data_adjusted$Log2.Coverage.Ratio + adjustment
        ## classify segments
        B1 <- overall_classify(Z.data = zygosity_data,
                               C.data = coverage_data_adjusted,
                               C.data.colname = coverage_data_colname,
                               seg.data = seg_data_adjusted,
                               chrinfo = chrinfo,
                               lower.limit.tc = tclower,
                               msg = paste0("Adjustment = ",adjustment,": "),
                               exceptional_regions = exceptional_regions,
                               gender = gender)
        ## append predictions
        A1 <- append_covzyg_predictions(classifs = B1$classifs,
                                        tc = B1$tc,
                                        C1 = coverage_data_adjusted,
                                        Z1 = zygosity_data)
        ## compute approximate ploidy
        ploidy <- sum(B1$classifs$Copy.Number*B1$classifs$Targeted.Exons,na.rm=TRUE)/sum(B1$classifs$Targeted.Exons,na.rm=TRUE)
        ## compute fit
        fit1 <- compute_fit(zygosity = A1$zygosity, coverage = A1$coverage,
                            coverage_colname = coverage_data_colname, exceptional_regions = exceptional_regions)
        return(list(B1 = B1, adjustment = adjustment, fit1 = fit1, ploidy = ploidy))
    },mc.cores=numproc)
    return(fitlist)
}

#' models2fits
#'
#' Compute summary stats concerning goodness of fit for a set of models
#' @param models List of models
#' @param adjustments vector of adjustments corresponding to the models
#' @return data.frame containing combined fits, ploidy, and adjustments for each model
#' @export
models2fits <- function(models, adjustments) {
    fits_12.list <- lapply(models,function(p) p$fit1)
    fits <- data.frame(Adjustment=adjustments,
                           Fit1=unlist(lapply(fits_12.list,function(x) x[1])),
                           Fit2=unlist(lapply(fits_12.list,function(x) x[2])),
                           stringsAsFactors=F)
    fits$Fit1 <- ifelse(is.na(fits$Fit1),1,fits$Fit1)
    fits$Fit2 <- ifelse(is.na(fits$Fit2),1,fits$Fit2)
    fits$Adjustment <- round(fits$Adjustment,2)
    fits <- fits[which(!duplicated(fits$Adjustment)),]
    fits <- fits[order(fits$Adjustment),]
    fits$Fit_Combined <- fits$Fit1*fits$Fit2
    fits$Fit_Combined_Penalized <- fits$Fit_Combined + abs(fits$Adjustment)
    fits$Ploidy <- unlist(lapply(models,function(p) p$ploidy))
    return(fits)
}

#' filter_to_plausible_adjustments
#'
#' Filter a set of possible adjustments to the most plausible ones based on the shape of the adjustment/fit function
#'
#' @param fits.df data.frame of fits
#' @return vector of plausible adjustments
#' @export
filter_to_plausible_adjustments <- function(fits.df) {
    ## consider local minima + global minimum + (adjustment=0) as possible adjustments
    w.lt.next <- which(diff(fits.df$Fit_Combined_Penalized) >= 0)
    w.lt.prev <- 1+which(diff(fits.df$Fit_Combined_Penalized) <= 0)
    w.local.minima <- intersect(w.lt.next,w.lt.prev)
    w.global.minimum <- which.min(fits.df$Fit_Combined_Penalized)
    WW <- sort(unique(c(w.local.minima,w.global.minimum,which(fits.df$Adjustment==0),
                        which.min(abs(fits.df$Ploidy-2)))))
    ## find adjustments that are within 0.10 units and whose fits differ by < 5% - in these cases choose
    ## the one with minimum value of Fit_Combined_Penalized
    M <- matrix(0,length(WW),length(WW))
    for (i in 1:length(WW)) {
        for (j in 1:length(WW)) {
            if (abs(fits.df$Adjustment[WW[i]]-fits.df$Adjustment[WW[j]]) < 0.10) {
                f1 <- sort(c(fits.df$Fit_Combined_Penalized[WW[i]],fits.df$Fit_Combined_Penalized[WW[j]]))
                if (f1[2]/f1[1] < 1.05) {
                    M[i,j] <- 1
                }
            }
        }
    }
    clusters <- list()
    for (j in 1:ncol(M)) {
        for (i in which(M[,j]==1)) {
            if (i %in% unlist(clusters) | j %in% unlist(clusters)) {
                w <- which(unlist(lapply(clusters,function(a) (i %in% a) | (j %in% a))))
                if (length(w) > 1) {
                    clusters[[w[1]]] <- unlist(clusters[w])
                    clusters[setdiff(w,w[1])] <- NA
                }
                clusters[[w[1]]] <- c(clusters[[w[1]]],j)
            } else {
                clusters[[length(clusters)+1]] <- j
            }
        }
    }
    clusters <- lapply(clusters,unique)
    fits.df_clustered <- lapply(clusters, function(x) fits.df[WW[x],])
    adjustments <- unlist(lapply(fits.df_clustered,function(p) p$Adjustment[which.min(p$Fit_Combined_Penalized)]))
    return(adjustments)
}

#' get_classifications_from_models
#'
#' retrieve classifications from a set of models
#' @param models list of models
#' @param adj_df data.frame of adjustments
#' @param cov_df data.frame containing unadjusted coverage values
#' @param seg_df data.frame containing unadjusted segmentation values
#' @param colname column name to use
#' @param numproc number of processors to use
#' @return list of classifications, one per model
#' @export
get_classifications_from_models <- function(models, adj_df, cov_df, seg_df, colname, numproc) {
    suppressPackageStartupMessages(require(parallel))
    C1 <- mclapply(1:nrow(adj_df),
                   function(k) {
                       adjustment.selected <- adj_df$Adjustment[k]
                       cov_df_tmp <- cov_df; seg_df_tmp <- seg_df
                       cov_df_tmp[,colname] <- cov_df_tmp[,colname] + adjustment.selected
                       seg_df_tmp$Log2.Coverage.Ratio <- seg_df_tmp$Log2.Coverage.Ratio + adjustment.selected
                       w <- which(round(unlist(lapply(models,function(p) p$adjustment)),2) ==
                                  round(adjustment.selected,2))
                       B1TMP <- models[[w]]$B1
                       return(B1TMP)
                   },mc.cores=numproc)
    return(C1)
}

get_summary_stats_from_classifications <- function(classification.list) {
    ## summary stats for homozygous deletions (proportion of genome, max segment)
    del_list <- lapply(classification.list,function(B1TMP) {
        x <- y <- 0
        dels0 <- subset(B1TMP$classifs,Copy.Number==0)
        if (nrow(dels0) > 0) {
            x <- sum(dels0$Targeted.Exons,na.rm=T)/sum(B1TMP$classifs$Targeted.Exons,na.rm=T)
            y <- max(dels0$Targeted.Exons,na.rm=T)
        }
        return(list(x=x,y=y))
    })
    proportion_genome_zerocopy <- unlist(lapply(del_list, function(a) a$x))
    max_exons_zerocopy <- unlist(lapply(del_list, function(a) a$y))
    ## summary stat for proportion of genome altered
    alt_list <- lapply(classification.list,function(B1TMP) {
        x <- 0
        alt0 <- subset(B1TMP$classifs,GainLoss != "Neutral")
        if (nrow(alt0) > 0) {
            x <- sum(alt0$Targeted.Exons,na.rm=T)/sum(B1TMP$classifs$Targeted.Exons,na.rm=T)
        }
        return(x)
    })
    proportion_genome_altered <- unlist(alt_list)
    ## estimated purity and ploidy
    estimated_purity <- unlist(lapply(classification.list, function(l1) l1$tc))
    estimated_ploidy <- unlist(lapply(classification.list, function(l1) {
        p <- l1$classif
        sum(p$Copy.Number * p$Targeted.Exons, na.rm=TRUE)/sum(p$Targeted.Exons, na.rm=TRUE)
    }))
    ## average zygosity shift for classified-neutral segments
    mean_zyg_shift_neutral <- unlist(lapply(classification.list,function(B1TMP) {
        stmp <- subset(B1TMP$classifs,GainLoss=="Neutral" & Probability >= 0.9 & Heterozygous.SNPs >= 10)
        if (nrow(stmp) > 0) a <- sum(stmp$Mean.Zygosity.Deviation*stmp$Heterozygous.SNPs,na.rm=T)/
            sum(stmp$Heterozygous.SNPs,na.rm=T) else a <- 0
        return(a)
    }))
    ## penalty for classified-1-copy segments that show no zygosity shift
    mean_zyg_shift_1copy <- unlist(lapply(classification.list,function(B1TMP) {
        stmp <- subset(B1TMP$classifs,Copy.Number==1 & Probability >= 0.9 & Heterozygous.SNPs >= 10)
        if (nrow(stmp) > 0) a <- sum(stmp$Mean.Zygosity.Deviation*stmp$Heterozygous.SNPs,na.rm=T)/
            sum(stmp$Heterozygous.SNPs,na.rm=T) else a <- NA
        return(a)
    }))
    out.df <- data.frame(Proportion_Homozygous = proportion_genome_zerocopy,
                         Max_Segment_Homozygous = max_exons_zerocopy,
                         Proportion_Altered = proportion_genome_altered,
                         Estimated_Purity = estimated_purity,
                         Estimated_Ploidy = estimated_ploidy,
                         Mean_Zygosity_Shift_Neutral = mean_zyg_shift_neutral,
                         Mean_Zygosity_Shift_OneCopy = mean_zyg_shift_1copy,
                         stringsAsFactors=F)
    return(out.df)
}

#' label_suspicious_adjustments
#'
#' Label a set of adjustments that appear suspicious for biological reasons
#'
#' @param df_in input data.frame containing summary statistics for a set of adjustments
#' @param max_homozygous_prop_genome Maximum proportion of the genome called as homozygously deleted before an adjustment is tagged as 'very suspicious' (default 0.05)
#' @param max_exons_seg_homozygous Maximum number of exons in a homozygously deleted segment before an adjustment is tagged as 'very suspicious' (default 200)
#' @param max_avg_zyg_shift_neutral Maximum average zygosity shift for segments called as neutral before an adjustment is tagged as 'suspicious' (default 0.10)
#' @param min_avg_zyg_shift_onecopy Maximum difference between predicted value of zygosity shift for 1-copy regions and observed mean value (default 0.10)
#' @return the same data.frame with an additional column named Suspicious taking values 0, 1, or 2.  1 means suspicious; 2 means very suspicious.
#' @export
label_suspicious_adjustments <- function(df_in, max_homozygous_prop_genome = 0.05,
                                         max_exons_seg_homozygous = 200,
                                         max_avg_zyg_shift_neutral = 0.1,
                                         min_avg_zyg_shift_onecopy = 0.1) {
    df_out <- df_in
    df_out$Suspicious <- rep(0,nrow(df_out))
    df_out$Suspicious[which(df_out$Mean_Zygosity_Shift_Neutral > max_avg_zyg_shift_neutral)] <- 1
    thresh1 <- 1/(2-df_out$Estimated_Purity)-0.5
    df_out$Suspicious[which(df_out$Mean_Zygosity_Shift_OneCopy <= thresh1 - min_avg_zyg_shift_onecopy)] <- 2
    df_out$Suspicious[which(df_out$Max_Segment_Homozygous > max_exons_seg_homozygous)] <- 2
    df_out$Suspicious[which(df_out$Proportion_Homozygous > max_homozygous_prop_genome)] <- 2
    return(df_out)
}

#' get_model_info
#'
#' Get information for a set of models
#' @param models list of models
#' @param adj_df data.frame containing adjustment information
#' @param coverage_df data.frame containing coverage data
#' @param seg_data data.frame containing segmented data
#' @param numproc number of processors to use
#' @param colname column name to use
#' @return list of classifications and summary stats
#' @export
get_model_info <- function(models, adj_df, coverage_df, seg_data, numproc, colname) {
    all_fits <- models2fits(models = models, adjustments = adj_df$Adjustment)
    ## filter to plausible adjustments
    possible_adjustments <- unique(c(filter_to_plausible_adjustments(all_fits),0))
    S1 <- subset(all_fits, all_fits$Adjustment %in% possible_adjustments)
    S1 <- rbind(S1,subset(S1,Adjustment==0)); rownames(S1) <- NULL
    classif.list <- get_classifications_from_models(models = models,
                                                    adj_df = S1,
                                                    cov_df = coverage_df,
                                                    seg_df = seg_data,
                                                    colname = colname,
                                                    numproc = numproc)
    S1 <- cbind(S1, get_summary_stats_from_classifications(classif.list))
    S1 <- label_suspicious_adjustments(S1)
    ## second version of 0-adjusted model is all-neutral classifications; change classif.list to reflect this
    atmp <- classif.list[[length(classif.list)]]$classifs
    atmp$Classification <- rep("2 Copies, Normal Genotype",nrow(atmp))
    atmp$Probability <- rep(1,nrow(atmp))
    atmp$Other.States <- rep(NA,nrow(atmp))
    atmp$Copy.Number <- rep(2,nrow(atmp))
    atmp$GainLoss <- rep("Neutral",nrow(atmp))
    classif.list[[length(classif.list)]]$classifs <- atmp
    ## All_Neutral column is all zero except for this one model
    S1$All_Neutral <- rep(0,nrow(S1)); S1$All_Neutral[nrow(S1)] <- 1
    for (cname in setdiff(colnames(S1),c("Adjustment","Ploidy","All_Neutral"))) {
        S1[which(S1$All_Neutral==1),cname] <- NA
    }
    ## order the adjustments (least-suspicious first, then ordered by increasing Fit_Combined_Penalized.
    ## all-neutral is always last)
    o1 <- order(S1$All_Neutral, S1$Suspicious, S1$Fit_Combined_Penalized, decreasing=FALSE)
    S1 <- S1[o1,]; rownames(S1) <- NULL; classif.list <- classif.list[o1]
    S1$Model_Number <- 1:nrow(S1); S1$Model_Name <- paste0("model_",S1$Model_Number)
    return(list(classif.list = classif.list, summary_stats = S1))
}

#' output_model_info
#'
#' Output model information
#' @param coverage_data data.frame containing coverage
#' @param coverage_data_col column name to use for coverage
#' @param zygosity_data data.frame containing zygosity
#' @param adjustment ploidy adjustment
#' @param chrinfo chromosome information data.frame
#' @param classification_df data.frame of classifications
#' @param model_outdir directory to which to write output files for model
#' @param prefix sample prefix to use for all result files
#' @param tc tumor content value for the model
#' @param locus_df data.frame of gene annotations to use for annotating the results (default NULL - no annotations)
#' @param exceptional_regions data.frame of exceptional regions
#' @param gender gender of sample (default female)
#' @export
output_model_info <- function(coverage_data, coverage_data_col, zygosity_data, adjustment,
                              chrinfo, classification_df, model_outdir, prefix, tc, locus_df = NULL,
                              exceptional_regions = NULL, gender = "female") {
    suppressPackageStartupMessages(require(IRanges))
    system(paste0("mkdir -p ",model_outdir))
    cat(adjustment,sep="\n",file=paste0(model_outdir,"/",prefix,".adjustment.txt"))
    coverage_data_adj <- coverage_data
    coverage_data_adj[,coverage_data_col] <- coverage_data_adj[,coverage_data_col] + adjustment
    ## non-PAR regions on X and Y for males have half the reported number of copies, modify this here
    w.modify <- integer(0)
    if (gender == "male") {
        w.modify <- which(classification_df$Chr %in% c("X","Y"))
        data(PAR_coords_hg19)
        w.par <- integer(0)
        w.par <- c(w.par, which(classification_df$Chr=="X" & classification_df$End <= PAR_coords_hg19$End[1]))
        w.par <- c(w.par, which(classification_df$Chr=="X" & classification_df$Start >= PAR_coords_hg19$Start[2]))
        w.modify <- setdiff(w.modify, w.par)
    }
    if (length(w.modify) > 0) {
        s1 <- strsplit(as.character(classification_df$Other.States[w.modify]),";")
        classification_df$Other.States[w.modify] <- unlist(lapply(s1, function(x) {
            ## num copies
            y <- unlist(lapply(strsplit(x,", "), function(a) a[1]))
            y2 <- unlist(lapply(strsplit(y," "), function(a) as.numeric(a[1])/2))
            y3 <- unlist(lapply(strsplit(y," "), function(a) a[2]))
            ## probs
            z <- unlist(lapply(strsplit(x," "),function(a) a[length(a)]))
            paste(paste(y2,y3,z),collapse=";")
        }))
        classification_df$Copy.Number[w.modify] <- classification_df$Copy.Number[w.modify]/2
        classification_df$GainLoss[w.modify] <- ifelse(classification_df$Copy.Number[w.modify] > 1, "Gain", 
                          ifelse(classification_df$Copy.Number[w.modify] < 1, "Loss", "Neutral"))
        classification_df$Classification[w.modify] <-
            unlist(lapply(strsplit(classification_df$Classification[w.modify]," "),
                          function(a) paste0(as.numeric(a[1])/2," Copies")))
        classification_df$Classification[which(classification_df$Classification == "1 Copies")] <- "1 Copy"
    }
    coverage_data_adj <- append_colors(df_in = coverage_data_adj,
                                       classif_df = classification_df,
                                       type = "coverage",
                                       exceptional_regions = exceptional_regions)
    zygosity_data <- append_colors(df_in = zygosity_data,
                                   classif_df = classification_df,
                                   type = "zygosity",
                                   exceptional_regions = exceptional_regions)
    cat(ifelse(is.na(tc),0,tc), sep="\n", file = paste0(model_outdir,"/",prefix,".tumorcontent.txt"))
    plot.file.0 <- paste0(model_outdir,"/",prefix,".coverage_zygosity_plot.png")
    if (is.na(tc)) { tc.title <- "Low Tumor Content: No Aberrations Called" }
    if (!is.na(tc)) { tc.title <- paste0("Estimated Tumor Content: ",round(100*tc),"%") }
    png(plot.file.0,width=1500,height=600,type="cairo")
    par(mfrow=c(2,1))
    plot_copynumber(coverage_data_adj, chrinfo, "coverage", coverage_data_col, NULL, pointsize=1.25)
    plot_copynumber(zygosity_data, chrinfo, "zygosity", NULL, tc.title, pointsize=1.25)
    dev.off()
    write.table(coverage_data_adj, file = paste0(model_outdir,"/",prefix,".coverage_data_classified.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(zygosity_data, file = paste0(model_outdir,"/",prefix,".zygosity_data_classified.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(classification_df, file = paste0(model_outdir,"/",prefix,".classified-segments.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    if (!is.null(locus_df)) {
        locus_df$Ind <- 1:nrow(locus_df)
        cn_by_locus <- do.call(rbind,
                               lapply(unique(classification_df$Chr), function(chrname) {
                                   wA <- which(classification_df$Chr==chrname)
                                   wR <- which(locus_df$Chr==chrname)
                                   IA <- IRanges(start=classification_df$Start[wA],end=classification_df$End[wA])
                                   IR <- IRanges(start=locus_df$Start[wR],end=locus_df$End[wR])
                                   FRA <- findOverlaps(query = IR, subject = IA)
                                   q1 <- queryHits(FRA)
                                   s1 <- subjectHits(FRA)
                                   data.frame(Locus = locus_df$Locus[wR[q1]],
                                              Copy_Number = classification_df$Copy.Number[wA[s1]],
                                              stringsAsFactors = FALSE)
                               }))
        ## average multiple values and add missing ones
        t1 <- tapply(cn_by_locus$Copy_Number, cn_by_locus$Locus, mean)
        cn_by_locus <- data.frame(Locus = names(t1), Copy_Number = t1, stringsAsFactors = F)
        cn_by_locus <- merge(locus_df, cn_by_locus, by = "Locus", all.x=T, all.y=F)
        cn_by_locus <- cn_by_locus[order(cn_by_locus$Ind),]
        cn_by_locus$Ind <- NULL
        rownames(cn_by_locus) <- NULL
        write.table(cn_by_locus, file = paste0(model_outdir,"/",prefix,".copynumber_by_locus.txt"),
                    quote = FALSE, sep = "\t", row.names = FALSE)
    }
}

#' integrated_classification
#'
#' Ploidy-adjusted classification across a set of models
#' @param covdata data.frame containing coverage
#' @param zygdata data.frame containing zygosity
#' @param segdata data.frame containing segmentation results
#' @param outdir output directory
#' @param refsample reference sample ('normal' or 'pool')
#' @param prefix sample prefix to use for result files
#' @param numproc Number of processors to use
#' @param locus_annotation_file data.frame containing gene structure for annotation purposes
#' @param exceptional_regions regions not to classify
#' @export
integrated_classification <- function(covdata, zygdata, segdata, outdir, refsample,
                                      prefix = "Sample", numproc = 2, locus_annotation_file = NULL,
                                      exceptional_regions = NULL, gender = "female") {
    ChrInfo <- get_chrinfo(chrnames = c(1:22,"X","Y"), numproc = numproc)
    if (refsample == "normal") {
        cname <- "LR.Corrected.TN"; subdir <- paste0(outdir,"/classification_models_TN"); cname2 <- "MeanCovN"
    }
    if (refsample == "pool") {
        cname <- "LR.Corrected.TP"; subdir <- paste0(outdir,"/classification_models_TP"); cname2 <- "MeanCovT"
    }
    system(paste0("mkdir -p ",subdir))
    g <- get_tc_lower_limit(coverage_df = covdata, colname_LR = cname, colname_cov = cname2)
    qc_stats <- data.frame(Lower_Limit=g[1],Sigma=g[2],stringsAsFactors=F)
    ws <- ifelse(refsample=="normal","TN","TP")
    qc_file <- paste0(outdir,"/classification_models_",ws,"/",prefix,".QC.",ws,".txt")
    write.table(qc_stats, file = qc_file, quote = FALSE, sep = "\t", row.names = FALSE)
    adjustment.df <- initialize_adjustments(seg_data = segdata, tc_lower = g[1])
    all_models <- fit_models(adj_df = adjustment.df,
                             coverage_data = covdata,
                             coverage_data_colname = cname,
                             seg_data = segdata,
                             zygosity_data = zygdata,
                             numproc = numproc,
                             chrinfo = ChrInfo,
                             exceptional_regions = exceptional_regions,
                             gender = gender)
    model_info_list <- get_model_info(models = all_models, adj_df = adjustment.df, coverage_df = covdata,
                                      seg_data = segdata, numproc = numproc, colname = cname)
    ## output classifications, tumor content estimates, and plots for each model.
    S1 <- model_info_list$summary_stats
    write.table(S1, file = paste0(subdir,"/model-summary-stats.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    if (is.null(locus_annotation_file)) {
        locus_df <- NULL
    } else {
        locus_df <- read.delim(locus_annotation_file, stringsAsFactors=F)
    }
    for (j in 1:nrow(S1)) {
        output_model_info(coverage_data = covdata,
                          coverage_data_col = cname,
                          zygosity_data = zygdata,
                          adjustment = S1$Adjustment[j],
                          chrinfo = ChrInfo,
                          classification_df = model_info_list$classif.list[[j]]$classifs,
                          model_outdir = paste0(subdir,"/",S1$Model_Name[j]),
                          prefix = prefix,
                          tc = S1$Estimated_Purity[j],
                          locus_df = locus_df,
                          exceptional_regions = exceptional_regions,
                          gender = gender)
    }
}
