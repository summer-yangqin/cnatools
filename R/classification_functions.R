#' create.statelist
#'
#' Create a list of copy number/zygosity states for a given tumor purity value
#'
#' @param theta Tumor purity value (between 0 and 1)
#' @param copyrange Vector of integer copy number values to consider (default 0:4)
#' @param lr.mean Log2 coverage ratio mean (default 0)
#' @param vf.mean Variant fraction mean (default 0)
#' @return data.frame containing columns State, Copies, Num.Var, Normal.Genotype, Tumor.Genotype, LR, and VF
#'
#' @export
create.statelist <- function(theta, copyrange=0:4, lr.mean=0, vf.mean=0.5) {
    cn.vect <- genotype.vect <- NULL
    for (cn in copyrange) {
        x <- (ceiling(cn/2)):cn
        cn.vect <- c(cn.vect, rep(cn, length(x)))
        genotype.vect <- c(genotype.vect, x)
    }
    D1 <- data.frame(State=1:length(cn.vect), Copies=cn.vect, Num.Var=genotype.vect, 
                     Normal.Genotype="A/B", Tumor.Genotype=NA, LR=NA, VF=NA, stringsAsFactors=F)
    D1$LR <- log2(theta*D1$Copies+(1-theta)*2)-1-lr.mean
    D1$VF <- (theta*D1$Num.Var+(1-theta)*1)/(theta*D1$Copies + (1-theta)*2)-(0.5-vf.mean)
    for (i in 1:nrow(D1)) D1$Tumor.Genotype[i] <- paste(c(rep("A", D1$Copies[i]-D1$Num.Var[i]), 
                                                          rep("B", D1$Num.Var[i])), collapse="/")
    return(D1)
}

#' get.fb.lik
#'
#' Retrieve likelihood for zygosity data
#'
#' @param x Number of variant reads (binomial number of successes)
#' @param n Amount of coverage (binomial number of experiments)
#' @param p Zygosity mean (binomial probability)
#' @return Likelihood using the beta-binomial function
#'
#' @export
get.fb.lik <- function(x, n, p) {
    suppressPackageStartupMessages(require(VGAM))
    dbetabinom(x, size=n, prob=p, rho=0.005)+(ifelse(x==n/2, 0, 1))*dbetabinom(n-x, size=n, prob=p, rho=0.005)
}

#' e.step
#'
#' Perform E-step of EM algorithm
#'
#' @param theta Tumor content estimate
#' @param segs data.frame containing segments
#' @param snplist List of snp data per segment, in the same order as segs data.frame
#' @param coveragelist List of coverage data per segment, in the same order as segs data.frame
#' @param coverage_colname Column name to use from coverage data.frame
#' @param stdev Global standard deviation of log2 coverage ratios
#' @param copyrange Vector of integer copy number states to consider (default 0:4)
#' @param lr.mean Global Log2 coverage ratio mean (default 0)
#' @param vf.mean Global Variant Fraction mean (default 0.5)
#' @return a list with two components, Probs and States
#'
#' @export
e.step <- function(theta, segs, snplist, coveragelist, coverage_colname, stdev,
                   copyrange=0:4, lr.mean=0, vf.mean=0.5) {
    States <- create.statelist(theta=theta, copyrange=copyrange, lr.mean=lr.mean, vf.mean=vf.mean)
    M <- matrix(NA, nrow=nrow(segs), ncol=nrow(States)); colnames(M) <- paste("Class.", States$State, sep="")
    if (nrow(segs) > 0) {
        for (i in 1:nrow(segs)) {
            ## compute log-likelihood over each of the states		
            ## variant fraction log-likelihood
            t1 <- tapply(States$VF, States$State, 
                         function(vf) {
                             aaa <- get.fb.lik(x=pmax(snplist[[i]]$Var.T, snplist[[i]]$Tot.T-snplist[[i]]$Var.T), 
                                               n=snplist[[i]]$Tot.T, p=vf)
                             return(sum(log(aaa[which(aaa != 0)])))
                         })
            ## coverage ratio log-likelihood
            t2 <- tapply(States$LR, States$State, function(lr) {
                aaa <- dnorm(coveragelist[[i]][, coverage_colname], mean=lr, sd=stdev)
                return(sum(log(aaa[which(aaa != 0)])))
            })
            t1 <- ifelse(t1==0, NA, t1)
            t2 <- ifelse(t2==0, NA, t2)
            total.loglik <- t1 + t2
            M[i, ] <- total.loglik
        }
        w <- which(apply(is.na(M), 1, sum) < ncol(M))
        M2 <- matrix(NA, nrow=nrow(M), ncol=ncol(M)); colnames(M2) <- colnames(M)
        if (length(w)==1) {
            M2[w, ] <- exp(M[w, ] - max(M[w,], na.rm=T))
        } else {
            M2[w, ] <- exp(t(apply(M[w, ], 1, function(x) x-max(x, na.rm=T))))
        }
        M <- ifelse(is.na(M2), 0, M2)
        ## monotonicity constraints:
        ## segments with Log-Ratios < 0 have all probabilities of 4 or more copies == 0.
        wseg <- which(segs$Log2.Coverage.Ratio < 0)
        wstate <- which(States$Copies >= 4)
        for (w1 in wseg) {
            for (w2 in wstate) {
                M[w1, w2] <- 0
            }
        }
        ## segments with Log-Ratios > 1 have probabilities of < 2 copies == 0.
        wseg <- which(segs$Log2.Coverage.Ratio > 1)
        wstate <- which(States$Copies < 2)
        for (w1 in wseg) {
            for (w2 in wstate) {
                M[w1, w2] <- 0
            }
        }
        ## renormalize probabilities to sum to 1
        M <- round(sweep(M, 1, apply(M, 1, sum), FUN="/"), 3)
        ## if probabilities are missing, classify by distance from closest state
        w <- which(apply(M, 1, function(x) sum(is.na(x))) > 0)
        if (length(w) > 0) {
            for (j in 1:length(w)) {
                snplist.problematic <- snplist[[w[j]]]
                crlist.problematic <- coveragelist[[w[j]]]
                mean.zyg.frac <- mean(0.5+abs(snplist.problematic$Var.T/snplist.problematic$Tot.T-0.5), na.rm=T)
                mean.l2.cr <- mean(crlist.problematic[, coverage_colname], na.rm=T)
                Atemp <- matrix(c(mean.zyg.frac, mean.l2.cr), nrow(States), 2, byrow=T)			
                atemp <- apply(States[, 6:7]-Atemp, 1, function(x) sum(x^2))
                wmin <- which.min(atemp)
                M[w[j], ] <- rep(0, ncol(M)); M[w[j], wmin] <- 1
            }
        }
    }
    return(list(Probs=M, States=States))
}

#' m.step
#'
#' Perform M-Step of EM algorithm
#'
#' @param E.step.results List representing output of e.step function
#' @param theta.seq Vector of tumor content values to consider
#' @param snplist List of snp data per segment, in the same order as segs data.frame
#' @param coveragelist List of coverage data per segment, in the same order as segs data.frame
#' @param coverage_colname Column name to use from coverage data.frame
#' @param stdev Global standard deviation of log2 coverage ratios
#' @param copyrange Vector of integer copy number states to consider (default 0:4)
#' @return Value of theta (tumor purity) maximizing the likelihood
#'
#' @export
m.step <- function(E.step.results, theta.seq, snplist, coveragelist, coverage_colname, stdev, copyrange=0:4) {
    loglik <- NULL
    for (j in 1:length(theta.seq)) {
        classif.vect <- apply(E.step.results$Probs, 1, which.max)
        States <- create.statelist(theta=theta.seq[j], copyrange=copyrange)
        vf.loglik <- 0
        if (length(snplist) > 0) {
            for (i in 1:length(snplist)) {
                aaa <- get.fb.lik(x=pmax(snplist[[i]]$Var.T, snplist[[i]]$Tot.T-snplist[[i]]$Var.T), 
                                  n=snplist[[i]]$Tot.T, 
                                  p=States$VF[classif.vect[i]])
                vf.loglik <- vf.loglik+sum(log(aaa[which(aaa != 0)]))
            }
        }
        lr.loglik <- 0
        if (length(coveragelist) > 0) {
            for (i in 1:length(coveragelist)) {
                aaa <- dnorm(coveragelist[[i]][, coverage_colname], 
                             mean=States$LR[classif.vect[i]], sd=stdev)
                lr.loglik <- lr.loglik+sum(log(aaa[which(aaa != 0)]))
            }
        }
        loglik[j] <- vf.loglik+lr.loglik
    }
    return(theta.seq[which.max(loglik)])
}

#' overall_classify
#'
#' Overall classification function
#'
#' @param Z.data zygosity data.frame
#' @param C.data coverage data.frame
#' @param C.data.colname Column name to use from coverage data.frame 
#' @param seg.data data.frame containing segments
#' @param chrinfo data.frame containing chromosome info
#' @param lower.limit.tc Lower limit to consider for tumor content estimation (default = 0.3)
#' @param msg Message to print (default "")
#' @param exceptional_regions A set of regions that do not get classified (data.frame, columns Chr, Start, End. default NULL)
#' @param gender gender of sample (default female)
#' @return a list with two components, classifs (data.frame of classifications) and tc (final inferred tumor content)
#'
#' @export
overall_classify <- function(Z.data, C.data, C.data.colname, seg.data, chrinfo, lower.limit.tc=0.3,
                             msg="", exceptional_regions = NULL, gender = "female") {
    suppressPackageStartupMessages(require(DNAcopy))
    suppressPackageStartupMessages(require(IRanges))
    ## PREPARE DATA FOR EM ALGORITHM:
    C.data$Midpoint <- round(0.5*(C.data$Start+C.data$End), 0)
    ## use 50X coverage as a minimum threshold:
    C.data <- subset(C.data, MeanCovN >= 50)
    Z.data <- subset(Z.data, Tot.N >= 50)
    ## for log2 coverage ratio, assume gaussian distribution (might change this) with known standard deviation
    overall.LR.sd <- sd(diff(C.data[, C.data.colname]), na.rm=T)/sqrt(2)
    ## limit to segments with enough data
    EM.segs <- subset(seg.data, Num.Agg.SNPs >= 20 & Num.Exons.NM >= 20)
    if (!is.null(exceptional_regions)) {
        w.exceptional <- integer(0)
        if (nrow(exceptional_regions) > 0) {
            for (i in 1:nrow(exceptional_regions)) {
                w <- which(EM.segs$Chr == exceptional_regions$Chr[i] &
                           EM.segs$Start >= exceptional_regions$Start[i] &
                           EM.segs$End <= exceptional_regions$End[i])
                if (length(w) > 0) w.exceptional <- c(w.exceptional, w)
            }
        }
        if (nrow(EM.segs) > 0) EM.segs <- EM.segs[setdiff(1:nrow(EM.segs), w.exceptional),]
    }
    if (nrow(EM.segs) > 0) EM.segs$EM.Segment <- 1:nrow(EM.segs) else EM.segs$EM.Segment <- numeric(0)
    ## get the zygosity and coverage data on each segment - store in lists called snplist and coveragelist
    mu.z <- sum(Z.data$Var.T)/sum(Z.data$Tot.T)
    snplist <- coveragelist <- list()
    if (nrow(EM.segs) > 0) {
        for (i in 1:nrow(EM.segs)) {
            snplist[[i]] <- subset(Z.data, Chr==EM.segs$Chr[i] & Pos >= EM.segs$Start[i] & 
                                          Pos <= EM.segs$End[i])[, c("Pos", "Tot.T", "Var.T")]
            coveragelist[[i]] <- subset(C.data, Chr==EM.segs$Chr[i] & Start >= EM.segs$Start[i] & 
                                               End <= EM.segs$End[i])[, c("Midpoint", C.data.colname)]
            coveragelist[[i]] <- coveragelist[[i]][which(!is.na(coveragelist[[i]][, 2])), ]
        }
    }
    ## INITIALIZE EM ALGORITHM (theta.init)
    ## nominate regions of 1-copy loss
    ## log2 coverage ratio <= -0.2, zygosity shift of at least 0.08, >= 50 snps
    S1 <- subset(seg.data, Dev.Zyg.Seg >= 0.08 & Log2.Coverage.Ratio <= -0.2 & 
                          Num.Agg.SNPs >= 10 & Num.Exons >= 10)
    w1 <- which.max(S1$Dev.Zyg.Seg)[1]
    ## initial estimate of tumor content, as derived from this segment
    mu1 <- 0.5+S1$Dev.Zyg.Seg[w1]
    theta1 <- 2-(1/mu1)
    ## nominate copy neutral LOH regions
    ## abs(log2 coverage ratio) < 0.10, zygosity shift of at least 0.15, >= 50 snps
    S2 <- subset(seg.data, Dev.Zyg.Seg >= 0.15 & abs(Log2.Coverage.Ratio) <= 0.05 & 
                          Num.Agg.SNPs >= 10 & Num.Exons >= 10)
    w2 <- which.max(S2$Dev.Zyg.Seg)[1]
    ## initial estimate of tumor content, as derived from this segment
    mu2 <- 0.5+S2$Dev.Zyg.Seg[w2]
    theta2 <- (mu2-0.5)/.5
    ## nominate regions of 1-copy gain
    ## log2 coverage ratio >= 0.2, zygosity shift of at least 0.065 and at most 0.16666, >= 50 snps
    S3 <- subset(seg.data, Dev.Zyg.Seg >= 0.065 & Dev.Zyg.Seg <= 1/6 & 
                          Log2.Coverage.Ratio >= 0.2 & Num.Agg.SNPs >= 50)
    w3 <- which.max(S3$Dev.Zyg.Seg)[1] # change this???
    ## initial estimate of tumor content, as derived from this segment
    mu3 <- 0.5+S3$Dev.Zyg.Seg[w3]
    theta3 <- (mu3-0.5)/.5
    theta3 <- (2*mu3-1)/(1-mu3)
    ## nominate regions of 2-copy gain
    ## log2 coverage ratio >= 0.6, zygosity shift of more than 0.166666 but at most 0.25, >= 50 snps
    S4 <- subset(seg.data, Dev.Zyg.Seg > 1/6 & Dev.Zyg.Seg <= 0.25 &
                          Log2.Coverage.Ratio >= 0.6 & Num.Agg.SNPs >= 50)
    w4 <- which.max(S4$Dev.Zyg.Seg)[1]
    ## initial estimate of tumor content, as derived from this segment
    mu4 <- 0.5+S4$Dev.Zyg.Seg[w4]
    theta4 <- (2*mu4-1)/(2-2*mu4)
    theta.init.all <- c(theta1, theta2, theta3, theta4)
    theta.type <- c("1-copy loss", "copy neutral LOH", "1-copy gain", "2-copy gain")
    w <- which(!is.na(theta.init.all))
    if (length(w) > 0) {
	init.type <- theta.type[min(w)]; theta.init <- round(theta.init.all[min(w)], 2)
    } else {
	init.type <- "no empirical data"; theta.init <- NA
    }
    if (is.na(theta.init)) theta.init <- 0.4
    ## iteratively apply E step and M step until convergence
    if (!is.na(theta.init)) {
	if (theta.init >= 0.99) theta.init <- 0.99
	tol <- 0.01; num.iter <- 5
	theta.path <- c(-1, theta.init); theta.current <- theta.init
	N <- length(theta.path)
        while (length(theta.path) < num.iter + 1) {
            E.current <- e.step(theta=theta.current, segs=EM.segs, 
                                snplist=snplist, coveragelist=coveragelist, 
                                coverage_colname = C.data.colname, 
                                stdev=overall.LR.sd, copyrange=0:8)
            theta.current <- m.step(E.current, theta.seq=seq(lower.limit.tc, 0.99, by=0.01), 
                                    snplist=snplist, coveragelist=coveragelist, 
                                    coverage_colname = C.data.colname, 
                                    stdev=overall.LR.sd, copyrange=0:8)
            theta.current <- round(theta.current, 2)
            if (abs(theta.current-theta.path[length(theta.path)]) < tol) {
                theta.path <- c(theta.path, rep(theta.current, num.iter+1-length(theta.path)))
            } else {
                theta.path <- c(theta.path, theta.current)
            }
            N <- length(theta.path)
        }
        if (theta.path[N-1]==theta.path[N]) {
            theta.final <- theta.path[N]
        } else {
            theta.final <- theta.path[N] ##NA
        }
    } else {
	theta.path <- NA
	N <- 1
	theta.final <- NA
    }
    cat(msg, "tumor content convergence path: ", paste(theta.path[which(theta.path != -1)], collapse=" "),
        "\n", sep="")
    ## derive list of all segments together with classifications
    ## first, limit to segments with at least 1 exon with a non-missing coverage quantification, and
    ## with at least one aggregated snp
    seg.data <- subset(seg.data, Num.Agg.SNPs > 0 | Num.Exons.NM > 0); rownames(seg.data) <- NULL
    ## retrieve the zygosity and coverage data on each segment (not just the ones used for EM algorithm)
    mu.z <- sum(Z.data$Var.T)/sum(Z.data$Tot.T)
    snplist.all <- coveragelist.all <- list()
    if (nrow(seg.data) > 0) {
        for (i in 1:nrow(seg.data)) {
            snplist.all[[i]] <- subset(Z.data, Chr==seg.data$Chr[i] & Pos >= seg.data$Start[i] & 
                                              Pos <= seg.data$End[i])[, c("Pos", "Tot.T", "Var.T")]
            coveragelist.all[[i]] <- subset(C.data, Chr==seg.data$Chr[i] & Start >= seg.data$Start[i] & 
                                                   End <= seg.data$End[i])[, c("Midpoint", C.data.colname)]
            coveragelist.all[[i]] <- coveragelist.all[[i]][which(!is.na(coveragelist.all[[i]][, 2])), ]
            ## make sure that number of snps/exons are in agreement between snplist.all/coveragelist.all
            ## and seg.data
            seg.data$Num.SNPs[i] <- nrow(snplist.all[[i]])
            seg.data$Num.Exons.NM[i] <- nrow(coveragelist.all[[i]])
        }
    }
    ## classify all segments containing both snps and exons with coverage
    w.both <- which(seg.data$Num.SNPs > 0 & seg.data$Num.Exons.NM > 0)
    w.nosnps <- which(seg.data$Num.SNPs == 0 & seg.data$Num.Exons.NM > 0)
    w.nocoverage <- which(seg.data$Num.SNPs > 0 & seg.data$Num.Exons.NM == 0)
    ## E-step should give posterior probabilities for these segments
    E.final <- e.step(theta=theta.final, segs=seg.data[w.both, ], snplist=snplist.all[w.both], 
                      coveragelist=coveragelist.all[w.both], 
                      coverage_colname = C.data.colname, 
                      stdev=overall.LR.sd, copyrange=0:20)
    S.final <- cbind(seg.data, matrix(NA, nrow=nrow(seg.data), ncol=ncol(E.final$Probs)))
    colnames(S.final)[(ncol(seg.data)+1):ncol(S.final)] <- colnames(E.final$Probs)
    if (length(w.both) > 0) S.final[w.both, (ncol(seg.data)+1):ncol(S.final)] <- E.final$Probs
    ## for segments with coverage data but no snp data, use a likelihood-based approach using only the coverage data
    ## add a constraint that segments with mean lr-coverage < -0.2 can only get classified into states with LR <= 0 (sim for pos)
    for (i in w.nosnps) {
        c1 <- mean(coveragelist.all[[i]][, C.data.colname], na.rm=T)
        sign1 <- ifelse(c1 < -1, -2, ifelse(c1 < -0.2, -1, ifelse(c1 > 1, 2, ifelse(c1 > 0.2, 1, 0))))
        wsign <- 1:nrow(E.final$States)
        if (sign1 == -1) wsign <- which(E.final$States$LR <= 0)
        if (sign1 == 1) wsign <- which(E.final$States$LR >= 0)
        if (sign1 == -2) wsign <- which(E.final$States$LR < 0)
        if (sign1 == 2) wsign <- which(E.final$States$LR > 0)
        loglik <- rep(0, nrow(E.final$States))
	for (j in wsign) {
            aaa <- dnorm(coveragelist.all[[i]][, C.data.colname], 
                         mean=E.final$States$LR[j], sd=overall.LR.sd)
            loglik[j] <- sum(log(aaa[which(aaa != 0)]))
	}
        ## if sample is male and chromosome is 'X' or 'Y', only allow states with an even number of copies
        ## unless the segment is contained within PAR. (will halve number of copies at the end)
        if (seg.data$Chr[i] %in% c("X","Y") & gender == "male") {
            par <- 0
            data(PAR_coords_hg19)
            for (k in 1:nrow(PAR_coords_hg19)) {
                if (seg.data$Chr[i]==PAR_coords_hg19$Chr[k] &
                    seg.data$Start[i] >= PAR_coords_hg19$Start[k] &
                    seg.data$End[i] <= PAR_coords_hg19$End[k]) {
                    par <- 1
                }
            }
            if (par == 0) loglik[which(round(floor(E.final$States$Copies/2),1) != E.final$States$Copies/2)] <- 0
        }
	w.loglik <- which(loglik != 0)
	probs.tmp <- exp(loglik[w.loglik]-max(loglik[w.loglik]))
	probs.tmp <- round(probs.tmp/sum(probs.tmp), 3)
        probs <- rep(0, nrow(E.final$States)); probs[w.loglik] <- probs.tmp
        ## if there are no snps in the region, we can make no inference about genotype,
        ## so collapse multiple states with the same copy number and choose the first
	probs.cn <- tapply(probs, E.final$States$Copies, sum)
	m <- match(as.numeric(names(probs.cn)), E.final$States$Copies)
	S.final[i, paste("Class.", E.final$States$State, sep="")[m]] <- probs.cn
    }
    ## for segments with snp data but no coverage data, assign to neutral state. (no information on copy number)
    for (i in w.nocoverage) {
        ww <- which(E.final$States$Tumor.Genotype=="A/B")
        state.neutral <- paste0("Class.", E.final$States$State[ww])
        S.final[i, grep("Class.", colnames(S.final))] <- 0
        S.final[i, state.neutral] <- 1
    }
    ## classes with >= 5% posterior probability
    Ltemp.prelim <- apply(S.final[, grep("Class.", colnames(S.final))], 1, function(x) x[which(x >= 0.05)])
    ## if a segment exists with no probabilities of at least 5%, then reclassify using 0:6 rather than 0:20.
    w.not5 <- which(unlist(lapply(Ltemp.prelim,length))==0)
    if (length(w.not5)==0) {
        Ltemp <- Ltemp.prelim
    } else {
        for (i in w.not5) {
            ## note: fit E.final.2 with indices duplicated because of problems with apply and simplify
            E.final.2 <- e.step(theta=theta.final, segs=seg.data[c(i,i), ],
                                snplist=snplist.all[c(i,i)], 
                                coveragelist=coveragelist.all[c(i,i)], 
                                coverage_colname = C.data.colname, 
                                stdev=overall.LR.sd, copyrange=0:6)
            ## correspondence between states
            m <- match(paste(E.final.2$States$Copies,E.final.2$States$Num.Var,
                             E.final.2$States$Normal.Genotype,E.final.2$States$Tumor.Genotype),
                       paste(E.final$States$Copies,E.final$States$Num.Var,
                             E.final$States$Normal.Genotype,E.final$States$Tumor.Genotype))
            old.new.key <- data.frame(New.State = E.final.2$States$State,
                                      Old.State = E.final$States$State[m],stringsAsFactors=F)
            old.new.key$New.State <- paste0(rep("Class.",nrow(old.new.key)),old.new.key$New.State)
            old.new.key$Old.State <- paste0(rep("Class.",nrow(old.new.key)),old.new.key$Old.State)
            ## clear previous probabilities
            S.final[i,which(substr(colnames(S.final),1,6)=="Class.")] <- 0
            S.final[i,old.new.key$Old.State[match(colnames(E.final.2$Probs),old.new.key$New.State)]] <-
                E.final.2$Probs[1,]
        }
        Ltemp <- apply(S.final[, grep("Class.", colnames(S.final))], 1, function(x) x[which(x >= 0.05)])
    }
    if (class(Ltemp) != "list") {
        Ltemp2 <- as.list(Ltemp)
        if (length(Ltemp2) > 0) {
            for (i in 1:length(Ltemp2)) {
                ind <- which(S.final[i, grep("Class.", colnames(S.final))]==Ltemp2[[i]])
                names(Ltemp2[[i]]) <- colnames(S.final)[grep("Class.", colnames(S.final))][ind]
            }
        }
        Ltemp <- Ltemp2
    }
    Y.out <- S.final[, 1:10]; Y.out$Classification <- NA; Y.out$Probability <- NA
    Y.out$Other.States <- NA
    w <- which(unlist(lapply(Ltemp, length)) > 0)
    for (w0 in w) {
	probs <- Ltemp[[w0]]
	m <- match(names(probs), paste("Class.", E.final$States$State, sep=""))
	numcopies <- E.final$States$Copies[m]
	genotype <- E.final$States$Tumor.Genotype[m]
	numcopies.label <- paste(numcopies, ifelse(numcopies==1, "Copy", "Copies"), sep=" ")
        genotype.label <- ifelse(numcopies==0, "", 
                          ifelse(numcopies==2 & genotype=="B/B", "LOH", 
                          ifelse(numcopies==2 & genotype=="A/B", "Normal Genotype", 
                          ifelse(numcopies==1, "", 
                                 paste("Genotype", genotype, sep=" ")))))
	overall.label <- ifelse(genotype.label=="", numcopies.label, 
                                paste(numcopies.label, genotype.label, sep=", "))
	wmax <- which.max(probs)
	Y.out$Classification[w0] <- overall.label[wmax]
	Y.out$Probability[w0] <- probs[wmax]
	if (length(probs) > 1) {
            others <- paste(overall.label[-wmax], " (Prob=", probs[-wmax], ")", sep="")
            Y.out$Other.States[w0] <- paste(others, collapse=";")
	}
    }
    Y.out$Classification[which(is.na(Y.out$Classification))] <- "Ambiguous"
    w.end <- which(Y.out$End==300000001)
    m <- match(Y.out$Chr[w.end], chrinfo$ChrName)
    Y.out$End[w.end] <- chrinfo$Length_Orig[m]
    Y.out$Start[which(Y.out$Start==0)] <- 1
    Y.out$Region <- paste("chr", Y.out$Chr, ":", Y.out$Start, "-", Y.out$End, sep="")
    Y.out$Segment <- Y.out$Num.Agg.SNPs <- Y.out$Num.Exons <- NULL
    colnames(Y.out)[which(colnames(Y.out)=="Num.SNPs")] <- "Heterozygous.SNPs"
    colnames(Y.out)[which(colnames(Y.out)=="Dev.Zyg.Seg")] <- "Mean.Zygosity.Deviation"
    colnames(Y.out)[which(colnames(Y.out)=="Num.Exons.NM")] <- "Targeted.Exons"
    xtmp <- unlist(lapply(strsplit(Y.out$Classification, " Cop"), function(x) x[1]))
    Y.out$Copy.Number <- rep(NA, nrow(Y.out))
    Y.out$Copy.Number[which(xtmp != "Ambiguous")] <- as.numeric(xtmp[which(xtmp != "Ambiguous")])
    Y.out$GainLoss <- ifelse(Y.out$Copy.Number > 2, "Gain", 
                      ifelse(Y.out$Copy.Number < 2, "Loss", "Neutral"))
    Y.out$GainLoss[grep("LOH", Y.out$Classification)] <- "Loss"
    Y.out$Mean.Zygosity.Deviation <- round(Y.out$Mean.Zygosity.Deviation, 3)
    Y.out$Log2.Coverage.Ratio <- round(Y.out$Log2.Coverage.Ratio, 3)
    ## if there are exceptional regions, set their classifications to neutral
    if (!is.null(exceptional_regions)) {
        if (nrow(exceptional_regions) > 0) {
            for (i in 1:nrow(exceptional_regions)) {
                w <- which(Y.out$Chr == exceptional_regions$Chr[i] &
                           Y.out$Start >= exceptional_regions$Start[i] &
                           Y.out$End <= exceptional_regions$End[i])
                if (length(w) > 0) {
                    Y.out$Classification[w] <- "2 Copies, Normal Genotype"
                    Y.out$Probability[w] <- 1
                    Y.out$Other.States[w] <- NA
                    Y.out$Copy.Number[w] <- 2
                    Y.out$GainLoss[w] <- "Neutral"
                }
            }
        }
    }
    ## if (!is.na(theta.final)) {
    ##     final.classifications <- apply(S.final[, -(1:ncol(seg.data))], 1, which.max)
    ##     classification.probs <- apply(S.final[, -(1:ncol(seg.data))], 1, max)
    ##     final.classifications[which(classification.probs < .5)] <- 1000
    ##     cn.classification.probs <- t(apply(S.final[, -(1:ncol(seg.data))], 1, tapply, 
    ##                                        E.final$States$Copies, sum, na.rm=T))
    ##     cn.classifications <- E.final$States$Copies[apply(cn.classification.probs, 1, which.max)]
    ## } else {
    ##     final.classifications <- rep(E.final$States$State[which(E.final$States$Copies==2 &
    ##                                                             E.final$States$Num.Var==1)], nrow(S.final))
    ##     classification.probs <- rep(1, nrow(S.final))
    ##     cn.classifications <- rep(2, nrow(S.final)); cn.classification.probs <- rep(1, nrow(S.final))
    ## }
    ## ## if final.classifications is empty, add NA
    ## w <- which(unlist(lapply(final.classifications, length))==0)
    ## if (length(w) > 0) final.classifications[w] <- NA
    ## ## whole-genome coverage / zygosity plot:
    ## if (is.na(theta.final)) color.by.seg <- rep("black", nrow(seg.data))
    ## if (!is.na(theta.final)) {
    ##     ##numcopies.by.seg <- cn.classifications#rep(E.final$States$Copies[final.classifications])
    ##     numcopies.by.seg <- E.final$States$Copies[unlist(final.classifications)]
    ##     numvar.by.seg <- E.final$States$Num.Var[unlist(final.classifications)]
    ##     color.by.seg <- ifelse(numcopies.by.seg > 2, "red", ifelse(numcopies.by.seg < 2, "blue", 
    ##                                                       ifelse(numvar.by.seg==2, "green", "black")))
    ##     color.by.seg[which(is.na(color.by.seg))] <- "grey"
    ## }
    return(list(classifs=Y.out, tc=theta.final))
}

#' germline_classification
#'
#' Classification of germline variants
#' @param covdata coverage data.frame
#' @param segdata data.frame of segments
#' @param outdir output directory (results go to 'germline' subdirectory of outdir)
#' @param prefix sample prefix to append to beginning of each file
#' @param chrinfo chromosome information as generated by get_chrinfo function
#' @param gene_df data.frame containing genes with which to annotate output files (default NULL, in which case results are not annotated)
#' @export
germline_classification <- function(covdata, segdata, outdir, prefix, chrinfo, gene_df = NULL) {
    system(paste0("mkdir -p ",outdir,"/germline"))
    segdata$GainLoss <- ifelse(segdata$Log2_Coverage_Ratio >= 0.5,"Gain",
                     ifelse(segdata$Log2_Coverage_Ratio <= -0.85,"Loss","Neutral"))
    ## put classifications into per-exon coverage data.frame
    S2 <- covdata; S2$GainLoss <- rep(NA,nrow(S2))
    for (chrname in union(S2$Chr,segdata$chrom)) {
        wS2 <- which(S2$Chr==chrname)
        wS1 <- which(segdata$Chr==chrname)
        IS1 <- IRanges(start=segdata$Start[wS1],end=segdata$End[wS1])
        IS2 <- IRanges(start=S2$Start[wS2],end=S2$End[wS2])
        F1 <- findOverlaps(query=IS1,subject=IS2)
        S2$GainLoss[wS2[subjectHits(F1)]] <- segdata$GainLoss[wS1[queryHits(F1)]]
    }
    S2 <- subset(S2,Polymorphic==0 & !is.na(GainLoss))
    S2$Color <- ifelse(S2$GainLoss=="Gain","red",
                ifelse(S2$GainLoss=="Loss","blue","black"))
    write.table(S2, file = paste0(outdir,"/germline/",prefix,".germline_copynumber_by_exon.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    png(paste0(outdir,"/germline/",prefix,".germline_copynumber.png"),
        width=1500,height=300,type="cairo")
    plot_copynumber(S2, chrinfo, "coverage", "LR.Corrected.NP",
                    title_text="Germline Copy Number vs. Pooled Normal",pointsize=1.25)
    dev.off()
    segdata$Copy.Number <- ifelse(segdata$GainLoss=="Neutral",2,round(2^(segdata$Log2_Coverage_Ratio+1),0))
    segdata$Genes <- rep(NA,nrow(segdata))
    segdata <- subset(segdata, !is.na(GainLoss))
    if (!is.null(gene_df)) {
        for (chrname in unique(segdata$Chr)) {
            wA <- which(segdata$Chr==chrname)
            wR <- which(gene_df$Chr==chrname)
            IA <- IRanges(start=segdata$Start[wA],end=segdata$End[wA])
            IR <- IRanges(start=gene_df$Start[wR],end=gene_df$End[wR])
            FAR <- findOverlaps(query = IA, subject = IR)
            for (q1 in unique(queryHits(FAR))) {
                w <- which(queryHits(FAR)==q1)
                segdata$Genes[wA[q1]] <- paste(unique(gene_df$Gene[wR[subjectHits(FAR)[w]]]),collapse=";")
            }
        }
        L <- lapply(1:nrow(segdata), function(i) {
            data.frame(Copy.Number = segdata$Copy.Number[i],
                       GainLoss = segdata$GainLoss[i],
                       Genes = unlist(strsplit(segdata$Genes[i],";")),
                       Exons = segdata$Exons[i],
                       stringsAsFactors=F)
        })
        B <- do.call(rbind,L)
        B <- subset(B,!is.na(Genes))
        B <- subset(B,!is.na(Copy.Number))
        B <- B[order(abs(B$Copy.Number-2),-B$Copy.Number,B$Exons,decreasing=TRUE),]
        B <- B[which(!duplicated(B$Genes)),]
        colnames(B) <- ifelse(colnames(B)=="Genes","Gene",colnames(B))
        B <- B[,c("Gene","Copy.Number","GainLoss","Exons")]
        colnames(B) <- ifelse(colnames(B)=="Exons","Exons_In_Segment",colnames(B))
        B <- B[order(B$Gene),]
        rownames(B) <- NULL
    }
    write.table(segdata, file = paste0(outdir,"/germline/",
                                       prefix,".germline.classified-segments.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    if (!is.null(gene_df)) {
        write.table(B, file = paste0(outdir,"/germline/",
                                     prefix,".germline.copynumber_by_gene.txt"),
                    quote = FALSE, sep = "\t", row.names = FALSE)
    }
}
