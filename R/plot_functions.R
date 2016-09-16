#' get_plotting_coords
#'
#' @param position_vec vector containing chromosomal positions
#' @param chr_vec vector containing chromosome names (from 1-22,X,Y)
#' @param numproc number of processors to use (default 2)
#' @return a vector of correctly-spaced coordinates to use for plotting along the whole genome, in units of 10kb
#'
#' @export
get_plotting_coords <- function(position_vec, chr_vec, numproc = 2) {
    suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg19))
    suppressPackageStartupMessages(require(parallel))
    options(scipen=99)
    ChrInfo <- data.frame(ChrName=c(1:22,"X","Y"),stringsAsFactors=F)
    m1 <- mclapply(as.list(ChrInfo$ChrName), function(chrname) {
        chrname_full <- paste0("chr",chrname)
        return(length(Hsapiens[[chrname_full]]))
    }, mc.cores = numproc)
    ChrInfo$Length <- as.numeric(unlist(m1))
    ChrInfo$Length <- ChrInfo$Length + 20000000
    ChrInfo$Length <- ChrInfo$Length/10000
    ChrInfo$PrevLengths <- c(0,cumsum(ChrInfo$Length[-nrow(ChrInfo)]))
    N <- subset(ChrInfo,ChrName=="Y")$PrevLengths+subset(ChrInfo,ChrName=="Y")$Length
    a <- 0.5*(N+1)*(1-1/1.06)
    b <- 0.5*(N+1)*(1+1/1.06)
    end_to_end <- position_vec/10000 + ChrInfo$PrevLengths[match(chr_vec, ChrInfo$ChrName)]
    return(end_to_end)
}

#' plot_copynumber
#'
#' Make a copy number plot
#'
#' @param df_in input data.frame
#' @param chrinfo chromosome information data.frame
#' @param type 'coverage' or 'zygosity'
#' @param coverage_colname column name to use for coverage
#' @param title_text text to use for the title of the plot
#' @param pointsize size of points in plot, default 1 (cex units)
#' @export
plot_copynumber <- function(df_in, chrinfo, type = "coverage", coverage_colname = NULL,
                            title_text = NULL, pointsize=1) {
    par(mar=c(1.6,4.1,4.1,4.1))
    if (type == "coverage") {
        y <- df_in[,coverage_colname]
        YLIM <- c(-4,4)
        YLIM_EXP <- c(-4.32,4.32)
        AXAT <- seq(-4,4,by=2)
        YLABS <- c("Log2 Coverage Ratio","Coverage Ratio")
        AXLAB2 <- AXAT; AXLAB4 <- 2^(AXAT)
    }
    if (type == "zygosity") {
        y <- df_in$Var.T/df_in$Tot.T
        YLIM <- c(0,1)
        YLIM_EXP <- c(-0.04,1.04)
        AXAT <- seq(0,1,by=0.2)
        AXLAB2 <- AXLAB4 <- paste0(seq(0,100,by=20),"%")
        YLABS <- rep("Tumor Allelic Fraction",2)
    }
    if (type == "coverage") y.modified <- ifelse(y >= 4.15,4.15,ifelse(y <= -4.15,-4.15,y))
    if (type == "zygosity") y.modified <- ifelse(y >= 1,1,ifelse(y <= 0,0,y))
    x <- df_in$Plot_Coord
    N <- subset(chrinfo,ChrName=="Y")$PrevLengths+subset(chrinfo,ChrName=="Y")$Length
    a <- 0.5*(N+1)*(1-1/1.06)
    b <- 0.5*(N+1)*(1+1/1.06)
    plot(x,y.modified,type="n",xlim=c(a,b),ylab=NA,xlab=NA,xaxt="n",yaxt="n",ylim=YLIM,las=2)
    for (i in (1:12)*2) rect(chrinfo$PrevLengths[i],YLIM_EXP[1],chrinfo$PrevLengths[i]+chrinfo$Length[i]-2000,
                             YLIM_EXP[2],col="grey90",border=NA)
    points(x,y.modified,pch=16,cex=0.4*pointsize,col=df_in$Color)
    axis(4,at=AXAT,labels=AXLAB4,las=2)
    axis(2,at=AXAT,labels=AXLAB2,las=2)
    mtext(YLABS[1],side=2,line=3,font=2); mtext(YLABS[2],side=4,line=3,font=2)
    mtext(chrinfo$ChrName,at=0.5*chrinfo$Length+chrinfo$PrevLengths-1000,cex=1,font=2,col="black",line=0.5)
    if (!is.null(title_text)) mtext(title_text, side = 3, line = 2.5, font = 2, cex = 1.25)
    box()
}

#' append_colors
#'
#' Append colors to a data.frame of classifications
#' @param df_in input data.frame
#' @param classif_df data.frame of classifications
#' @param type 'coverage' or 'zygosity'
#' @param exceptional_regions data.frame containing exceptional regions (Chr, Start, End)
#' @param exceptional_regions_col color to use for exceptional regions, default orange
#' @return data.frame with additional column 'Color' (navyblue, blue, green, black, red2, red4)
#' @export
append_colors <- function(df_in, classif_df, type = "coverage",
                          exceptional_regions = NULL, exceptional_regions_col = "orange") {
    suppressPackageStartupMessages(require(IRanges))
    df_out <- df_in; df_out$Color <- rep(NA,nrow(df_out))
    for (chrom in unique(df_out$Chr)) {
        wC1 <- which(df_out$Chr==chrom)
        wclass1 <- which(classif_df$Chr==chrom)
        if (type == "coverage") IC <- IRanges(start=df_out$Start[wC1],end=df_out$End[wC1])
        if (type == "zygosity") IC <- IRanges(start=df_out$Pos[wC1],end=df_out$Pos[wC1])
        IB <- IRanges(start=classif_df$Start[wclass1],end=classif_df$End[wclass1])
        F1 <- findOverlaps(query=IC,subject=IB)
        cn1 <- classif_df$Copy.Number[wclass1[subjectHits(F1)]]
        cols1 <- rep("black",length(wC1[queryHits(F1)]))
        cols1[which(cn1==0)] <- "navyblue"
        ##cols1[which(cn1==1)] <- "blue"
        cols1[which(cn1==1 & classif_df$GainLoss[wclass1[subjectHits(F1)]]=="Loss")] <- "blue"
        cols1[which(cn1==1 & classif_df$GainLoss[wclass1[subjectHits(F1)]]=="Neutral")] <- "black"
        cols1[which(cn1==2 & classif_df$GainLoss[wclass1[subjectHits(F1)]]=="Loss")] <- "green"
        cols1[which(cn1==2 & classif_df$GainLoss[wclass1[subjectHits(F1)]]=="Neutral")] <- "black"
        cols1[which(cn1==2 & classif_df$GainLoss[wclass1[subjectHits(F1)]]=="Gain")] <- "red2"
        cols1[which(cn1 %in% c(3,4))] <- "red2"
        cols1[which(cn1 >= 5)] <- "red4"
        df_out$Color[wC1[queryHits(F1)]] <- cols1
    }
    if (!is.null(exceptional_regions)) {
        if (nrow(exceptional_regions) > 0) {
            for (i in 1:nrow(exceptional_regions)) {
                if (type == "coverage") {
                    wi <- which(df_out$Chr == exceptional_regions$Chr[i] &
                                df_out$End <= exceptional_regions$End[i] &
                                df_out$Start >= exceptional_regions$Start[i])
                }
                if (type == "zygosity") {
                    wi <- which(df_out$Chr == exceptional_regions$Chr[i] &
                                df_out$Pos <= exceptional_regions$End[i] &
                                df_out$Pos >= exceptional_regions$Start[i])
                }
                if (length(wi) > 0) {
                    df_out$Color[wi] <- exceptional_regions_col
                }
            }
        }
    }
    return(df_out)
}

#' make_germline_zygosity_plot
#'
#' Make a germline zygosity plot
#' @param zygdata Zygosity data.frame
#' @param polymorphic_regionsfile file containing polymorphic regions to exclude
#' @param outdir output directory to which to write plot (which goes to /germline subdirectory)
#' @param prefix sample prefix to use
#' @param chrinfo chromosome information data.frame
#' @export
make_germline_zygosity_plot <- function(zygdata, polymorphic_regionsfile, outdir, prefix, chrinfo) {
    suppressPackageStartupMessages(require(IRanges))
    system(paste0("mkdir -p ",outdir,"/germline"))
    if (!is.null(polymorphic_regionsfile)) {
        poly_regions <- read.delim(polymorphic_regionsfile, stringsAsFactors=F)
        inds_remove <- NULL
        for (chrname in c(1:22,"X","Y")) {
            wC <- which(zygdata$Chr==chrname)
            wG <- which(poly_regions$Chr==chrname)
            IC <- IRanges(start=zygdata$Pos[wC],end=zygdata$Pos[wC])
            IG <- IRanges(start=poly_regions$Start[wG],poly_regions$End[wG])
            F1 <- findOverlaps(query=IG,subject=IC)
            inds_remove <- c(inds_remove,wC[subjectHits(F1)])
        }
        zygdata$Polymorphic <- rep(0,nrow(zygdata))
        zygdata$Polymorphic[inds_remove] <- 1
        zygdata <- subset(zygdata, Polymorphic == 0)
    }
    zygdata$Color <- rep("black",nrow(zygdata))
    zygdata <- subset(zygdata,Pop_Freq > 0 & Tot.N >= 30)
    zygdata$Var.T <- zygdata$Var.N; zygdata$Tot.T <- zygdata$Tot.N
    png(paste0(outdir,"/germline/",prefix,".germline_zygosity.png"),
        width=1500,height=300,type="cairo")
    plot_copynumber(zygdata, chrinfo, "zygosity",
                    title_text="Germline Zygosity Data", pointsize=1.25)
    dev.off()
}
