#!/usr/bin/env Rscript

# hka_test2.R <mafs file> <max p-value for SNP> <min within-morph MAF for site to be variable> <outfile name>

# Performs HKA test for baned, striped, and banded+striped
# Update of 'hka_test.R' to discard all sites segregating in outgroup

## load libraries
library(snpStats) # for genomic control with qq.chisq

## parse data
cat("Reading in data\n")

args <- commandArgs(trailingOnly=TRUE)

mafs_file <- as.character(args[1])
maxp <- as.numeric(args[2])
minmaf <- as.numeric(args[3])
out_file <- as.character(args[4])

mafs <- read.table(mafs_file, head=TRUE) # read in data

## define functions

count <- function(snps = NULL) {
   # reassign rownames to be able to use indexes
   rownames(snps) <- 1:nrow(snps)
   snp.idx <- 1:nrow(snps)

   # sites variable among models
   vari.s <- which(snps$variabilis.maf %% 1 != 0)
   summer.s <- which(snps$summersi.maf %% 1 != 0)
   model.fixed.diff <- which((snps$variabilis.maf == 0 & snps$summersi.maf == 1) | (snps$summersi.maf == 0 & snps$variabilis.maf == 1))
   model.s <- sort(unique(c(vari.s, summer.s, model.fixed.diff)))

   # find sites that are variable within banded and striped morphs
   band.s <- which(snps$band.maf %% 1 != 0) # fixed sites should have allele frequency 0 or 1
   band.s <- band.s[which(band.s %in% model.s == FALSE)] # mask out sites segregating in outgroup
   band.n.s <- length(band.s)
   stripe.s <- which(snps$stripe.maf %% 1 != 0)
   stripe.s <- stripe.s[which(stripe.s %in% model.s == FALSE)] # mask out sites segregating in outgroup
   stripe.n.s <- length(stripe.s)

   # find fixed differences between banded and striped
   morph.fixed.diff <- which((snps$band.maf == 0 & snps$stripe.maf == 1) | (snps$stripe.maf == 0 & snps$band.maf == 1))
   morph.fixed.diff <- morph.fixed.diff[which(morph.fixed.diff %in% model.s == FALSE)] # remove sites polymorphic outgroup

   # sites variable within imitator
   imi.s <- sort(unique(c(band.s, stripe.s, morph.fixed.diff)))
   imi.n.s <- length(imi.s)

   # fixed differences between each respective imitator morph and the models
   band.fixed.same <- which(snps$band.maf %% 1 == 0 & (snps$band.maf == snps$variabilis.maf | snps$band.maf == snps$summersi.maf)) # indices for sites fixed for same allele in banded and models
   stripe.fixed.same <- which(snps$stripe.maf %% 1 == 0 & (snps$stripe.maf == snps$variabilis.maf | snps$stripe.maf == snps$summersi.maf)) # indices for sites fixed for same allele in striped and models

   band.fixed.diff <- which(snp.idx %in% sort(unique(c(band.s, model.s, band.fixed.same))) == FALSE)
   band.n.f <- length(band.fixed.diff)

   stripe.fixed.diff <- which(snp.idx %in% sort(unique(c(stripe.s, model.s, stripe.fixed.same))) == FALSE)
   stripe.n.f <- length(stripe.fixed.diff)

   # mask sites where at least one imitator morph and a model shares a fixed allele
   imi.same <- sort(unique(c(band.fixed.same, stripe.fixed.same))) # indices for sites where band is fixed and matches a model OR striped is fixed and matches a model

   # fixed differences between imitator and models
   imi.fixed.diff <- which(1:nrow(snps) %in% sort(unique(c(imi.s, model.s, imi.same))) == FALSE)
   imi.n.f <- length(imi.fixed.diff)

   # effective number of sites, i.e. sites variable within imitator OR fixed differences between imitator and models
   band.n.effective <- band.n.s + band.n.f
   stripe.n.effective <- stripe.n.s + stripe.n.f
   imi.n.effective <- imi.n.s + imi.n.f

   return(c(imi.n.s, imi.n.f, imi.n.effective, band.n.s, band.n.f, band.n.effective, stripe.n.s, stripe.n.f, stripe.n.effective))
}

gof <- function(obs.s, obs.f, exp.s, exp.f) {
   if (is.na(exp.s) || is.na(exp.f) || exp.s == 0 || exp.f == 0) return(NA)
   x = (obs.s - exp.s)^2/exp.s + (obs.f - exp.f)^2/exp.f
   return(x)
}

gof_yates <- function(obs.s, obs.f, exp.s, exp.f) {
   if (is.na(exp.s) || is.na(exp.f) || exp.s == 0 || exp.f == 0) return(NA)
   x = (abs(obs.s - exp.s) - 0.5)^2/exp.s + (abs(obs.f - exp.f) - 0.5)^2/exp.f
   return(x)
}

## extract SNPs that are segregating across banded+striped+models
cat(paste0("Extracting SNPs with p-value <= ", maxp,"\n"))

snps <- mafs[which(mafs$meta.snp.pval <= maxp),]
regions <- unique(snps$chromo) # all contigs with SNPs

# recode allele frequencies based on MAF cutoff
maxmaf = 1 - minmaf
snps$band.maf[which(snps$band.maf <= minmaf)] = 0
snps$band.maf[which(snps$band.maf >= maxmaf)] = 1
snps$stripe.maf[which(snps$stripe.maf <= minmaf)] = 0
snps$stripe.maf[which(snps$stripe.maf >= maxmaf)] = 1
snps$variabilis.maf[which(snps$variabilis.maf <= minmaf)] = 0
snps$variabilis.maf[which(snps$variabilis.maf >= maxmaf)] = 1
snps$summersi.maf[which(snps$summersi.maf <= minmaf)] = 0
snps$summersi.maf[which(snps$summersi.maf >= maxmaf)] = 1

## Calculate genome-wide rates of segregating sites and fixed differences
cat(paste0("Calculating genome-wide rates using MAF >= ", minmaf,"\n"))
genome.counts <- count(snps=snps)

# rates are [segregating, fixed_difference]
imi.rates = NULL
if (genome.counts[3] > 0) {
   imi.rates = c(genome.counts[1]/genome.counts[3], genome.counts[2]/genome.counts[3])
   cat("\nimitator.n.effective\timitator.n.seg\timitator.n.fixed\n")
   cat(paste0(genome.counts[3],"\t",genome.counts[1],"\t",genome.counts[2],"\n"))
} else cat("\nWARNING: No genome-wide effective sites for imitator-wide -> skipping HKA calculations\n\n")

band.rates = NULL
if (genome.counts[6] > 0) {
   band.rates = c(genome.counts[4]/genome.counts[6], genome.counts[5]/genome.counts[6])
   cat("\nband.n.effective\tband.n.seg\tband.n.fixed\n")
   cat(paste0(genome.counts[6],"\t",genome.counts[4],"\t",genome.counts[5],"\n"))
} else cat("\nWARNING: No genome-wide effective sites for banded -> skipping HKA calculations\n\n")

stripe.rates = NULL
if (genome.counts[9] > 0) {
   stripe.rates = c(genome.counts[7]/genome.counts[9], genome.counts[8]/genome.counts[9])
   cat("\nstripe.n.effective\tstripe.n.seg\tstripe.n.fixed\n")
   cat(paste0(genome.counts[9],"\t",genome.counts[7],"\t",genome.counts[8],"\n"))
} else cat("\nWARNING: No genome-wide effective sites for striped -> skipping HKA calculations\n\n")

## Caculate expected and observed numbers of segregating sites and polymorphisms for each region
cat("\nCalculating observed and expected counts for each region\n")

hka.df <- data.frame(matrix(NA,nrow=length(regions),ncol=19))
colnames(hka.df) <- c("region", "imi.obs.seg", "imi.obs.fix", "imi.exp.seg", "imi.exp.fix", "imi.hka", "imi.hka.yates", "band.obs.seg", "band.obs.fix", "band.exp.seg", "band.exp.fix", "band.hka", "band.hka.yates", "stripe.obs.seg", "stripe.obs.fix", "stripe.exp.seg", "stripe.exp.fix", "stripe.hka", "stripe.hka.yates")
hka.df$region = regions

for (i in 1:nrow(hka.df)) {
   c <- count(snps[which(snps$chromo == hka.df$region[i]),])
   if (!is.null(imi.rates)) {
      hka.df[i,2] = c[1] # imi obs S
      hka.df[i,3] = c[2] # imi obs fixed diff
      hka.df[i,4] = c[3] * imi.rates[1] # imi expected S
      hka.df[i,5] = c[3] * imi.rates[2] # imi expected fixed diff
      hka.df[i,6] = gof(obs.s = c[1], obs.f = c[2], exp.s = hka.df[i,4], exp.f = hka.df[i,5]) # imi GOF
      hka.df[i,7] = gof_yates(obs.s = c[1], obs.f = c[2], exp.s = hka.df[i,4], exp.f = hka.df[i,5]) # imi GOF with Yates correction
   }

   if (!is.null(band.rates)) {
      hka.df[i,8] = c[4] # band obs S
      hka.df[i,9] = c[5] # band obs fixed diff
      hka.df[i,10] = c[6] * band.rates[1] # band expected S
      hka.df[i,11] = c[6] * band.rates[2] # band expected fixed diff
      hka.df[i,12] = gof(obs.s = c[4], obs.f = c[5], exp.s = hka.df[i,10], exp.f = hka.df[i,11]) # band GOF
      hka.df[i,13] = gof_yates(obs.s = c[4], obs.f = c[5], exp.s = hka.df[i,10], exp.f = hka.df[i,11]) # band GOF with Yates correction
   }

   if (!is.null(stripe.rates)) {
      hka.df[i,14] = c[7] # stripe obs S
      hka.df[i,15] = c[8] # stripe obs fixed diff
      hka.df[i,16] = c[9] * stripe.rates[1] # stripe expected S
      hka.df[i,17] = c[9] * stripe.rates[2] # stripe expected fixed diff
      hka.df[i,18] = gof(obs.s = c[7], obs.f = c[8], exp.s = hka.df[i,16], exp.f = hka.df[i,17]) # band GOF
      hka.df[i,19] = gof_yates(obs.s = c[7], obs.f = c[8], exp.s = hka.df[i,16], exp.f = hka.df[i,17]) # band GOF with Yates correction
   }

   if (i %% 1000 == 0) cat("Completed",i, "regions\n")
}

## Calculate p-values for HKA statistic and performing Genomic Control
cat("Performing Genomic Control and calculating p-values\n")

# HKA stat is distributed according to a chi-square(df = 1)
hka.df$imi.hka.pval <- pchisq(hka.df$imi.hka, df=1, lower.tail=FALSE)
hka.df$imi.hka.yates.pval <- pchisq(hka.df$imi.hka.yates, df=1, lower.tail=FALSE)

hka.df$band.hka.pval <- pchisq(hka.df$band.hka, df=1, lower.tail=FALSE)
hka.df$band.hka.yates.pval <- pchisq(hka.df$band.hka.yates, df=1, lower.tail=FALSE)

hka.df$stripe.hka.pval <- pchisq(hka.df$stripe.hka, df=1, lower.tail=FALSE)
hka.df$stripe.hka.yates.pval <- pchisq(hka.df$stripe.hka.yates, df=1, lower.tail=FALSE)

# Estimate inflation factors and apply Genomic control (see https://www.rdocumentation.org/packages/snpStats/versions/1.22.0/topics/qq.chisq)

# imitator-wide hka inflation factors
png(paste0(out_file,"_imi_hka_qq.png"))
imi.hka.lambda = qq.chisq(x = hka.df$imi.hka, df=1, overdisp=TRUE, trim=0.5, thin=NA, slope.one=TRUE, slope.lambda=TRUE, main="Imitator HKA")[3]
invisible(dev.off())
cat(paste0("Inflation factor for imitator-wide HKA stat: ", imi.hka.lambda, "\n"))

png(paste0(out_file,"_imi_hka_yates_qq.png"))
imi.hka.yates.lambda = qq.chisq(x = hka.df$imi.hka.yates, df=1, overdisp=TRUE, trim=0.5, thin=NA, slope.one=TRUE, slope.lambda=TRUE, main="Imitator Yates HKA")[3]
invisible(dev.off())
cat(paste0("Inflation factor for imitator-wide Yates-adjusted HKA stat: ", imi.hka.yates.lambda, "\n"))

# band hka inflation factors
png(paste0(out_file,"_band_hka_qq.png"))
band.hka.lambda = qq.chisq(x = hka.df$band.hka, df=1, overdisp=TRUE, trim=0.5, thin=NA, slope.one=TRUE, slope.lambda=TRUE, main="Banded HKA")[3]
invisible(dev.off())
cat(paste0("Inflation factor for banded HKA stat: ", band.hka.lambda, "\n"))

png(paste0(out_file,"_band_hka_yates_qq.png"))
band.hka.yates.lambda = qq.chisq(x = hka.df$band.hka.yates, df=1, overdisp=TRUE, trim=0.5, thin=NA, slope.one=TRUE, slope.lambda=TRUE, main="Banded Yates HKA")[3]
invisible(dev.off())
cat(paste0("Inflation factor for banded Yates-adjusted HKA stat: ", band.hka.yates.lambda, "\n"))

# stripe hka inflation factors
png(paste0(out_file,"_stripe_hka_qq.png"))
stripe.hka.lambda = qq.chisq(x = hka.df$stripe.hka, df=1, overdisp=TRUE, trim=0.5, thin=NA, slope.one=TRUE, slope.lambda=TRUE, main="Striped HKA")[3]
invisible(dev.off())
cat(paste0("Inflation factor for striped HKA stat: ", stripe.hka.lambda, "\n"))

png(paste0(out_file,"_stripe_hka_yates_qq.png"))
stripe.hka.yates.lambda = qq.chisq(x = hka.df$stripe.hka.yates, df=1, overdisp=TRUE, trim=0.5, thin=NA, slope.one=TRUE, slope.lambda=TRUE, main="Striped Yates HKA")[3]
invisible(dev.off())
cat(paste0("Inflation factor for striped Yates-adjusted HKA stat: ", stripe.hka.yates.lambda, "\n"))

# adjust HKA statistics by inflation factors
hka.df$imi.hka.gc <- hka.df$imi.hka/imi.hka.lambda
hka.df$imi.hka.yates.gc <- hka.df$imi.hka.yates/imi.hka.yates.lambda

hka.df$band.hka.gc <- hka.df$band.hka/band.hka.lambda
hka.df$band.hka.yates.gc <- hka.df$band.hka.yates/band.hka.yates.lambda

hka.df$stripe.hka.gc <- hka.df$stripe.hka/stripe.hka.lambda
hka.df$stripe.hka.yates.gc <- hka.df$stripe.hka.yates/stripe.hka.yates.lambda

# calculated p-values for genomic-controlled HKA stats

hka.df$imi.hka.gc.pval <- pchisq(hka.df$imi.hka.gc, df=1, lower.tail=FALSE)
hka.df$imi.hka.yates.gc.pval <- pchisq(hka.df$imi.hka.yates.gc, df=1, lower.tail=FALSE)

hka.df$band.hka.gc.pval <- pchisq(hka.df$band.hka.gc, df=1, lower.tail=FALSE)
hka.df$band.hka.yates.gc.pval <- pchisq(hka.df$band.hka.yates.gc, df=1, lower.tail=FALSE)

hka.df$stripe.hka.gc.pval <- pchisq(hka.df$stripe.hka.gc, df=1, lower.tail=FALSE)
hka.df$stripe.hka.yates.gc.pval <- pchisq(hka.df$stripe.hka.yates.gc, df=1, lower.tail=FALSE)

## Output qq-plot for p-values and write results
cat("Writing results\n")

# p-value qqplots

imi.exp.p <- ppoints(length(which(!is.na(hka.df$imi.hka.gc.pval))))
png(file=paste0(out_file,"_imitator_hka_pval_qq.png"),width=800,height=800,res=120)
par(mar=c(5,5,5,4))
plot(sort(-log10(imi.exp.p)), sort(-log10(hka.df$imi.hka.gc.pval[-which(is.na(hka.df$imi.hka.gc.pval))])), xlab=expression("-log"[10]*"(expected p-value)"), ylab=expression("-log"[10]*"(observed p-value)"), main="Imitator HKA GC", cex.axis=1.25, cex.lab=1.25)
abline(0,1, col="red",lty=2)
invisible(dev.off())

png(file=paste0(out_file,"_imitator_hka_yates_pval_qq.png"),width=800,height=800,res=120)
par(mar=c(5,5,5,4))
plot(sort(-log10(imi.exp.p)), sort(-log10(hka.df$imi.hka.yates.gc.pval[-which(is.na(hka.df$imi.hka.yates.gc.pval))])), xlab=expression("-log"[10]*"(expected p-value)"), ylab=expression("-log"[10]*"(observed p-value)"), main="Imitator Yates HKA GC", cex.axis=1.25, cex.lab=1.25)
abline(0,1, col="red",lty=2)
invisible(dev.off())

band.exp.p <- ppoints(length(which(!is.na(hka.df$band.hka.gc.pval))))
png(file=paste0(out_file,"_band_hka_pval_qq.png"),width=800,height=800,res=120)
par(mar=c(5,5,5,4))
plot(sort(-log10(band.exp.p)), sort(-log10(hka.df$band.hka.gc.pval[-which(is.na(hka.df$band.hka.gc.pval))])), xlab=expression("-log"[10]*"(expected p-value)"), ylab=expression("-log"[10]*"(observed p-value)"), main="Banded HKA GC", cex.axis=1.25, cex.lab=1.25)
abline(0,1, col="red",lty=2)
invisible(dev.off())

png(file=paste0(out_file,"_band_hka_yates_pval_qq.png"),width=800,height=800,res=120)
par(mar=c(5,5,5,4))
plot(sort(-log10(band.exp.p)), sort(-log10(hka.df$band.hka.yates.gc.pval[-which(is.na(hka.df$band.hka.yates.gc.pval))])), xlab=expression("-log"[10]*"(expected p-value)"), ylab=expression("-log"[10]*"(observed p-value)"), main="Banded Yates HKA GC", cex.axis=1.25, cex.lab=1.25)
abline(0,1, col="red",lty=2)
invisible(dev.off())

stripe.exp.p <- ppoints(length(which(!is.na(hka.df$stripe.hka.gc.pval))))
png(file=paste0(out_file,"_stripe_hka_pval_qq.png"),width=800,height=800,res=120)
par(mar=c(5,5,5,4))
plot(sort(-log10(stripe.exp.p)), sort(-log10(hka.df$stripe.hka.gc.pval[-which(is.na(hka.df$stripe.hka.gc.pval))])), xlab=expression("-log"[10]*"(expected p-value)"), ylab=expression("-log"[10]*"(observed p-value)"), main="Striped HKA GC", cex.axis=1.25, cex.lab=1.25)
abline(0,1, col="red",lty=2)
invisible(dev.off())

stripe.exp.p <- ppoints(length(which(!is.na(hka.df$stripe.hka.gc.pval))))
png(file=paste0(out_file,"_stripe_hka_yates_pval_qq.png"),width=800,height=800,res=120)
par(mar=c(5,5,5,4))
plot(sort(-log10(stripe.exp.p)), sort(-log10(hka.df$stripe.hka.yates.gc.pval[-which(is.na(hka.df$stripe.hka.yates.gc.pval))])), xlab=expression("-log"[10]*"(expected p-value)"), ylab=expression("-log"[10]*"(observed p-value)"), main="Striped Yates HKA GC", cex.axis=1.25, cex.lab=1.25)
abline(0,1, col="red",lty=2)
invisible(dev.off())

# reformat dataframe for pretty output
hka.df2 <- hka.df[,c(1,2,3,4,5,6,20,7,21,26,32,27,33,8,9,10,11,12,22,13,23,28,34,29,35,14,15,16,17,18,24,19,25,30,36,31,37)]
write.table(hka.df2, file=out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
