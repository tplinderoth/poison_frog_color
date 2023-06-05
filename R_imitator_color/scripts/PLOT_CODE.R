### PLOT_CODE.R ###

## SFS
workdir='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/sfs'

all1 <- scan('all_imi_allsites_set1.sfs')
all1 <- all1[2:length(all1)]
all1.p <- all1/sum(all1)

all1s <- scan('all_imi_allsites_set1strict.sfs')
all1s <- all1s[2:length(all1s)]
all1s.p <- all1s/sum(all1s)

# plot difference in strictly filtered and more liberally filtered SFS
sfsdiff <- abs(all1.p - all1s.p)
pdf(file="coverage_filter_SFS_comparison.pdf", width=12, height=14)
#x11(width=12, height=14)
par(mfrow=c(3,1), mar=c(5, 6.5, 4, 2)+1)
#barplot(all1, xlab="MAF", ylab="Number of SNPs", main=paste0("all_imi_allsites_set1.sfs","\n",">14 individuals with >2X coverage in each population"))
barplot(all1, xlab="MAF", ylab="Number of sites", main="", cex.lab=1.8, cex.axis=1.6, names=1:length(all1), cex.names=1.5)
text(x=70,y=150000,labels="> 14 individuals from each population with > 2X depth",cex=2)
#barplot(all1s, xlab="MAF", ylab="Number of SNPs", main=paste0("all_imi_allsites_set1strict.sfs","\n",">95% of individuals with >3X coverage in each population"))
barplot(all1s, xlab="MAF", ylab="Number of sites", main="", cex.lab=1.8, cex.axis=1.6, names=1:length(all1s), cex.names=1.5)
text(x=70,y=60000,labels="> 95% of individuals from each population with > 3X depth",cex=2)
#plot(x=1:length(sfsdiff), y=sfsdiff, type="b", pch=16, xlab="MAF", ylab="Absolute difference in the proportion of SNPs", main="all_imi_allsites_set1.sfs vs all_imi_allsites_set1strict.sfs")
plot(x=1:length(sfsdiff), y=sfsdiff, type="b", pch=16, xlab="MAF", ylab="", main="", cex.lab=1.8, cex.axis=1.6)
title(ylab=paste("Absolute difference in the","\n","proportion of sites","\n",sep=""), cex.lab=1.8, line=2)
abline(v=5, col="red", lty=2)
text(x=15,y=0.03,"2% MAF cutoff",cex=1.8)
#legend('topright', legend="MAF cutoff", col="red", lty=2, bty='n')
dev.off()

band1 <- scan('banded_allsites_set1.sfs')
band1 <- band1[2:length(band1)]
stripe1 <- scan('striped_allsites_set1.sfs')
stripe1 <- stripe1[2:length(stripe1)]
admix1 <- scan('admixed_allsites_set1.sfs')
admix1 <- admix1[2:length(admix1)]

pdf(file='imi_sfs_qc_set1.pdf',width=12,height=14)
par(mfrow=c(3,1),cex=1.1)
barplot(band1,main="Banded (n=33) QC set1 SFS", xlab="MAF", ylab="Number SNPs", names.arg=c(1:length(band1)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2,col="slateblue4")
barplot(stripe1,main="Striped (n=33) QC set1 SFS", xlab="MAF", ylab="Number SNPs", names.arg=c(1:length(stripe1)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2,col="gold2")
barplot(admix1,main="Admixed (n=58) QC set1 SFS", xlab="MAF", ylab="Number SNPs", names.arg=c(1:length(admix1)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2,col="cyan4")
dev.off()

band2 <- scan('banded_allsites_set2.sfs')
band2 <- band2[2:length(band2)]
stripe2 <- scan('striped_allsites_set2.sfs')
stripe2 <- stripe2[2:length(stripe2)]
admix2 <- scan('admixed_allsites_set2.sfs')
admix2 <- admix2[2:length(admix2)]

pdf(file='imi_sfs_qc_set2.pdf',width=12,height=14)
par(mfrow=c(3,1),cex=1.1)
barplot(band2,main="Banded (n=33) QC set2 SFS", xlab="MAF", ylab="Number SNPs",names.arg=c(1:length(band2)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2,col="slateblue4")
barplot(stripe2,main="Striped (n=33) QC set2 SFS", xlab="MAF", ylab="Number SNPs",names.arg=c(1:length(stripe2)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2,col="gold2")
barplot(admix2,main="Admixed (n=58) QC set2 SFS", xlab="MAF", ylab="Number SNPs",names.arg=c(1:length(admix2)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2,col="cyan4")
dev.off()

band3 <- scan('banded_allsites_set3.sfs')
band3 <- band3[2:length(band3)]
stripe3 <- scan('striped_allsites_set3.sfs')
stripe3 <- stripe3[2:length(stripe3)]
admix3 <- scan('admixed_allsites_set3.sfs')
admix3 <- admix3[2:length(admix3)]

pdf(file='imi_sfs_qc_set3.pdf',width=12,height=14)
par(mfrow=c(3,1),cex=1.1)
barplot(band3,main="Banded (n=33) QC set3 SFS", xlab="MAF", ylab="Number SNPs",names.arg=c(1:length(band3)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2, col="slateblue4")
barplot(stripe3,main="Striped (n=33) QC set3 SFS", xlab="MAF", ylab="Number SNPs",names.arg=c(1:length(stripe3)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2,col="gold2")
barplot(admix3,main="Admixed (n=58) QC set3 SFS", xlab="MAF", ylab="Number SNPs",names.arg=c(1:length(admix3)), cex.axis=1.3, cex.lab=1.3, cex.names=1.2,col="cyan4")
dev.off()

#
#
#

## Trees
wordir='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/tree'

R
library(ggtree)
library(scales)

morphs <- read.table('band_stripe_model_morph_map.txt',head=FALSE)

# exome-wide tree
exome <- read.tree('band_stripe_models.raxml.support')
tipcol.exome <- as.character(morphs$V2[match(exome$tip.label, morphs$V1)])
tipcol.exome <- replace(tipcol.exome, which(tipcol.exome == "band"), "slateblue4")
tipcol.exome <- replace(tipcol.exome, which(tipcol.exome == "stripe"), "gold2")
tipcol.exome <- replace(tipcol.exome, which(tipcol.exome == "variabilis"), "olivedrab4")
tipcol.exome <- replace(tipcol.exome, which(tipcol.exome == "summersi"), "maroon1")
#png(file='exome_all_species_tree.png', bg="transparent")
exome_tree <- ggtree(exome, layout="equal_angle") + geom_tippoint(fill=alpha(tipcol.exome,0.4), color=tipcol.exome, shape=21, size=5, stroke=1) + theme(plot.background = ggplot2::element_rect(fill="transparent", colour=NA), panel.background = ggplot2::element_rect(fill="transparent", colour=NA))
ggsave("exome_all_species_tree.png", exome_tree, bg="transparent")
#dev.off()


# mc1r tree
mc1r <- read.tree('band_stripe_models_mc1r.nwk')
tipcol.mc1r <- as.character(morphs$V2[match(mc1r$tip.label, morphs$V1)])
tipcol.mc1r <- replace(tipcol.mc1r, which(tipcol.mc1r == "band"), "slateblue4")
tipcol.mc1r <- replace(tipcol.mc1r, which(tipcol.mc1r == "stripe"), "gold2")
tipcol.mc1r <- replace(tipcol.mc1r, which(tipcol.mc1r == "variabilis"), "olivedrab4")
tipcol.mc1r <- replace(tipcol.mc1r, which(tipcol.mc1r == "summersi"), "maroon1")
#png(file='mc1r_all_species_tree.png', bg="transparent")
mc1r_tree <- ggtree(mc1r, layout="equal_angle") + geom_tippoint(fill=alpha(tipcol.mc1r,0.4), color=tipcol.mc1r, shape=21, size=5, stroke=1) + theme(plot.background = ggplot2::element_rect(fill="transparent", colour=NA), panel.background = ggplot2::element_rect(fill="transparent", colour=NA))
ggsave("mc1r_all_species_tree.png", mc1r_tree, bg="transparent")
#dev.off()

# asip tree
read.tree('band_stripe_models_asip.nwk')
tipcol.asip <- as.character(morphs$V2[match(asip$tip.label, morphs$V1)])
tipcol.asip <- replace(tipcol.asip, which(tipcol.asip == "band"), "slateblue4")
tipcol.asip <- replace(tipcol.asip, which(tipcol.asip == "stripe"), "gold2")
tipcol.asip <- replace(tipcol.asip, which(tipcol.asip == "variabilis"), "olivedrab4")
tipcol.asip <- replace(tipcol.asip, which(tipcol.asip == "summersi"), "maroon1")
#png(file='asip_all_species_tree.png', bg="transparent")
asip_tree <- ggtree(asip, layout="equal_angle") + geom_tippoint(fill=alpha(tipcol.asip,0.4), color=tipcol.asip, shape=21, size=5, stroke=1) + theme(plot.background = ggplot2::element_rect(fill="transparent", colour=NA), panel.background = ggplot2::element_rect(fill="transparent", colour=NA))
ggsave('asip_all_species_tree.png', asip_tree, bg="transparent")
#dev.off()

# bsn tree
bsn <- read.tree('band_stripe_models_bsn.nwk')
tipcol.bsn <- as.character(morphs$V2[match(bsn$tip.label, morphs$V1)])
tipcol.bsn <- replace(tipcol.bsn, which(tipcol.bsn == "band"), "slateblue4")
tipcol.bsn <- replace(tipcol.bsn, which(tipcol.bsn == "stripe"), "gold2")
tipcol.bsn <- replace(tipcol.bsn, which(tipcol.bsn == "variabilis"), "olivedrab4")
tipcol.bsn <- replace(tipcol.bsn, which(tipcol.bsn == "summersi"), "maroon1")
#png(file='bsn_all_species_tree.png', bg="transparent")
bsn_tree <- ggtree(bsn, layout="equal_angle") + geom_tippoint(fill=alpha(tipcol.bsn,0.4), color=tipcol.bsn, shape=21, size=5, stroke=1) + theme(plot.background = ggplot2::element_rect(fill="transparent", colour=NA), panel.background = ggplot2::element_rect(fill="transparent", colour=NA))
ggsave('bsn_all_species_tree.png', bsn_tree, bg="transparent")
#dev.off()

# krt8.2 tree
krt8 <- read.tree('band_stripe_models_krt8.2.nwk')
tipcol.krt8 <- as.character(morphs$V2[match(krt8$tip.label, morphs$V1)])
tipcol.krt8 <- replace(tipcol.krt8, which(tipcol.krt8 == "band"), "slateblue4")
tipcol.krt8 <- replace(tipcol.krt8, which(tipcol.krt8 == "stripe"), "gold2")
tipcol.krt8 <- replace(tipcol.krt8, which(tipcol.krt8 == "variabilis"), "olivedrab4")
tipcol.krt8 <- replace(tipcol.krt8, which(tipcol.krt8 == "summersi"), "maroon1")
#png(file='krt8.2_all_species_tree.png', bg="transparent")
krt8_tree <- ggtree(krt8, layout="equal_angle") + geom_tippoint(fill=alpha(tipcol.krt8,0.4), color=tipcol.krt8, shape=21, size=5, stroke=1) + theme(plot.background = ggplot2::element_rect(fill="transparent", colour=NA), panel.background = ggplot2::element_rect(fill="transparent", colour=NA))
ggsave('krt8_all_species_tree.png', krt8_tree, bg="transparent")
#dev.off()

# retsat tree
retsat <- read.tree('band_stripe_models_retsat.nwk')
tipcol.retsat <- as.character(morphs$V2[match(retsat$tip.label, morphs$V1)])
tipcol.retsat <- replace(tipcol.retsat, which(tipcol.retsat == "band"), "slateblue4")
tipcol.retsat <- replace(tipcol.retsat, which(tipcol.retsat == "stripe"), "gold2")
tipcol.retsat <- replace(tipcol.retsat, which(tipcol.retsat == "variabilis"), "olivedrab4")
tipcol.retsat <- replace(tipcol.retsat, which(tipcol.retsat == "summersi"), "maroon1")
#png(file='retsat_all_species_tree.png', bg="transparent")
retsat_tree <- ggtree(retsat, layout="equal_angle") + geom_tippoint(fill=alpha(tipcol.retsat,0.4), color=tipcol.retsat, shape=21, size=5, stroke=1) + theme(plot.background = ggplot2::element_rect(fill="transparent", colour=NA), panel.background = ggplot2::element_rect(fill="transparent", colour=NA))
ggsave('retsat_all_species_tree.png', retsat_tree, bg="transparent")
#dev.off()

# export all trees on the same scale
x=range(exome_tree$data$x, mc1r_tree$data$x, asip_tree$data$x, bsn_tree$data$x, krt8_tree$data$x, retsat_tree$data$x)
y=range(exome_tree$data$y, mc1r_tree$data$y, asip_tree$data$y, bsn_tree$data$y, krt8_tree$data$y, retsat_tree$data$y)

#xlim(x[1],x[2])
#<ScaleContinuousPosition>
# Range:  
# Limits: -1.27 -- 0.244
#ylim(y[1],y[2])
#<ScaleContinuousPosition>
# Range:  
# Limits: -0.081 -- 0.251

exome_tree.scale <- exome_tree + xlim(x[1],x[2]) + ylim(y[1],y[2])
ggsave('exome_all_species_tree_scale.png', exome_tree.scale, bg="transparent")
ggsave('exome_all_species_tree_scale_bar.png', exome_tree.scale + geom_treescale(x=-0.5, y=0, width=0.1, offset=0.001), bg="transparent")
ggsave('exome_all_species_tree_scale_bar2.png', exome_tree.scale + geom_treescale(x=-0.8, y=0, width=0.5, offset=0.005, linesize=2, fontsize=10), bg="transparent")

mc1r_tree.scale <- mc1r_tree + xlim(x[1],x[2]) + ylim(y[1],y[2])
ggsave('mc1r_all_species_tree_scale.png', mc1r_tree.scale, bg="transparent")
ggsave('mc1r_all_species_tree_scale_bar.png', mc1r_tree.scale + geom_treescale(x=-0.5, y=0, width=0.1, offset=0.001), bg="transparent")

asip_tree.scale <- asip_tree + xlim(x[1],x[2]) + ylim(y[1],y[2])
ggsave('asip_all_species_tree_scale.png', asip_tree.scale, bg="transparent")
ggsave('asip_all_species_tree_scale_bar.png', asip_tree.scale + geom_treescale(x=-0.5, y=0, width=0.1, offset=0.001), bg="transparent")

bsn_tree.scale <- bsn_tree + xlim(x[1],x[2]) + ylim(y[1],y[2])
ggsave('bsn_all_species_tree_scale.png', bsn_tree.scale, bg="transparent")
ggsave('bsn_all_species_tree_scale_bar.png', bsn_tree.scale + geom_treescale(x=-0.5, y=0, width=0.1, offset=0.001), bg="transparent")

krt8_tree.scale <- krt8_tree + xlim(x[1],x[2]) + ylim(y[1],y[2])
ggsave('krt8_all_species_tree_scale.png', krt8_tree.scale, bg="transparent")
ggsave('krt8_all_species_tree_scale_bar.png', krt8_tree.scale + geom_treescale(x=-0.5, y=0, width=0.1, offset=0.001), bg="transparent")

retsat_tree.scale <- retsat_tree + xlim(x[1],x[2]) + ylim(y[1],y[2])
ggsave('retsat_all_species_tree_scale.png', retsat_tree.scale, bg="transparent")
ggsave('retsat_all_species_tree_scale_bar.png', retsat_tree.scale + geom_treescale(x=-0.5, y=0, width=0.1, offset=0.001), bg="transparent")

#
#
#

## PCA
library(scales)

workdir='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/pca'

# define population colors
band.col="slateblue4"
stripe.col="gold2"
admix.col="cyan4"
summersi.col="slategray3"
variabilis.col="lightcoral"

# imitator and model species

covarmat <- as.matrix(read.table('imi_model_minmaf2_set2.covar',head=FALSE))
eig <- eigen(covarmat)
pcvar <- (eig$values/sum(eig$values))*100 # percent genetic variance explained by each PC
popcol <- c(rep(band.col,33),rep(stripe.col,33),rep(admix.col,58),rep(variabilis.col,3),rep(summersi.col,2)) # define colors for plotting populations
poppch.outline <- c(rep(1,124), rep(0,3), rep(2,2))
poppch.fill <- c(rep(16,124), rep(15,3), rep(17,3))

# PC1 vs PC2
#pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/imi_model_PC_1_2.pdf',width=7.5,height=7.5) # original
#pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/imi_model_PC_1_2_rev1.pdf')
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/imi_model_PC_1_2_rev1.png',res=150,width=7,height=7,units='in')
par(mar=c(5.1,5.1,4.1,2.1),mgp=c(3.5,1.3,0))
plot(x=eig$vectors[,1], y=eig$vectors[,2], xlab=paste0("PC1 (",sprintf("%.2f",pcvar[1]), "%)"), ylab=paste0("PC2 (",sprintf("%.2f",pcvar[2]),"%)"), main="", pch=poppch.outline, col=alpha(popcol,0.9), cex=2.2, cex.lab=1.9, cex.axis=1.7, lwd=2)
points(x=eig$vectors[,1], y=eig$vectors[,2]
, col=alpha(popcol,0.8), pch=poppch.fill, cex=2.2)
# uncomment two lines below for legend
legend(x=-0.04,y=0.64,legend=c("band","stripe","admix"),col="black",pt.bg=c(band.col,stripe.col,admix.col),bty='n',pt.cex=1.8,pch=21,pt.lwd=1,cex=1.5)
legend(x=0.1,y=0.64,legend=c(expression(italic("R. variabilis")), expression(italic("R. summersi"))),col="black",pt.bg=c(variabilis.col, summersi.col),bty='n',pt.cex=1.8,pch=c(22,24),pt.lwd=1,cex=1.5)
text(x=0.03,y=0.64,substitute(paste(italic('R. imitator'))), cex=1.5)
rect(xleft=-0.03, ybottom=0.61, xright=0.09, ytop=0.61,lwd=1.3)
text(x=0.205,y=0.64,'model species', cex=1.5)
rect(xleft=0.11, ybottom=0.61, xright=0.3, ytop=0.61,lwd=1.3)
dev.off()

# multipanel plot of the first 8 PCs
pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/imi_model_multipanel_PCA.pdf')
par(mfrow=c(2,2))
plot(x=eig$vectors[,1], y=eig$vectors[,2], xlab=paste0("PC1 (",sprintf("%.2f",pcvar[1]), "%)"), ylab=paste0("PC2 (",sprintf("%.2f",pcvar[2]),"%)"), main="", pch=poppch.outline, col=alpha(popcol,0.8), cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(x=eig$vectors[,1], y=eig$vectors[,2], col=alpha(popcol,0.6), pch=poppch.fill, cex=1.3)
legend(x=-0.04,y=0.7,legend=c("band","stripe","admix"),col="black",pt.bg=c(band.col,stripe.col,admix.col),bty='n',pt.cex=1.1,pch=21,pt.lwd=1,cex=1.1)
legend(x=0.1,y=0.7,legend=c(expression(italic("R. variabilis")), expression(italic("R. summersi"))),col="black",pt.bg=c(variabilis.col, summersi.col),bty='n',pt.cex=1.1,pch=c(22,24),pt.lwd=1,cex=1.1)

# PC3 vs PC4
plot(x=eig$vectors[,3], y=eig$vectors[,4], xlab=paste0("PC3 (",sprintf("%.2f",pcvar[3]), "%)"), ylab=paste0("PC4 (",sprintf("%.2f",pcvar[4]),"%)"), main="", pch=poppch.outline, col=alpha(popcol,0.8), cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(x=eig$vectors[,3], y=eig$vectors[,4], col=alpha(popcol,0.6), pch=poppch.fill, cex=1.3)

# PC5 vs PC6
plot(x=eig$vectors[,5], y=eig$vectors[,6], xlab=paste0("PC5 (",sprintf("%.2f",pcvar[5]), "%)"), ylab=paste0("PC6 (",sprintf("%.2f",pcvar[6]),"%)"), main="", pch=poppch.outline, col=alpha(popcol,0.8), cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(x=eig$vectors[,5], y=eig$vectors[,6], col=alpha(popcol,0.6), pch=poppch.fill, cex=1.3)

# PC7 vs PC8
plot(x=eig$vectors[,7], y=eig$vectors[,8], xlab=paste0("PC7 (",sprintf("%.2f",pcvar[7]), "%)"), ylab=paste0("PC8 (",sprintf("%.2f",pcvar[8]),"%)"), main="", pch=poppch.outline, col=alpha(popcol,0.8), cex=1.3, cex.lab=1.3, cex.axis=1.3)
points(x=eig$vectors[,7], y=eig$vectors[,8], col=alpha(popcol,0.6), pch=poppch.fill, cex=1.3)

dev.off()

# pca imitator only

covarmat.imi <- as.matrix(read.table('imi_minmaf2_set2.covar',head=FALSE))
eig.imi <- eigen(covarmat.imi)
pcvar.imi <- (eig.imi$values/sum(eig.imi$values))*100 # percent genetic variance explained by each PC
popcol.imi <- c(rep(band.col,33),rep(stripe.col,33),rep(admix.col,58)) # define colors for plotting populations

# PC1 vs PC2
#pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/imi_PC_1_2.pdf',width=7.5,height=7.5)
#pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/imi_PC_1_2_rev1.pdf')
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/imi_PC_1_2_rev1.png',width=7,height=7,units='in',res=150)
par(mar=c(5.1,5.1,4.1,2.1),mgp=c(3.5,1.3,0))
plot(x=eig.imi$vectors[,1], y=eig.imi$vectors[,2], xlab=paste0("PC1 (",sprintf("%.2f",pcvar.imi[1]), "%)"), ylab=paste0("PC2 (",sprintf("%.2f",pcvar.imi[2]),"%)"), main="", pch=1, col=alpha(popcol,0.8), cex=2, cex.lab=1.9, cex.axis=1.7)
points(x=eig.imi$vectors[,1], y=eig.imi$vectors[,2], col=alpha(popcol,0.6), pch=16, cex=2)
text(x=-0.1,y=0.15,expression(italic("R. imitator")),cex=1.5)
dev.off()

# pca imitator only with localities

library(ggplot2)

pc.df <- read.table('~/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/pca/imi_minmaf2_set2_PC.txt',head=TRUE,sep="\t")
pc.df$rivcol <- replace(as.character(pc.df$river), which(as.character(pc.df$river) == "North"), "black")
pc.df$rivcol <- replace(pc.df$rivcol, which(pc.df$rivcol == "South"), "red")

# define population colors
popcols = c("#473C8B99", "#EEC90099", "#008B8B99", "#008B8B99") # band, stripe, admix_north, admix_south

# calculate convex hulls
pc.band <- pc.df[which(pc.df$pop == "band"),]
hulldat <- pc.band[chull(x=pc.band$PC1, y=pc.band$PC2),]

pc.stripe <- pc.df[which(pc.df$pop == "stripe"),]
hulldat <- rbind(hulldat, pc.stripe[chull(x=pc.stripe$PC1, y=pc.stripe$PC2),])

pc.admix.north <- pc.df[which(pc.df$pop == "admix" & pc.df$river == "North"),]
pc.admix.north$pop = "admix_north"
hulldat <- rbind(hulldat, pc.admix.north[chull(x=pc.admix.north$PC1, y=pc.admix.north$PC2),])

pc.admix.south <- pc.df[which(pc.df$pop == "admix" & pc.df$river == "South"),]
pc.admix.south$pop = "admix_south"
hulldat <- rbind(hulldat, pc.admix.south[chull(x=pc.admix.south$PC1, y=pc.admix.south$PC2),])

#plot
pcplot <- ggplot(aes(x = PC1, y = PC2), data=pc.df) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size=15), axis.title =
element_text(size=17)) + scale_y_continuous(breaks=seq(-0.2,0.2,by=0.05))
# add population polygon
pcplot <- pcplot + geom_polygon(data=hulldat, aes(group=pop, fill=pop, colour=pop, alpha=0.3), show.legend=FALSE) + scale_fill_manual(values=popcols) + scale_color_manual(values=popcols)
# add points
pcplot <- pcplot + geom_point(aes(x = PC1, y = PC2), color=pc.df$rivcol, shape=1, lwd=2, stroke=1)
# modify axes
pcplot + xlab("PC1 (6.15%)") + ylab("PC2 (1.97%)")

# save plot
ggsave(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/pca/imi_minmaf2_set2_sampling_pca.png')

#
#
#

## Admixture

# K2 nice plot

workdir='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/admixture'
band.col="slateblue4"
stripe.col="gold2"

admixp <- read.table('band_stripe_admix_minmaf2_set2_k2.qopt',head=FALSE)

band <- admixp[1:33,] #band
band <- band[order(band$V1),]

stripe <- admixp[34:66,] #stripe
stripe <- stripe[order(stripe$V1),]

admix <- admixp[67:124,] # admix
admix <- admix[order(admix$V1),]

admixp.sort <- rbind(band,admix,stripe) # rearrange for banded, admixed, stripe order
admixp.mat <- t(as.matrix(admixp.sort))
pops <- c(rep("band",33),rep("admix",58),rep("stripe",33))

#x11(width=20, height=6)
#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/band_stripe_admix_minmaf2_set2_k2.png',width=1200,height=360)
#pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/band_stripe_admix_minmaf2_set2_k2.pdf',width=20,height=6) # original
pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/band_stripe_admix_minmaf2_set2_k2_rev1.pdf',width=20,height=3.5) # science revised
par(mar=c(5.1,5.7,4.1,2.1),mgp=c(3.8,1,0))
#barplot(admixp.mat,col=c(stripe.col,band.col),space=0,xlab="Individuals",ylab="Ancestry proportion", axisnames=FALSE, cex.lab=2.3, cex.axis=2.2) # original plot
barplot(admixp.mat,col=c(stripe.col,band.col),space=0,xlab="Individuals",ylab="Ancestry", axisnames=FALSE, cex.lab=2.3, cex.axis=2.2) # revised Science plot
lines(x<-c(33,33),y<-c(0,1),col="grey100",lwd=2, lty=1)
lines(x<-c(53,53),y<-c(0,1),col="grey100",lwd=2, lty=1)
lines(x<-c(91,91),y<-c(0,1),col="grey100",lwd=2, lty=1)
#text(c(16.5,62,107.5), -0.07, unique(pops),xpd=T,cex=2.3) # used in original plot
dev.off()

# K 3-5

admixp.3 <- read.table('band_stripe_admix_minmaf2_set2_k3.qopt',head=FALSE)
admixp.4 <- read.table('band_stripe_admix_minmaf2_set2_k4.qopt',head=FALSE)
admixp.5 <- read.table('band_stripe_admix_minmaf2_set2_k5.qopt',head=FALSE)

# band
band.q3 <- admixp.3[1:33,]
band.q4 <- admixp.4[1:33,]
band.q5 <- admixp.5[1:33,]

# stripe
stripe.q3 <- admixp.3[34:66,]
stripe.q4 <- admixp.4[34:66,]
stripe.q5 <- admixp.5[34:66,]

# admixed - sort by North and South of Huallaga river
meta <- read.csv('/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/admixture/band_stripe_admix_minmaf2_set2_k2_admixp.txt', head=TRUE, sep="\t")
meta.admix <- meta[which(meta$pop == "admix"),1:7]
geoidx <- c(which(meta.admix$river == "South"), which(meta.admix$river == "North"))

admix.q3 <- admixp.3[67:124,]
admix.q3 <- admix.q3[geoidx,]
admix.q4 <- admixp.4[67:124,]
admix.q4 <- admix.q4[geoidx,]
admix.q5 <- admixp.5[67:124,]
admix.q5 <- admix.q5[geoidx,]

# assemble matrix of where each population has been sorted by admixture proportion
# order by banded, admixed, striped

admixp3.sort <- rbind(band.q3,admix.q3,stripe.q3)
admixp3.mat <- t(as.matrix(admixp3.sort))

admixp4.sort <- rbind(band.q4,admix.q4,stripe.q4)
admixp4.mat <- t(as.matrix(admixp4.sort))

admixp5.sort <- rbind(band.q5,admix.q5,stripe.q5)
admixp5.mat <- t(as.matrix(admixp5.sort))

# plot k3-5

col1="gold2"
col2="slateblue4"
col3="cyan4"
col4="deeppink4"
col5="rosybrown"

# k=3
#x11(width=20, height=6)
pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/admixture/plots/band_stripe_admix_minmaf2_set2_k3.pdf',width=20,height=6)
par(mar=c(5.1,5.7,4.1,2.1),mgp=c(2,-1,-2))
barplot(admixp3.mat,col=c(col2,col3,col1),space=0,xlab="",ylab="Ancestry proportion", axisnames=FALSE, cex.lab=2.3, cex.axis=2.2)
lines(x<-c(33,33),y<-c(0,1),col="grey100",lwd=2, lty=1)
lines(x<-c(91,91),y<-c(0,1),col="grey100",lwd=2, lty=1)
lines(x<-c(53,53),y<-c(0,1),col="grey100",lwd=2, lty=1)
text(c(16.5,43,72,107.5), 1.07, c("band","admix south", "admix north","stripe"),xpd=T,cex=2.3)
text(62.5, -0.125, "Individuals",xpd=T,cex=2.3)
text(128,0.5,"K = 3",xpd=T,cex=2.3)
dev.off()

# k=4
#x11(width=20, height=6)
pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/admixture/plots/band_stripe_admix_minmaf2_set2_k4.pdf',width=20,height=6)
par(mar=c(5.1,5.7,4.1,2.1),mgp=c(2,-1,-2))
barplot(admixp4.mat,col=c(col2,col3,col4,col1),space=0,xlab="",ylab="Ancestry proportion", axisnames=FALSE, cex.lab=2.3, cex.axis=2.2)
lines(x<-c(33,33),y<-c(0,1),col="grey100",lwd=2, lty=1)
lines(x<-c(91,91),y<-c(0,1),col="grey100",lwd=2, lty=1)
lines(x<-c(53,53),y<-c(0,1),col="grey100",lwd=2, lty=1)
text(c(16.5,43,72,107.5), 1.07, c("band","admix south", "admix north","stripe"),xpd=T,cex=2.3)
text(62.5, -0.125, "Individuals",xpd=T,cex=2.3)
text(128,0.5,"K = 4",xpd=T,cex=2.3)
dev.off()

# k=5
#x11(width=20, height=6)
pdf(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/admixture/plots/band_stripe_admix_minmaf2_set2_k5.pdf',width=20,height=6)
par(mar=c(5.1,5.7,4.1,2.1),mgp=c(2,-1,-2))
barplot(admixp5.mat,col=c(col3,col2,col1,col5,col4),space=0,xlab="",ylab="Ancestry proportion", axisnames=FALSE, cex.lab=2.3, cex.axis=2.2)
lines(x<-c(33,33),y<-c(0,1),col="grey100",lwd=2, lty=1)
lines(x<-c(91,91),y<-c(0,1),col="grey100",lwd=2, lty=1)
lines(x<-c(53,53),y<-c(0,1),col="grey100",lwd=2, lty=1)
text(c(16.5,43,72,107.5), 1.07, c("band","admix south", "admix north","stripe"),xpd=T,cex=2.3)
text(62.5, -0.125, "Individuals",xpd=T,cex=2.3)
text(128,0.5,"K = 5",xpd=T,cex=2.3)
dev.off()

#
#
#

## GWAS Manhattan

library(scales)

chrcol <- function(chrvec=NULL) {
	# chrvec: vector of contig names

	#col1="gray35" # made all light grey for revision - uncomment to get different shades for alternating contigs
	col1="gray65"
	col2="gray65"

	colvec=rep(NA,length(chrvec))
	lastcontig=chrvec[1]
	lastcol=col1
	for (i in 1:length(chrvec)) {
		if (chrvec[i] != lastcontig) {
			lastcol = ifelse(lastcol == col1, col2, col1)
			lastcontig = chrvec[i]
		}
		colvec[i] = lastcol
	}
	return(colvec)
}

## === COMBINED DIVERGENCE/ADMIXTURE === ##

mc1r_col = "darkturquoise"
asip_col = "tomato1"
bsn_col = "olivedrab2"
retsat_col = "plum"
krt8_col = "dodgerblue"
#cex.cand = 1.8
cex.cand = 2 # revised size

## PATTERN

pat <- read.table('/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/combined/pattern_combined.txt',head=TRUE)
# exclude mitochondrial sites and the following arg2 contig because assembly/mapping appears unreliable based on genome assembly
pat <- pat[-which(pat$chr == "contig8307_combined" | pat$chr == "contig13298_combined_CYTB"),]

# check for candidates among the top 20
pat20 = head(pat[order(-pat$X2),], n=20)
#write.table(pat20,'/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/combined/pattern_top20.txt',col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# The following candidate gene SNPs are in the top 20:
# contig13293_combined_mc1r:1662
# contig13293_combined_mc1r:1813
# contig6404_combined_asip:1582
# contig6404_combined_asip:219
# contig12059_combined_bsn:1764

# log transform p-values
pat$logp <- -log10(pat$X2.pval.gc)

# assign point colors and sizes to candidate SNPs
pat$pntcol <- chrcol(as.character(pat$chr))

# assign colors to candidate gene SNPs
pat$pntcol[which(pat$chr == "contig13293_combined_mc1r")] = mc1r_col
pat$pntcol[which(pat$chr == "contig6404_combined_asip")] = asip_col
pat$pntcol[which(pat$chr == "contig12059_combined_bsn")] = bsn_col

# assign colors to specific candidate SNPs
pat$pntcol[which(pat$chr == "contig13293_combined_mc1r" & pat$pos == 1662)] = mc1r_col
pat$pntcol[which(pat$chr == "contig13293_combined_mc1r" & pat$pos == 1813)] = mc1r_col
pat$pntcol[which(pat$chr == "contig6404_combined_asip" & pat$pos == 1582)] = asip_col
pat$pntcol[which(pat$chr == "contig6404_combined_asip" & pat$pos == 219)] = asip_col
pat$pntcol[which(pat$chr == "contig12059_combined_bsn" & pat$pos == 1764)] = bsn_col

pat$pntsz = 1
pat$pntsz[which(pat$chr == "contig13293_combined_mc1r" & pat$pos == 1662)] = cex.cand
pat$pntsz[which(pat$chr == "contig13293_combined_mc1r" & pat$pos == 1813)] = cex.cand
pat$pntsz[which(pat$chr == "contig6404_combined_asip" & pat$pos == 1582)] = cex.cand
pat$pntsz[which(pat$chr == "contig6404_combined_asip" & pat$pos == 219)] = cex.cand
pat$pntsz[which(pat$chr == "contig12059_combined_bsn" & pat$pos == 1764)] = cex.cand

# make dark outline circles for candidate SNPs
snp.x = which(pat$chr == "contig13293_combined_mc1r" & pat$pos == 1662 | pat$chr == "contig13293_combined_mc1r" & pat$pos == 1813 | pat$chr == "contig6404_combined_asip" & pat$pos == 1582 | 
pat$chr == "contig6404_combined_asip" & pat$pos == 219 | pat$chr == "contig12059_combined_bsn" & pat$pos == 1764)

snp.y=pat$logp[snp.x]

pat.top20line = -log10((pat$X2.pval.gc[order(pat$X2.pval.gc)])[21])
pat.bfline = -log10(0.05/nrow(pat))
#pat.ylim = c(0,max(-log10(pat$X2.pval.gc),pat.bfline)) # full plot
pat.ylim = c(1.5,max(pat$logp,pat.bfline)) # trimmed plot
#pat.ylim = c(1.5,6) # trimmed plot 2, max(pat$logp) = 5.247534

#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/pattern_manhattan.png',width=1200,height=600) #original
#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/pattern_manhattan_trim.png',width=1200,height=450) #science revision
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/pattern_manhattan_trim_noline.png',width=1200,height=450) #science revision2
par(mar=c(5,5.2,6,5))
plot(y=pat$logp, x=1:nrow(pat), pch=21, col=alpha("black",0.8), bg=pat$pntcol, cex=pat$pntsz, lwd=1, 
xlab="Genomic position", ylab=expression(-log[10](p-value)), xaxt='n', cex.axis=2, cex.lab=2, ylim=pat.ylim)
points(x=snp.x, y=snp.y, pch=1, cex=cex.cand+0.5, lwd=1.5)
#abline(h=pat.top20line, col=alpha("red3",0.4), lty=1, lwd=1.5)
#abline(h=pat.bfline, col=alpha("red3",0.4), lty=2, lwd=1.5)
dev.off()

## DORSAL COLOR

dorsum <- read.table('/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/combined/dorsalA_combined.txt', head=TRUE)
# exclude mitochondrial sites and the following arg2 contig because assembly/mapping appears unreliable based on genome assembly
dorsum <- dorsum[-which(dorsum$chr == "contig8307_combined" | dorsum$chr == "contig13298_combined_CYTB"),]

# check for candidates among the top 20
dorsum20 = head(dorsum[order(-dorsum$X2),], n=20)
#write.table(dorsum20,'/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/combined/dorsum_top20.txt',col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# The following candidate gene SNPs are in the top 20:
# contig13293_combined_mc1r:1662
# contig13293_combined_mc1r:1813
# contig6404_combined_asip:219
# contig12059_combined_bsn:6952
# contig12398_combined_retsat:6779

# log transform p-values
dorsum$logp <- -log10(dorsum$X2.pval.gc)

# assign colors and point sizes
dorsum$pntcol <- chrcol(as.character(dorsum$chr))

# assign colors to candidate gene SNPs
dorsum$pntcol[which(dorsum$chr == "contig13293_combined_mc1r")] = mc1r_col
dorsum$pntcol[which(dorsum$chr == "contig6404_combined_asip")] = asip_col
dorsum$pntcol[which(dorsum$chr == "contig12059_combined_bsn")] = bsn_col
dorsum$pntcol[which(dorsum$chr == "contig12398_combined_retsat")] = retsat_col

# assign point colors to specific candidate gene SNPs
dorsum$pntcol[which(dorsum$chr == "contig13293_combined_mc1r" & dorsum$pos == 1662)] = mc1r_col
dorsum$pntcol[which(dorsum$chr == "contig13293_combined_mc1r" & dorsum$pos == 1813)] = mc1r_col
dorsum$pntcol[which(dorsum$chr == "contig6404_combined_asip" & dorsum$pos == 219)] = asip_col
dorsum$pntcol[which(dorsum$chr == "contig12059_combined_bsn" & dorsum$pos == 6952)] = bsn_col
dorsum$pntcol[which(dorsum$chr == "contig12398_combined_retsat" & dorsum$pos == 6779)] = retsat_col

dorsum$pntsz = 1
dorsum$pntsz[which(dorsum$chr == "contig13293_combined_mc1r" & dorsum$pos == 1662)] = cex.cand
dorsum$pntsz[which(dorsum$chr == "contig13293_combined_mc1r" & dorsum$pos == 1813)] = cex.cand
dorsum$pntsz[which(dorsum$chr == "contig6404_combined_asip" & dorsum$pos == 219)] = cex.cand
dorsum$pntsz[which(dorsum$chr == "contig12059_combined_bsn" & dorsum$pos == 6952)] = cex.cand
dorsum$pntsz[which(dorsum$chr == "contig12398_combined_retsat" & dorsum$pos == 6779)] = cex.cand

# make dark outline circles for candidate SNPs
snp.dx = which(dorsum$chr == "contig13293_combined_mc1r" & dorsum$pos == 1662 | dorsum$chr == "contig13293_combined_mc1r" & dorsum$pos == 1813 | 
dorsum$chr == "contig6404_combined_asip" & dorsum$pos == 219 | dorsum$chr == "contig12059_combined_bsn" & dorsum$pos == 6952 | dorsum$chr == "contig12398_combined_retsat" & dorsum$pos == 6779)

snp.dy=dorsum$logp[snp.dx]

dorsum.top20line = -log10((dorsum$X2.pval.gc[order(dorsum$X2.pval.gc)])[21])
dorsum.bfline = -log10(0.05/nrow(dorsum))
#dorsum.ylim = c(0,max(-log10(dorsum$X2.pval.gc),dorsum.bfline))
dorsum.ylim = c(1.5,max(dorsum$logp,pat.bfline)) # trimmed plot

#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/dorsal_color_manhattan.png',width=1200,heigh=600)
#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/dorsal_manhattan_trim.png',width=1200,height=450) #science revision
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/dorsal_manhattan_trim_noline.png',width=1200,height=450) #science revision2
par(mar=c(5,5.2,6,5))
plot(y=-log10(dorsum$X2.pval.gc), x=1:nrow(dorsum), pch=21, col=alpha("black",0.8), bg=dorsum$pntcol, cex=dorsum$pntsz, lwd=1,
xlab="Genomic position", ylab=expression(-log[10](p-value)), xaxt='n', cex.axis=2, cex.lab=2, ylim=dorsum.ylim)
points(x=snp.dx, y=snp.dy, pch=1, cex=cex.cand+0.5, lwd=1.5)
#abline(h=dorsum.top20line, col=alpha("red3",0.4), lty=1, lwd=1.5)
#abline(h=dorsum.bfline, col=alpha("red3",0.4), lty=2, lwd=1.5)
dev.off()

## HINDLIMB COLOR

leg <- read.table('/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/combined/legA_combined.txt', head=TRUE)
# exclude mitochondrial sites and the following arg2 contig because assembly/mapping appears unreliable based on genome assembly
leg <- leg[-which(leg$chr == "contig8307_combined" | leg$chr == "contig13298_combined_CYTB"),]

# check for candidates among the top 20
leg20 = head(leg[order(-leg$X2),], n=20)
#write.table(leg20,'/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/combined/leg_top20.txt',col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# The following candidate gene SNPs are in the top 20:
# contig13293_combined_mc1r:1662
# contig13293_combined_mc1r:1813
# contig6404_combined_asip:219
# contig11350_combined_krt8.2:2306

# log transform p-values
leg$logp <- -log10(leg$X2.pval.gc)

# assign point colors and sizes
leg$pntcol <- chrcol(as.character(leg$chr))

# assign colors to candidate gene SNPs
leg$pntcol[which(leg$chr == "contig13293_combined_mc1r")] = mc1r_col
leg$pntcol[which(leg$chr == "contig6404_combined_asip")] = asip_col
leg$pntcol[which(leg$chr == "contig11350_combined_krt8.2")] = krt8_col

# assign colors to specific candidate SNPs
leg$pntcol[which(leg$chr == "contig13293_combined_mc1r" & leg$pos == 1662)] = mc1r_col
leg$pntcol[which(leg$chr == "contig13293_combined_mc1r" & leg$pos == 1813)] = mc1r_col
leg$pntcol[which(leg$chr == "contig6404_combined_asip" & leg$pos == 219)] = asip_col
leg$pntcol[which(leg$chr == "contig11350_combined_krt8.2" & leg$pos == 2306)] = krt8_col

leg$pntsz = 1
leg$pntsz[which(leg$chr == "contig13293_combined_mc1r" & leg$pos == 1662)] = cex.cand
leg$pntsz[which(leg$chr == "contig13293_combined_mc1r" & leg$pos == 1813)] = cex.cand
leg$pntsz[which(leg$chr == "contig6404_combined_asip" & leg$pos == 219)] = cex.cand
leg$pntsz[which(leg$chr == "contig11350_combined_krt8.2" & leg$pos == 2306)] = cex.cand

# make dark outline circles for candidate SNPs
snp.lx = which(leg$chr == "contig13293_combined_mc1r" & leg$pos == 1662 | leg$chr == "contig13293_combined_mc1r" & leg$pos == 1813 | leg$chr == "contig6404_combined_asip" & leg$pos == 219 |
leg$chr == "contig11350_combined_krt8.2" & leg$pos == 2306)

snp.ly=leg$logp[snp.lx]

leg.top20line = -log10((leg$X2.pval.gc[order(leg$X2.pval.gc)])[21])
leg.bfline = -log10(0.05/nrow(leg))
#leg.ylim = c(0,max(-log10(leg$X2.pval.gc),leg.bfline))
leg.ylim = c(1.5,max(leg$logp,pat.bfline)) # trimmed plot

#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/leg_color_manhattan.png',width=1200,heigh=600)
#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/leg_manhattan_trim.png',width=1200,height=450) #science revision
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/manhattan/leg_manhattan_trim_noline.png',width=1200,height=450) #science revision2
par(mar=c(5,5.2,6,5))
plot(y=-log10(leg$X2.pval.gc), x=1:nrow(leg), pch=21, col=alpha("black",0.8), bg=leg$pntcol, cex=leg$pntsz, lwd=1,
xlab="Genomic position", ylab=expression(-log[10](p-value)), xaxt='n', cex.axis=2, cex.lab=2, ylim=leg.ylim)
points(x=snp.lx, y=snp.ly, pch=1, cex=cex.cand+0.5, lwd=1.5)
#abline(h=leg.top20line, col=alpha("red3",0.4), lty=1, lwd=1.5)
#abline(h=leg.bfline, col=alpha("red3",0.4), lty=2, lwd=1.5)
dev.off()

#
#
#

### PHENOTYPE DISTRIBUTIONS & ALLELE FREQUENCY PIE CHARTS

library(ggplot2)
library(cowplot)
library(ggpubr)

### PHENOTYPE VIOLIN PLOTS

pheno <- read.table('/home/tyler/Dropbox/research/imitator/admixture_mapping/image_analysis/imitator_photo_phenotypes.txt',head=TRUE)
pheno$pop = factor(pheno$pop, levels=c("band","admix","stripe"),order=TRUE)

# plot2 colors
#fillcol = c("slateblue4", "maroon","gold2")
#bordcol = c("purple4", "red4", "yellow4")

# plot1 colors
fillcol = c("slateblue4", "cyan4","gold2")
bordcol = c("purple4", "darkslategray", "yellow4")

text.sz1=60
text.sz2=40
panelmar=c(1,3,0.5,0.5)

# rotation violin plot
rot <- ggplot(pheno, aes(x=pop, y=rotation)) + geom_violin(trim=FALSE, aes(fill=pop, color=pop)) + scale_fill_manual(values=fillcol) + scale_color_manual(values=bordcol)
# inlay histogram
rot <- rot + geom_boxplot(fill="white", alpha=0.6, width=0.1, color="black")
# make them transparent
#rot <- rot + theme_transparent() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=text.sz1, color="black"), axis.title.y = element_text(size=text.sz1, color="blacK"), legend.position = "none", plot.margin=unit(panelmar, "lines"), axis.line.y.left = element_line(linetype=1, color="black", size=1), axis.line.x.bottom = element_line(linetype=1, color="black", size=1)) + ylab("Pattern rotation") + xlab("Population") # 2023-03-19, changed Y-axis title from "Pattern rotation" to "Pattern orientation"

rot <- rot + theme_transparent() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=text.sz1, color="black"), axis.title.y = element_text(size=text.sz1, color="blacK"), legend.position = "none", plot.margin=unit(panelmar, "lines"), axis.line.y.left = element_line(linetype=1, color="black", linewidth=1), axis.line.x.bottom = element_line(linetype=1, color="black", linewidth=1)) + ylab("Pattern orientation") + xlab("Population")

# dorsal a-axis violin plot
dorsala <- ggplot(pheno, aes(x=pop, y=body.a)) + geom_violin(trim=FALSE, aes(fill=pop, color=pop)) + scale_fill_manual(values=fillcol) + scale_color_manual(values=bordcol)
# inlay histogram
dorsala <- dorsala + geom_boxplot(fill="white", alpha=0.6, width=0.1, color="black")
# make them transparent
#dorsala <- dorsala + theme_transparent() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=text.sz1, color="black"), axis.title.y = element_text(size=text.sz1, color="blacK"), legend.position = "none", plot.margin=unit(panelmar, "lines"), axis.line.y.left = element_line(linetype=1, color="black", size=1), axis.line.x.bottom = element_line(linetype=1, color="black", size=1)) + ylab("Dorsal body color") +xlab("Population")

dorsala <- dorsala + theme_transparent() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=text.sz1, color="black"), axis.title.y = element_text(size=text.sz1, color="blacK"), legend.position = "none", plot.margin=unit(panelmar, "lines"), axis.line.y.left = element_line(linetype=1, color="black", linewidth=1), axis.line.x.bottom = element_line(linetype=1, color="black", linewidth=1)) + ylab("Dorsal body color") +xlab("Population")

# leg a-axis violin plot
lega <- ggplot(pheno, aes(x=pop, y=leg.a)) + geom_violin(trim=FALSE, aes(fill=pop, color=pop)) + scale_fill_manual(values=fillcol) + scale_color_manual(values=bordcol)
# inlay histogram
lega <- lega + geom_boxplot(fill="white", alpha=0.6, width=0.09, color="black")
# make them transparent
#lega <- lega + theme_transparent() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=text.sz1, color="black"), axis.title.y = element_text(size=text.sz1, color="blacK"), legend.position = "none", plot.margin=unit(panelmar, "lines"), axis.line.y.left = element_line(linetype=1, color="black", size=1), axis.line.x.bottom = element_line(linetype=1, color="black",size=1)) + ylab("Hindlimb color") +xlab("Population")

lega <- lega + theme_transparent() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=text.sz1, color="black"), axis.title.y = element_text(size=text.sz1, color="blacK"), legend.position = "none", plot.margin=unit(panelmar, "lines"), axis.line.y.left = element_line(linetype=1, color="black", linewidth=1), axis.line.x.bottom = element_line(linetype=1, color="black",linewidth=1)) + ylab("Hindlimb color") +xlab("Population")


# make 3-panel plot
#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_17092020/paper_figures/pheno_box1.png',width=2000,height=600, bg="transparent")
#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_17092020/paper_figures/pheno_box2.png',width=2000,height=600, bg="transparent")
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pheno_box3.png',width=2000,height=600, bg="transparent")
plot_grid(plotlist=list(rot,dorsala,lega),nrow=1,ncol=3, rel_widths=c(1,1,1))
#dev.off()

### T-TEST AND MANN-WHITNEY U TEST FOR PHENOTYPIC DIFFERENCES BETWEEN POPS

# pattern: band vs stripe
t.test(x=pheno$rotation[pheno$pop == "band"], y=pheno$rotation[pheno$pop == "stripe"], alternative="two.sided")
wilcox.test(x=pheno$rotation[pheno$pop == "band"], y=pheno$rotation[pheno$pop == "stripe"], alternative="two.sided")

# pattern: band vs admix
t.test(x=pheno$rotation[pheno$pop == "band"], y=pheno$rotation[pheno$pop == "admix"], alternative="two.sided")
wilcox.test(x=pheno$rotation[pheno$pop == "band"], y=pheno$rotation[pheno$pop == "admix"], alternative="two.sided")

# pattern: stripe vs admix
t.test(x=pheno$rotation[pheno$pop == "stripe"], y=pheno$rotation[pheno$pop == "admix"], alternative="two.sided")
wilcox.test(x=pheno$rotation[pheno$pop == "stripe"], y=pheno$rotation[pheno$pop == "admix"], alternative="two.sided")

# dorsal body color: band vs stripe

t.test(x=pheno$body.a[pheno$pop == "band"], y=pheno$body.a[pheno$pop == "stripe"], alternative="two.sided")
wilcox.test(x=pheno$body.a[pheno$pop == "band"], y=pheno$body.a[pheno$pop == "stripe"], alternative="two.sided")

# dorsal body color: band vs admix
t.test(x=pheno$body.a[pheno$pop == "band"], y=pheno$body.a[pheno$pop == "admix"], alternative="two.sided")
wilcox.test(x=pheno$body.a[pheno$pop == "band"], y=pheno$body.a[pheno$pop == "admix"], alternative="two.sided")

# dorsal body color: stripe vs admix
t.test(x=pheno$body.a[pheno$pop == "stripe"], y=pheno$body.a[pheno$pop == "admix"], alternative="two.sided")
wilcox.test(x=pheno$body.a[pheno$pop == "stripe"], y=pheno$body.a[pheno$pop == "admix"], alternative="two.sided")

# leg color: band vs stripe
t.test(x=pheno$leg.a[pheno$pop == "band"], y=pheno$leg.a[pheno$pop == "stripe"], alternative="two.sided")
wilcox.test(x=pheno$leg.a[pheno$pop == "band"], y=pheno$leg.a[pheno$pop == "stripe"], alternative="two.sided")

# leg color: band vs admix
t.test(x=pheno$leg.a[pheno$pop == "band"], y=pheno$leg.a[pheno$pop == "admix"], alternative="two.sided")
wilcox.test(x=pheno$leg.a[pheno$pop == "band"], y=pheno$leg.a[pheno$pop == "admix"], alternative="two.sided")

# leg color: stripe vs admix
t.test(x=pheno$leg.a[pheno$pop == "stripe"], y=pheno$leg.a[pheno$pop == "admix"], alternative="two.sided")
wilcox.test(x=pheno$leg.a[pheno$pop == "stripe"], y=pheno$leg.a[pheno$pop == "admix"], alternative="two.sided")

### ALLELE FREQUENCY PIE CHARTS

af <- read.table('/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_17092020/allele_frequencies/high_confidence_candidate_gene_frequencies.txt', head=TRUE)
af$locus = paste(af$gene,af$snp,sep="_")
afu = af[sapply(unique(af$locus),function(x,df){which(df$locus == x)[1]},df=af),]
afu = afu[c(1,5,2,3,4,6,7,8),]

# stripe pie charts (best use ggplot for multipanel)
#par(mfrow=c(2,4), mar=c(0,0,0,0), oma=c(0,0,0,0))
#pie(c(afu$stripe_alt_freq[1], 1-afu$stripe_alt_freq[1]), labels=NA, col=c("darkturquoise","black")) # mc1r 1813
#pie(c(afu$stripe_alt_freq[2], 1-afu$stripe_alt_freq[2]), labels=NA, col=c("tomato1","black")) # asip 219
#pie(c(afu$stripe_alt_freq[3], 1-afu$stripe_alt_freq[3]), labels=NA, col=c("mediumpurple3","black")) # asip 1582
#pie(c(afu$stripe_alt_freq[4], 1-afu$stripe_alt_freq[4]), labels=NA, col=c("maroon2","black")) # arg2 2193
#pie(c(afu$stripe_alt_freq[5], 1-afu$stripe_alt_freq[5]), labels=NA, col=c("olivedrab","black")) # bsn 1746
#pie(c(afu$stripe_alt_freq[6], 1-afu$stripe_alt_freq[6]), labels=NA, col=c("khaki1","black")) # bsn 6952
#pie(c(afu$stripe_alt_freq[7], 1-afu$stripe_alt_freq[7]), labels=NA, col=c("plum1","black")) # retsat 6779
#pie(c(afu$stripe_alt_freq[8], 1-afu$stripe_alt_freq[8]), labels=NA, col=c("dodgerblue2","black")) # krt8.2 2306

piemar=c(0,0,0,0)

stripe.pie1 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$stripe_alt_freq[1], 1-afu$stripe_alt_freq[1])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("darkturquoise","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # mc1r 1813

stripe.pie2 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$stripe_alt_freq[2], 1-afu$stripe_alt_freq[2])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("tomato1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none",plot.margin=unit(piemar, "lines")) # asip 219

stripe.pie3 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$stripe_alt_freq[3], 1-afu$stripe_alt_freq[3])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("mediumpurple3","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # asip 1582

stripe.pie4 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$stripe_alt_freq[4], 1-afu$stripe_alt_freq[4])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("maroon2","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # arg2 2193

stripe.pie5 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$stripe_alt_freq[5], 1-afu$stripe_alt_freq[5])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("olivedrab","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # bsn 1746

stripe.pie6 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$stripe_alt_freq[6], 1-afu$stripe_alt_freq[6])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("khaki1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # bsn 6952

stripe.pie7 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$stripe_alt_freq[7], 1-afu$stripe_alt_freq[7])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("plum1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # retsat 6779

stripe.pie8 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$stripe_alt_freq[8], 1-afu$stripe_alt_freq[8])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("dodgerblue2","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # krt8.2 2306

png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_17092020/paper_figures/stripe_af_pie.png',width=1200,height=600, bg="transparent")
plot_grid(plotlist=list(stripe.pie1, stripe.pie2, stripe.pie3, stripe.pie4, stripe.pie5, stripe.pie6, stripe.pie7, stripe.pie8),nrow=2,ncol=4)
dev.off()

# band pie charts

piemar=c(0,0,0,0)

band.pie1 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$band_alt_freq[1], 1-afu$band_alt_freq[1])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("darkturquoise","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # mc1r 1813

band.pie2 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$band_alt_freq[2], 1-afu$band_alt_freq[2])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("tomato1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none",plot.margin=unit(piemar, "lines")) # asip 219

band.pie3 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$band_alt_freq[3], 1-afu$band_alt_freq[3])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("mediumpurple3","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # asip 1582

band.pie4 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$band_alt_freq[4], 1-afu$band_alt_freq[4])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("maroon2","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # arg2 2193

band.pie5 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$band_alt_freq[5], 1-afu$band_alt_freq[5])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("olivedrab","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # bsn 1746

band.pie6 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$band_alt_freq[6], 1-afu$band_alt_freq[6])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("khaki1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # bsn 6952

band.pie7 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$band_alt_freq[7], 1-afu$band_alt_freq[7])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("plum1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # retsat 6779

band.pie8 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$band_alt_freq[8], 1-afu$band_alt_freq[8])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("dodgerblue2","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # krt8.2 2306

png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_17092020/paper_figures/band_af_pie.png',width=1200,height=600, bg="transparent")
plot_grid(plotlist=list(band.pie1, band.pie2, band.pie3, band.pie4, band.pie5, band.pie6, band.pie7, band.pie8),nrow=2,ncol=4)
dev.off()

# admix pie charts

piemar=c(0,0,0,0)

admix.pie1 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$admix_alt_freq[1], 1-afu$admix_alt_freq[1])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("darkturquoise","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # mc1r 1813

admix.pie2 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$admix_alt_freq[2], 1-afu$admix_alt_freq[2])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("tomato1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none",plot.margin=unit(piemar, "lines")) # asip 219

admix.pie3 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$admix_alt_freq[3], 1-afu$admix_alt_freq[3])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("mediumpurple3","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # asip 1582

admix.pie4 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$admix_alt_freq[4], 1-afu$admix_alt_freq[4])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("maroon2","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # arg2 2193

admix.pie5 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$admix_alt_freq[5], 1-afu$admix_alt_freq[5])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("olivedrab","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # bsn 1746

admix.pie6 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$admix_alt_freq[6], 1-afu$admix_alt_freq[6])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("khaki1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # bsn 6952

admix.pie7 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$admix_alt_freq[7], 1-afu$admix_alt_freq[7])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("plum1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # retsat 6779

admix.pie8 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(afu$admix_alt_freq[8], 1-afu$admix_alt_freq[8])), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("dodgerblue2","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # krt8.2 2306

png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_17092020/paper_figures/admix_af_pie.png',width=1200,height=600, bg="transparent")
plot_grid(plotlist=list(admix.pie1, admix.pie2, admix.pie3, admix.pie4, admix.pie5, admix.pie6, admix.pie7, admix.pie8),nrow=2,ncol=4)
dev.off()

# make pie backgrounds

piemar=c(0,0,0,0)

pieback1 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(1, 0)), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("darkturquoise","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # mc1r 1813

pieback2 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(1, 0)), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("tomato1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none",plot.margin=unit(piemar, "lines")) # asip 219

pieback3 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(1, 0)), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("mediumpurple3","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # asip 1582

pieback4 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(1, 0)), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("maroon2","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # arg2 2193

pieback5 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(1, 0)), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("olivedrab","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # bsn 1746

pieback6 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(1, 0)), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("khaki1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # bsn 6952

pieback7 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(1, 0)), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("plum1","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # retsat 6779

pieback8 <- ggplot(data.frame(allele=factor(c("alt","ref")), freq=c(1, 0)), aes(x="", y=freq, fill=allele)) + geom_bar(stat="identity",width=1) + scale_fill_manual(values=c("dodgerblue2","black")) + coord_polar("y", start=0) + theme_transparent() + theme(legend.position="none", plot.margin=unit(piemar, "lines")) # krt8.2 2306

png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_17092020/paper_figures/af_pie_backgrounds.png',width=1200,height=600, bg="transparent")
plot_grid(plotlist=list(pieback1, pieback2, pieback3, pieback4, pieback5, pieback6, pieback7, pieback8),nrow=2,ncol=4, scale=1.1)
dev.off()

#
#
#

### PHENOTYPE DISTRIBUTIONS

## Correlations

# see 'pairs', 'cpairs' and 'chart.Correlation' of 'PerformanceAnalytics' package

# read phenotype measures for wild individuals from the admixture zone
pheno <- read.table('/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/admix/admix_phenotypes_all_bamorder.txt',head=TRUE)

# read in phenotype measures from pedigree individuals
ped.body <- read.table('/home/tyler/Dropbox/research/imitator/sanger/pedigree/results/ped_dorsum.txt',head=TRUE,sep="\t")
ped.leg <- read.table('/home/tyler/Dropbox/research/imitator/sanger/pedigree/results/ped_hindlimb.txt',head=TRUE,sep="\t")
ped.pheno <- cbind(ped.body[,c(1,2,3,4,14,15,16,17)],ped.leg[match(ped.body$Individual_ID,ped.leg$Individual_ID),c(15,16,17)])
colnames(ped.pheno) <- c("family","sample","parent1","parent2","rotation_score","body.B2","body.S1Y","body.S1R","leg.B2","leg.S1G","leg.S1R")
#write.table(ped.pheno, file='/home/tyler/Dropbox/research/imitator/sanger/pedigree/results/pedigree_phenotypes.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# define plotting functions (see https://r-coder.com/correlation-plot-r/)
#panel.hist <- function(x, ...) {
#   usr <- par("usr")
#   on.exit(par(usr))
#   par(usr = c(usr[1:2], 0, 1.5))
#   his <- hist(x, plot = FALSE)
#   breaks <- his$breaks
#   nB <- length(breaks)
#   y <- his$counts
#   y <- y/max(y)
#   rect(breaks[-nB], 0, breaks[-1], y, col = "slategray3", ...)
   # lines(density(x), col = 2, lwd = 2)
#}

# plot

# wild individuals
pdf(file="/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/admix/phenotype_plots/wild_admix_phenotype_correlations.pdf")
pairs(pheno[,c(2,3,4,5,7,8)], labels=c("Pattern", paste("Percent\nblack",sep=""), "Body L*", "Body a*", "Leg L*", "Leg a*"), upper.panel = NULL, cex.axis=1.2, col=rgb(123,6,44,maxColorValue=255))
dev.off()

# pedigree individuals
pdf(file="/home/tyler/Dropbox/research/imitator/sanger/pedigree/results/phenotype_plots/pedigree_phenotype_correlations.pdf")
pairs(ped.pheno[,c(5,6,7,8,9,10,11)], labels=c("Pattern",paste("Body\nB2",sep=""), paste("Body\nS1Y",sep=""), paste("Body\nS1R",sep="\n"), paste("Leg\nB2",sep=""), paste("Leg\nS1G",sep=""), paste("Leg\nS1R",sep="")), upper.panel=NULL, cex.axis=1.2, col=rgb(20,50,112,maxColorValue=255))
dev.off()

## Distributions

# wild individuals

pdf(file="/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/admix/phenotype_plots/wild_admix_phenotype_distributions.pdf")
par(mfrow=c(3,3),mar=c(4.5,4.5,3.1,1.1))
hist(pheno$pattern_rotation, breaks=10, col=rgb(123,6,44,maxColorValue=255),main=NULL, xlab="Pattern",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(pheno$black_proportion, breaks=5, col=rgb(123,6,44,maxColorValue=255),main=NULL, xlab="Percent black",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(pheno$body_L, breaks=10, col=rgb(123,6,44,maxColorValue=255),main=NULL, xlab="Body L*",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(pheno$body_a, breaks=10, col=rgb(123,6,44,maxColorValue=255),main=NULL, xlab="Body a*",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(pheno$leg_L, breaks=10, col=rgb(123,6,44,maxColorValue=255),main=NULL, xlab="Leg L*",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(pheno$leg_a, breaks=10, col=rgb(123,6,44,maxColorValue=255),main=NULL, xlab="Leg a*",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
dev.off()

# pedigree individuals

pdf(file="/home/tyler/Dropbox/research/imitator/sanger/pedigree/results/phenotype_plots/pedigree_phenotype_distributions.pdf")
par(mfrow=c(3,3),mar=c(4.5,4.5,3.1,1.1))
hist(ped.pheno$rotation_score, breaks=8, col=rgb(20,50,112,maxColorValue=255),main=NULL, xlab="Pattern",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
plot.new()
plot.new()
hist(ped.pheno$body.B2, breaks=10, col=rgb(20,50,112,maxColorValue=255),main=NULL, xlab="Body B2",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(ped.pheno$body.S1Y, breaks=5, col=rgb(20,50,112,maxColorValue=255),main=NULL, xlab="Body S1Y",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(ped.pheno$body.S1R, breaks=5, col=rgb(20,50,112,maxColorValue=255),main=NULL, xlab="Body S1R",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(ped.pheno$leg.B2, breaks=10, col=rgb(20,50,112,maxColorValue=255),main=NULL, xlab="Leg B2",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(ped.pheno$leg.S1G, breaks=10, col=rgb(20,50,112,maxColorValue=255),main=NULL, xlab="Leg S1G",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
hist(ped.pheno$leg.S1R, breaks=10, col=rgb(20,50,112,maxColorValue=255),main=NULL, xlab="Leg S1R",ylab="Frequency",cex.lab=1.6, cex.axis=1.4)
dev.off()

#
#
#

### PHENOTYPE VS GENOTYPE PLOTS FOR CANDIDATE LOCI

dat <- read.table('/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/gwas/candidate_variance/candidate_geno_pheno_data.txt',head=TRUE,sep="\t",quote="")

#
#
#

### PEDIGREE PHENOTYPE VS GENOTYPE BOXPLOTS

library(scales)

### === DORSUM === ###

## PARENTALS, F1s, F2s

dorsum <- read.table('/home/tyler/Dropbox/research/imitator/sanger/pedigree/results/ped_dorsum.txt',head=TRUE,sep="\t")
pattern_range = c(min(dorsum$rotation_score, na.rm=TRUE),max(dorsum$rotation_score, na.rm=TRUE))
dorsum_b2_range=c(min(dorsum$B2, na.rm=TRUE), max(dorsum$B2, na.rm=TRUE))
dorsum_s1y_range=c(min(dorsum$S1Y, na.rm=TRUE), max(dorsum$S1Y, na.rm=TRUE))

# define colors based on Parental, F1, or F2
f2col = "hotpink"
f1col = "steelblue2"
pcol = "gold2"
dorsum$gencol <- as.character(dorsum$Individual_ID)
dorsum$gencol[grep("-",dorsum$gencol)] <- f2col # F2
dorsum$gencol[grep("x",dorsum$gencol)] <- f1col # F1
dorsum$gencol[which(dorsum$gencol != f2col & dorsum$gencol != f1col)] <- pcol # Parental

## MC1R
dorsum$mc1r_geno = paste0(dorsum$mc1r1, dorsum$mc1r2)
dorsum$mc1r_code <- replace(dorsum$mc1r_geno, which(dorsum$mc1r_geno == "CC"),1)
dorsum$mc1r_code <- replace(dorsum$mc1r_code, which(dorsum$mc1r_code == "CT"),2)
dorsum$mc1r_code <- replace(dorsum$mc1r_code, which(dorsum$mc1r_code == "TT"),3)
dorsum$mc1r_code <- factor(dorsum$mc1r_code, levels=c("1","2","3"))

# mc1r:pattern
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/mc1r_pattern.png', width=800, height=600)
par(mfrow=c(1,2)) # this is just to make the plots with one and two panels match

par(mar=c(5,5,6,5))
boxplot(dorsum$rotation_score ~ dorsum$mc1r_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="Pattern rotation score", main="mc1r", 
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=pattern_range)
for (i in levels(dorsum$mc1r_code)) {
	idx = which(as.character(dorsum$mc1r_code) == i)
	xjitter = jitter(rep(which(levels(dorsum$mc1r_code) == i),length(idx)), amount=0.2)
	points(rotation_score[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
	#points(rotation_score[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
	points(rotation_score[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}
dev.off()

# mc1r:B2,S1Y
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/mc1r_dorsum_color.png', width=800,height=600)
par(mfrow=c(1,2))

par(mar=c(5,5,6,5)) # mc1r:B2
boxplot(dorsum$B2 ~ dorsum$mc1r_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="B2", main="mc1r", col="gray80", 
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=dorsum_b2_range)
for (i in levels(dorsum$mc1r_code)) {
        idx = which(as.character(dorsum$mc1r_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$mc1r_code) == i),length(idx)), amount=0.2)
        points(B2[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

par(mar=c(5,5,6,5)) # mc1r:S1Y
boxplot(dorsum$S1Y ~ dorsum$mc1r_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="S1Y", main="mc1r", col=alpha("yellow2",0.2), 
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=dorsum_s1y_range)
for (i in levels(dorsum$mc1r_code)) {
        idx = which(as.character(dorsum$mc1r_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$mc1r_code) == i),length(idx)), amount=0.2)
        points(S1Y[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(S1Y[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(S1Y[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

dev.off()

## ASIP
dorsum$asip_geno = paste0(dorsum$asip1,dorsum$asip2)
dorsum$asip_code <- replace(dorsum$asip_geno, which(dorsum$asip_geno == "CC"),1)
dorsum$asip_code <- replace(dorsum$asip_code, which(dorsum$asip_code == "CT"),2)
dorsum$asip_code <- replace(dorsum$asip_code, which(dorsum$asip_code == "TT"),3)
dorsum$asip_code <- factor(dorsum$asip_code, levels=c("1","2","3"))

# asip:pattern
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/asip_pattern.png', width=800,height=600)
par(mfrow=c(1,2))

par(mar=c(5,5,6,5))
boxplot(dorsum$rotation_score ~ dorsum$asip_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="Pattern rotation score", main="asip",
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=pattern_range)
for (i in levels(dorsum$asip_code)) {
        idx = which(as.character(dorsum$asip_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$asip_code) == i),length(idx)), amount=0.2)
        points(rotation_score[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(rotation_score[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(rotation_score[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

dev.off()

# asip:B2
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/asip_dorsum_color.png', width=800,height=600)
par(mfrow=c(1,2))

par(mar=c(5,5,6,5))
boxplot(dorsum$B2 ~ dorsum$asip_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="B2", main="asip", col="gray80",
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=dorsum_b2_range)
for (i in levels(dorsum$asip_code)) {
        idx = which(as.character(dorsum$asip_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$asip_code) == i),length(idx)), amount=0.2)
        points(B2[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

dev.off()

## BSN
dorsum$bsn_geno = paste0(dorsum$bsn2.1, dorsum$bsn2.2)
dorsum$bsn_code <- replace(dorsum$bsn_geno, which(dorsum$bsn_geno == "AA"),1)
dorsum$bsn_code <- replace(dorsum$bsn_code, which(dorsum$bsn_code == "AC"),2)
dorsum$bsn_code <- replace(dorsum$bsn_code, which(dorsum$bsn_code == "CC"),3)
dorsum$bsn_code <- factor(dorsum$bsn_code, levels=c("1","2","3"))

# bsn:pattern
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/bsn_pattern.png', width=800,height=600)
par(mfrow=c(1,2))

par(mar=c(5,5,6,5))
boxplot(dorsum$rotation_score ~ dorsum$bsn_code, outline=FALSE, names=c("AA", "AC", "CC"), xlab="", ylab="Pattern rotation score", main="bsn",
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=pattern_range)
for (i in levels(dorsum$bsn_code)) {
        idx = which(as.character(dorsum$bsn_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$bsn_code) == i),length(idx)), amount=0.2)
        points(rotation_score[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(rotation_score[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(rotation_score[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

dev.off()

# bsn:B2,S1Y
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/bsn_dorsum_color.png', width=800,height=600)
par(mfrow=c(1,2))

par(mar=c(5,5,6,5)) # bsn:B2
boxplot(dorsum$B2 ~ dorsum$bsn_code, outline=FALSE, names=c("AA", "AC", "CC"), xlab="", ylab="B2", main="bsn", col="gray80",
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=dorsum_b2_range)
for (i in levels(dorsum$bsn_code)) {
        idx = which(as.character(dorsum$bsn_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$bsn_code) == i),length(idx)), amount=0.2)
        points(B2[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

par(mar=c(5,5,6,5)) # bsn:S1Y
boxplot(dorsum$S1Y ~ dorsum$bsn_code, outline=FALSE, names=c("AA", "AC", "CC"), xlab="", ylab="S1Y", main="bsn", col=alpha("yellow2",0.2),
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=dorsum_s1y_range)
for (i in levels(dorsum$bsn_code)) {
        idx = which(as.character(dorsum$bsn_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$bsn_code) == i),length(idx)), amount=0.2)
        points(S1Y[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(S1Y[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(S1Y[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

dev.off()

## RETSAT
dorsum$retsat_geno = paste0(dorsum$retsat1, dorsum$retsat2)
dorsum$retsat_code <- replace(dorsum$retsat_geno, which(dorsum$retsat_geno == "CC"),1)
dorsum$retsat_code <- replace(dorsum$retsat_geno, which(dorsum$retsat_code == "CT"),2)
dorsum$retsat_code <- replace(dorsum$retsat_code, which(dorsum$retsat_code == "TT"),3)
dorsum$retsat_code <- factor(dorsum$retsat_code, levels=c("1","2","3"))

### !! Do this to remove family 3 because of genotype errors !! ###
dorsum <- dorsum[-which(dorsum$Fam == 3),]

# retsat:B2,S1Y
#png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/retsat_dorsum_color.png', width=800,height=600)
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/retsat_dorsum_color_nofam3.png', width=800,height=600)
par(mfrow=c(1,2))

par(mar=c(5,5,6,5)) # retsat:B2
b2_lim1=c(min(dorsum$B2[which(!is.na(dorsum$retsat_code))], na.rm=TRUE), max(dorsum$B2[which(!is.na(dorsum$retsat_code))], na.rm=TRUE))
boxplot(dorsum$B2 ~ dorsum$retsat_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="B2", main="retsat", col="gray80",
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=b2_lim1)
for (i in levels(dorsum$retsat_code)) {
        idx = which(as.character(dorsum$retsat_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$retsat_code) == i),length(idx)), amount=0.2)
        points(B2[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

par(mar=c(5,5,6,5)) # retsat:S1Y
boxplot(dorsum$S1Y ~ dorsum$retsat_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="S1Y", main="retsat", col=alpha("yellow2",0.2),
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=dorsum_s1y_range)
for (i in levels(dorsum$retsat_code)) {
        idx = which(as.character(dorsum$retsat_code) == i)
        xjitter = jitter(rep(which(levels(dorsum$retsat_code) == i),length(idx)), amount=0.2)
        points(S1Y[idx] ~ xjitter, pch=1, cex=1.2, data=dorsum)
        #points(S1Y[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=dorsum) # all maroon
        points(S1Y[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=dorsum) # color by generation
}

dev.off()

### === HINDLIMB === ###

leg <- read.table('/home/tyler/Dropbox/research/imitator/sanger/pedigree/results/ped_hindlimb.txt',head=TRUE,sep="\t")
leg_b2_range=c(min(leg$B2, na.rm=TRUE), max(leg$B2, na.rm=TRUE))
leg_s1g_range=c(min(leg$S1G, na.rm=TRUE), max(leg$S1G, na.rm=TRUE))

# define colors based on Parental, F1, F2s

leg$gencol <- as.character(leg$Individual_ID)
leg$gencol[grep("-",leg$gencol)] <- f2col # F2
leg$gencol[grep("x",leg$gencol)] <- f1col # F1
leg$gencol[which(leg$gencol != f2col & leg$gencol != f1col)] <- pcol # Parental

## MC1R
leg$mc1r_geno = paste0(leg$mc1r1, leg$mc1r2)
leg$mc1r_code <- replace(leg$mc1r_geno, which(leg$mc1r_geno == "CC"),1)
leg$mc1r_code <- replace(leg$mc1r_code, which(leg$mc1r_code == "CT"),2)
leg$mc1r_code <- replace(leg$mc1r_code, which(leg$mc1r_code == "TT"),3)
leg$mc1r_code <- factor(leg$mc1r_code, levels=c("1","2","3"))

# mc1r:B2
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/mc1r_leg_color.png', width=800,height=600)
par(mfrow=c(1,2))

par(mar=c(5,5,6,5))
boxplot(leg$B2 ~ leg$mc1r_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="B2", main="mc1r", col="gray80",
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=leg_b2_range)
for (i in levels(leg$mc1r_code)) {
        idx = which(as.character(leg$mc1r_code) == i)
        xjitter = jitter(rep(which(levels(leg$mc1r_code) == i),length(idx)), amount=0.2)
        points(B2[idx] ~ xjitter, pch=1, cex=1.2, data=leg)
        #points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=leg) # all maroon
        points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=leg) # color by generation
}

dev.off()

## ASIP
leg$asip_geno = paste0(leg$asip1,leg$asip2)
leg$asip_code <- replace(leg$asip_geno, which(leg$asip_geno == "CC"),1)
leg$asip_code <- replace(leg$asip_code, which(leg$asip_code == "CT"),2)
leg$asip_code <- replace(leg$asip_code, which(leg$asip_code == "TT"),3)
leg$asip_code <- factor(leg$asip_code, levels=c("1","2","3"))

# asip:B2,S1G
png(file='/home/tyler/Dropbox/research/imitator/demultiplex2/spades/round2/sites_290421/paper_figures/pedigree_boxplots/asip_leg_color.png', width=800,height=600)
par(mfrow=c(1,2))

par(mar=c(5,5,6,5)) # asip:B2
boxplot(leg$B2 ~ leg$asip_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="B2", main="asip", col="gray80",
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=leg_b2_range)
for (i in levels(leg$asip_code)) {
        idx = which(as.character(leg$asip_code) == i)
        xjitter = jitter(rep(which(levels(leg$asip_code) == i),length(idx)), amount=0.2)
        points(B2[idx] ~ xjitter, pch=1, cex=1.2, data=leg)
        #points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=leg) # all maroon
        points(B2[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=leg) # color by generation
}

par(mar=c(5,5,6,5)) # asip:S1G
boxplot(leg$S1G ~ leg$asip_code, outline=FALSE, names=c("CC", "CT", "TT"), xlab="", ylab="S1G", main="asip", col=alpha("palegreen3",0.4),
cex.axis=1.5, cex.lab=1.5, lwd=2.5, range=1.5, ylim=leg_s1g_range)
for (i in levels(leg$asip_code)) {
        idx = which(as.character(leg$asip_code) == i)
        xjitter = jitter(rep(which(levels(leg$asip_code) == i),length(idx)), amount=0.2)
        points(S1G[idx] ~ xjitter, pch=1, cex=1.2, data=leg)
        #points(S1G[idx] ~ xjitter, pch=16, cex=1.2, col=alpha("maroon",0.6), data=leg) # all maroon
        points(S1G[idx] ~ xjitter, pch=16, cex=1.2, col=alpha(gencol[idx],0.6), data=leg) # color by generation
}

dev.off()
