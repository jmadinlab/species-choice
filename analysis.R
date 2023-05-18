# Species choice

# LIBRARIES
library("vegan")
library("spatstat")
library("psych")
library("raster")
library("plotrix")

# FUNCTIONS
source("R/functions.R")

# DATA PREPARATION
source("R/data_prep.R")

# TRAIT SPACE
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat)

fit <- envfit(dat[c("PC1", "PC2")], dat[,cats]) # use envfit for drawing arrows
fit2 <- fit$vectors$arrows*-3 # drawing line segments opposite arrows

# FIGURE 1, PCAs

png("figs/fig_1.png", width = 6.5, height = 6, units = 'in', res = 300)


m <- rbind(c(1, 1, 1, 2, 2, 2), 
           c(1, 1, 1, 2, 2, 2),
           c(1, 1, 1, 2, 2, 2),
           c(3, 3, 4, 4, 5, 5),
           c(3, 3, 4, 4, 5, 5))

layout(m)
par(mar=c(3, 2, 2, 0), oma=c(3, 3, 0, 1))

# Loadings
plot(PC2 ~ PC1, dat, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("A. Trait loadings", line=-1, adj=0, cex=0.9)
mtext("PC1 (47.2%)", 1, 2, cex=0.8)
mtext("PC2 (18.6%)", 2, 2, cex=0.8)

r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
xx <- seq(xmn, xmx, length=6)[1:5]
yy <- seq(ymx, ymn, length=6)[1:5]
xxd <- diff(xx)[1]/2
yyd <- diff(yy)[1]/2
xx <- xx + xxd
yy <- yy + yyd
xxv <- rep(xx, length(yy))
yyv <- rep(yy, each=length(xx))
cells <- unique(cellFromXY(r, cbind(dat$PC1, dat$PC2)))
rect(xxv[cells] - xxd, yyv[cells] - yyd, xxv[cells] + xxd, yyv[cells] + yyd, border="lightblue")

plot(fit, cex=0.90, col=1, labels=list(vectors = c("GR","CW","MCS","SD","CH","SAV","R")))

## Schematic 

plot(PC2 ~ PC1, dat, pch=20, cex=1, col="lightgrey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
draw.ellipse(0.2, 0.2, 3.9, 2.7, 0, border=NA, col=rgb(0, 0, 0, 0.1))
draw.ellipse(2.3, 0.0, 1.8, 1.2, -85, border=NA, col=adjustcolor(cols[1], alpha.f = 0.3))
draw.ellipse(-1.2, -0.9, 2.2, 0.9, -25, border=NA, col=adjustcolor(cols[1], alpha.f = 0.3))
draw.ellipse(-0.8, 1.5, 2.1, 0.9, 15, border=NA, col=adjustcolor(cols[1], alpha.f = 0.3))
points(PC2 ~ PC1, dat, pch=20, cex=1, col="lightgrey")
text(-0.7, 0.3, "Coastal\nprotection", adj=0, cex=0.8)
text(1.4, 0.3, "Habitat provision\n(small fishes)", adj=0, cex=0.8)
text(-2.1, -1.0, "Habitat provision\n(large fishes)", adj=0, cex=0.8)
text(-2.5, 1.4, "Reef stability", adj=0, cex=0.8)
mtext("B. Ecosystem services", line=-1, adj=0, cex=0.9)

# Goreau

plot(PC2 ~ PC1, dat, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
oo <- ordiellipse(pca, dat$goreau, conf=0.95, col=cols[3])
for (i in 1:length(oo)) {
  text(oo[[i]][[2]][1], oo[[i]][[2]][2], names(oo[i]), cex=0.8)
}
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("C. Reef building", line=-1, adj=0, cex=0.9)

# Life history

plot(PC2 ~ PC1, dat, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
oo <- ordiellipse(pca, dat$life_history, conf=0.95, col=cols[1])
for (i in 1:length(oo)) {
  text(oo[[i]][[2]][1], oo[[i]][[2]][2], names(oo[i]), cex=0.8)
}
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("D. Life history strategies", line=-1, adj=0, cex=0.9)

# Phylogeny
plot(PC2 ~ PC1, dat, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
oo <- ordiellipse(pca, dat$family_molecules, conf=0.75, col=cols[2])
for (i in 1:length(oo)) {
  text(oo[[i]][[2]][1], oo[[i]][[2]][2], names(oo[i]), cex=0.75)
}

axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("E. Phylogenetic diversity", line=-1, adj=0, cex=0.9)

dev.off()

### FIGURE 2, OCCUPANCY

store <- read.csv("output/occupancy.csv", as.is=TRUE)

png("figs/fig_2.png", width = 5.5, height = 5.5, units = 'in', res = 300)

m <- rbind(c(1, 1, 1, 1), 
           c(1, 1, 1, 1),
           c(2, 2, 3, 3),
           c(2, 2, 3, 3))
layout(m)
par(mar=c(3, 2, 2, 0), oma=c(0, 0, 0, 1))

n5_max <- max(store$celln5_f)
m5_max <- max(store$cellm5_f)

par(mar=c(3, 5, 3, 1))

plot((store$n), (store$celln5_rand/n5_max), type="l", xlab="", ylab="Proportion cell occupancy", axes=FALSE, lty=1, xlim=c(0, 60), ylim=c(0, 1))
polygon(c((store$n), rev((store$n))), c((store$celln5_rand_lo/n5_max), rev((store$celln5_rand_hi/n5_max))), border=NA, col=rgb(0,0,0,0.1))

mtext("A. Hedging and marginal returns", 3, 0, adj=0, cex=0.9)

axis(2, las=2)
axis(1, las=1)

lines((store$n), (store$celln5_e/n5_max), col=cols[1])
lines((store$n), (store$celln5_f/n5_max), col=cols[2])
lines((store$n), (store$celln5_fe/n5_max), col=cols[3])

abline(v=11, lty=3)
abline(v=28, lty=3)

text(5, 1, "High", cex=0.8)
text(19, 1, "Medium", cex=0.8)
text(34, 1, "Low", cex=0.8)

points(sum(dat$lhi), (store$celln5_lhi[!is.na(store$celln5_lhi)]/n5_max), col="black", pch=8)
points(sum(dat$cnp), (store$celln5_cnp[!is.na(store$celln5_cnp)]/n5_max), col="black", pch=9)


legend("bottomright", legend=c("Random selection", "Persistence", "Trait diversity", "Persistence & traits", "Lord Howe Island", "Coral Nurture Program"), fill=c("lightgrey", NA, NA, NA, NA, NA), border=NA, lty=c(1, 1, 1, 1, NA, NA), pch=c(NA, NA, NA, NA, 8, 9), col=c("black", cols[1], cols[2], cols[3], "black", "black"), bty="n", cex=0.7, xpd=TRUE)
mtext("Number of species", 1, 2, cex=0.8)

s <- 15
lhi_choice <- choice(dat, s, vars=c("abund", "range", "bleach"), method="nearest", trait=TRUE)

plot(PC2 ~ PC1, dat, pch=20, cex=0.5, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE)
lhi <- dat[dat$lhi,]
oo <- ordiellipse(pca, dat$goreau, conf=0.95, col=cols[3])
points(PC2 ~ PC1, lhi, col="black", pch=20, cex=1)
points(PC2 ~ PC1, lhi_choice, pch=3, cex=1, col=cols[1])
mtext("B. Lord Howe Island", line=-1, adj=0, cex=0.9)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("PC1", 1, 2, cex=0.8)
mtext("PC2", 2, 2, cex=0.8)
text(-3, -2.5, "n=15", cex=0.8)

s <- 39
cnp_choice <- choice(dat, s, vars=c("abund", "range", "bleach"), method="nearest", trait=TRUE)

plot(PC2 ~ PC1, dat, pch=20, cex=0.5, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE)
cnp <- dat[dat$cnp,]
oo <- ordiellipse(pca, dat$goreau, conf=0.95, col=cols[3])
points(PC2 ~ PC1, cnp, col="black", pch=20, cex=1)
points(PC2 ~ PC1, cnp_choice, pch=3, cex=1, col=cols[1])
mtext("C. Coral Nurture Program", line=-1, adj=0, cex=0.9)
gfas <- dat[dat$cnp & dat$species == "Galaxea fascicularis",]
arrows(-1, -2.1, gfas$PC1, gfas$PC2-0.1, length=0.1, col=cols[3], lwd=2)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
text(-3, -2.5, "n=39", cex=0.8)

dev.off()

### FIGURE 3, Hedging stages

s = 20

choice_1 <- choice(dat, s, vars=c("abund", "range", "bleach"), trait=FALSE)[c("species_n", "value")]
names(choice_1) <- c("Species", "Ecological persistence")
choice_2 <- choice(dat, s, vars=c("abund", "range", "bleach"), trait=TRUE)[c("species_n", "value")]
names(choice_2) <- c("Species", "Persistence & traits")
choice_3 <- choice(dat, s, vars=c("abund", "range", "bleach", "restore"), trait=TRUE)[c("species_n", "value")]
names(choice_3) <- c("Species", "Persistence, traits & restoration")

mat <- merge(choice_1, choice_2, all=TRUE)
mat <- merge(mat, choice_3, all=TRUE)
mat <- mat[order(mat["Ecological persistence"], decreasing=FALSE),]
mat <- mat[order(mat["Persistence & traits"], decreasing=FALSE),]
mat <- mat[order(mat["Persistence, traits & restoration"], decreasing=FALSE),]

choice_4 <- choice(dat, s, vars=c("abund", "range", "bleach"), trait=FALSE, oppo=TRUE)[c("species_n", "value")]
names(choice_4) <- c("Species", "Vulnerable to bleaching")
choice_5 <- choice(dat, s, vars=c("abund", "range", "bleach"), trait=TRUE, oppo=TRUE)[c("species_n", "value")]
names(choice_5) <- c("Species", "Vulnerable & traits")
choice_6 <- choice(dat, s, vars=c("abund", "range", "bleach", "restore"), trait=TRUE, oppo=TRUE)[c("species_n", "value")]
names(choice_6) <- c("Species", "Vulnerable, traits & restoration")

mat2 <- merge(choice_4, choice_5, all=TRUE)
mat2 <- merge(mat2, choice_6, all=TRUE)
mat2 <- mat2[order(mat2[,2], decreasing=FALSE),]
mat2 <- mat2[order(mat2[,3], decreasing=FALSE),]
mat2 <- mat2[order(mat2[,4], decreasing=FALSE),]

png("figs/fig_3.png", width = 5, height = 7, units = 'in', res = 300)

par(mfrow=c(1,2), mar=c(1, 6, 8.2, 1))
image(t(as.matrix(mat[,2:4])), axes=FALSE, col=hcl.colors(20, "Heat", rev = TRUE))
mtext("A", line=4, cex=1.5, at=-1)

text(rep(0, nrow(mat[,2:4])), seq(0, 1, (1/(nrow(mat[,2:4])-1))), (round(as.matrix(mat[,2:4][,1]), 2)), cex=0.6)
text(rep(0.5, nrow(mat[,2:4])), seq(0, 1, (1/(nrow(mat[,2:4])-1))), (round(as.matrix(mat[,2:4][,2]), 2)), cex=0.6)
text(rep(1, nrow(mat[,2:4])), seq(0, 1, (1/(nrow(mat[,2:4])-1))), (round(as.matrix(mat[,2:4][,3]), 2)), cex=0.6)

axis(2, at=seq(0, 1, (1/(nrow(mat)-1))), labels=mat$Species, las=2, cex.axis=0.7, font=3)
axis(3, at=c(0, 0.5, 1), labels=c("Ecological persistence", "Persistence & traits", "Persistence, traits\n& restoration"), cex.axis=0.8, las=2)

image(t(as.matrix(mat2[,2:4])), axes=FALSE, col=hcl.colors(20, "Heat", rev = TRUE))
mtext("B", line=4, cex=1.5, at=-1)

text(rep(0, nrow(mat2[,2:4])), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,2]), 2)), cex=0.6)
text(rep(0.5, nrow(mat2[,2:4])), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,3]), 2)), cex=0.6)
text(rep(1, nrow(mat2[,2:4])), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,4]), 2)), cex=0.6)

axis(2, at=seq(0, 1, (1/(nrow(mat2)-1))), labels=mat2$Species, las=2, cex.axis=0.7, font=3)
axis(3, at=c(0, 0.5, 1), labels=c("Vulnerable to bleaching", "Vulnerable & traits", "Vulnerable, traits\n& restoration"), cex.axis=0.8, las=2)

dev.off()

## FIGURE S1, method comparison

mat_vv <- read.csv("output/mat_vv.csv", as.is=TRUE)
mat_nn <- read.csv("output/mat_nn.csv", as.is=TRUE)
mat_hh <- read.csv("output/mat_hh.csv", as.is=TRUE)

n <- 20
i <- dim(dat)[1] + 1 - n

png("figs/fig_S1.png", width = 4.5, height = 4, units = 'in', res = 300)
par(mar=c(3, 3, 0, 0))

plot(PC2 ~ PC1, dat, pch=20, cex=1, col="lightgrey", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
points(PC2 ~ PC1, dat[mat_nn[,i],], col=cols[1], pch=1, cex=1)
points(PC2 ~ PC1, dat[mat_vv[,i],], col=cols[2], pch=1, cex=1.5)
points(PC2 ~ PC1, dat[mat_hh[,i],], col=cols[3], pch=3, cex=1)
mtext("n = 20", line=-1, adj=0)
mtext("PC1 (47.2%)", 1, 2, cex=0.8)
mtext("PC2 (18.6%)", 2, 2, cex=0.8)

legend("topright", pch=c(1, 1, 3), legend=c("Nearest neighbor", "Voronoi area", "Hypervolume contribution"), bty="n", col=c(cols[1], cols[2], cols[3]), cex=0.9)

dev.off()

## FIGURE S2, 10 cell occupancy

n10_max <- max(store$celln10_f)
m10_max <- max(store$cellm10_f)

png("figs/fig_S2.png", width = 4, height = 6, units = 'in', res = 300)

par(mfrow=c(2, 1), oma=c(3, 0, 0, 0), mar=c(2, 4, 2, 1))

plot((store$n), (store$celln10_rand/n10_max), type="l", xlab="", ylab="Proportion cell occupancy", axes=FALSE, lty=1, xlim=c(0, 100), ylim=c(0, 1))
polygon(c((store$n), rev((store$n))), c((store$celln10_rand_lo/n10_max), rev((store$celln10_rand_hi/n10_max))), border=NA, col=rgb(0,0,0,0.1))

mtext("(a)", 3, 0, adj=0)

axis(2, las=2)
axis(1, at=seq(0, 100, 10), labels=NA, las=2)

lines((store$n), (store$celln10_e/n10_max), col=cols[1])
lines((store$n), (store$celln10_f/n10_max), col=cols[2])
lines((store$n), (store$celln10_fe/n10_max), col=cols[3])
abline(v=29, lty=3)
abline(v=55, lty=3)

text(12, 1, "High", cex=0.8)
text(42, 1, "Medium", cex=0.8)
text(66, 1, "Low", cex=0.8)

points((15), (store$celln10_lhi[!is.na(store$celln10_lhi)]/n10_max), col="black", pch=8)
points((39), (store$celln10_cnp[!is.na(store$celln10_cnp)]/n10_max), col="black", pch=9)

plot((store$n), (store$cellm10_rand), type="l", xlab="Number of species", ylab="Mean cell redundancy", axes=FALSE, lty=1, xlim=c(0, 100), ylim=c(1, 3))
polygon(c((store$n), rev((store$n))), c((store$cellm10_rand_lo), rev((store$cellm10_rand_hi))), border=NA, col=rgb(0,0,0,0.1))

mtext("(b)", 3, 0, adj=0)

axis(2, las=2)
axis(1, at=seq(0, 100, 10), las=2)

lines((store$n), (store$cellm10_e), col=cols[1])
lines((store$n), (store$cellm10_f), col=cols[2])
lines((store$n), (store$cellm10_fe), col=cols[3])
abline(v=29, lty=3)
abline(v=55, lty=3)

points((15), (store$cellm10_lhi[!is.na(store$cellm10_lhi)]), col="black", pch=8)
points((39), (store$cellm10_cnp[!is.na(store$cellm10_cnp)]), col="black", pch=9)

legend("topright", legend=c("Random selection", "Persistence", "Trait diversity", "Persistence & traits", "Lord Howe Island", "Coral Nurture Program"), lty=c(1, 1, 1, 1, NA, NA), pch=c(NA, NA, NA, NA, 8, 9), col=c("black", cols[1], cols[2], cols[3], "black", "black"), bty="n", cex=0.5, xpd=TRUE)

mtext("Number of species", 1, 1, outer=TRUE)

dev.off()

# FIGURE S3, CORRELATIONS AMONG VARS

png("figs/fig_S3.png", width = 6, height = 6, units = 'in', res = 300)
pairs.panels(dat[c("abund", "range", "bleach", "PC1", "PC2")], hist.col="lightgrey")
dev.off()

