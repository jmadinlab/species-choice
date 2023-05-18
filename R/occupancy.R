# Occupancy

### NUMBERS

ns <- seq(100, 3, -1)
store <- data.frame() 

r5 <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
r10 <- raster(extent(xmn, xmx, ymn, ymx), nrow = 10, ncol = 10)

for (s in ns) {
  ch_f <- choice(dat, s=s, vars=c(), method="nearest", trait=TRUE)
  ch_e <- choice(dat, s=s, vars=c("abund", "range", "bleach"), method="nearest", trait=FALSE)
  ch_fe <- choice(dat, s=s, vars=c("abund", "range", "bleach"), method="nearest", trait=TRUE)
  ch_fer <- choice(dat, s=s, vars=c("abund", "range", "bleach", "restore"), method="nearest", trait=TRUE)
  
  cells <- cellFromXY(r5, cbind(ch_f$PC1, ch_f$PC2))
  celln5_f <- length(unique(cells))
  cellm5_f <- mean(table(cells))
  cells <- cellFromXY(r10, cbind(ch_f$PC1, ch_f$PC2))
  celln10_f <- length(unique(cells))
  cellm10_f <- mean(table(cells))
  
  cells <- cellFromXY(r5, cbind(ch_e$PC1, ch_e$PC2))
  celln5_e <- length(unique(cells))
  cellm5_e <- mean(table(cells))
  cells <- cellFromXY(r10, cbind(ch_e$PC1, ch_e$PC2))
  celln10_e <- length(unique(cells))
  cellm10_e <- mean(table(cells))
  
  cells <- cellFromXY(r5, cbind(ch_fe$PC1, ch_fe$PC2))
  celln5_fe <- length(unique(cells))
  cellm5_fe <- mean(table(cells))
  cells <- cellFromXY(r10, cbind(ch_fe$PC1, ch_fe$PC2))
  celln10_fe <- length(unique(cells))
  cellm10_fe <- mean(table(cells))
  
  cells <- cellFromXY(r5, cbind(ch_fer$PC1, ch_fer$PC2))
  celln5_fer <- length(unique(cells))
  cellm5_fer <- mean(table(cells))
  cells <- cellFromXY(r10, cbind(ch_fer$PC1, ch_fer$PC2))
  celln10_fer <- length(unique(cells))
  cellm10_fer <- mean(table(cells))
  
  celln5_lhi<- cellm5_lhi<- celln10_lhi<- cellm10_lhi<- NA
  if (s == 15) {
    ch_lhi <- dat[dat$lhi,]
    cells <- cellFromXY(r5, cbind(ch_lhi$PC1, ch_lhi$PC2))
    celln5_lhi <- length(unique(cells))
    cellm5_lhi <- mean(table(cells))
    cells <- cellFromXY(r10, cbind(ch_lhi$PC1, ch_lhi$PC2))
    celln10_lhi <- length(unique(cells))
    cellm10_lhi <- mean(table(cells))
  }
  
  celln5_cnp <- cellm5_cnp <- celln10_cnp <- cellm10_cnp <- NA
  if (s == 39) {
    ch_cnp <- dat[dat$cnp,]
    cells <- cellFromXY(r5, cbind(ch_cnp$PC1, ch_cnp$PC2))
    celln5_cnp <- length(unique(cells))
    cellm5_cnp <- mean(table(cells))
    cells <- cellFromXY(r10, cbind(ch_cnp$PC1, ch_cnp$PC2))
    celln10_cnp <- length(unique(cells))
    cellm10_cnp <- mean(table(cells))
  }
  
  rcelln5_store <- c()
  rcelln10_store <- c()
  rcellm5_store <- c()
  rcellm10_store <- c()
  
  for (i in 1:1000) {
    samp <- dat[c("PC1", "PC2")][sample(1:nrow(dat), size=s, replace=FALSE),]
    
    x1 <- samp$PC1
    y1 <- samp$PC2
    cells5 <- cellFromXY(r5, cbind(samp$PC1, samp$PC2))
    cells10 <- cellFromXY(r10, cbind(samp$PC1, samp$PC2))
    rcelln5_store <- c(rcelln5_store, length(unique(cells5)))
    rcelln10_store <- c(rcelln10_store, length(unique(cells10)))
    rcellm5_store <- c(rcellm5_store, mean(table(cells5)))
    rcellm10_store <- c(rcellm10_store, mean(table(cells10)))
  }
  
  rcelln5_store <- sort(rcelln5_store)
  rcelln10_store <- sort(rcelln10_store)
  rcellm5_store <- sort(rcellm5_store)
  rcellm10_store <- sort(rcellm10_store)
  
  store <- rbind(store, data.frame(n=s, 
                                   celln5_f=celln5_f, celln5_e=celln5_e, celln5_fe=celln5_fe, celln5_fer=celln5_fer, celln5_lhi=celln5_lhi, celln5_cnp=celln5_cnp, 
                                   celln5_rand=rcelln5_store[500], celln5_rand_lo=rcelln5_store[50], celln5_rand_hi=rcelln5_store[950], 
                                   cellm5_f=cellm5_f, cellm5_e=cellm5_e, cellm5_fe=cellm5_fe, cellm5_fer=cellm5_fer, cellm5_lhi=cellm5_lhi, cellm5_cnp=cellm5_cnp, 
                                   cellm5_rand=rcellm5_store[500], cellm5_rand_lo=rcellm5_store[50], cellm5_rand_hi=rcellm5_store[950],
                                   
                                   celln10_f=celln10_f, celln10_e=celln10_e, celln10_fe=celln10_fe, celln10_fer=celln10_fer, celln10_lhi=celln10_lhi, celln10_cnp=celln10_cnp, 
                                   celln10_rand=rcelln10_store[500], celln10_rand_lo=rcelln10_store[50], celln10_rand_hi=rcelln10_store[950], 
                                   cellm10_f=cellm10_f, cellm10_e=cellm10_e, cellm10_fe=cellm10_fe, cellm10_fer=cellm10_fer, cellm10_lhi=cellm10_lhi, cellm10_cnp=cellm10_cnp, 
                                   cellm10_rand=rcellm10_store[500], cellm10_rand_lo=rcellm10_store[50], cellm10_rand_hi=rcellm10_store[950]))
  
}

write.csv(store, "output/occupancy.csv", row.names=FALSE)

### SUBS

png("figs/fig_sub5.png", width = 3.5, height = 4, units = 'in', res = 300)

plot(PC2 ~ PC1, dat, pch=20, cex=0.8, ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
xx <- seq(xmn, xmx, length=6)[1:5]
yy <- seq(ymx, ymn, length=6)[1:5]
xxd <- diff(xx)[1]/2
yyd <- diff(yy)[1]/2
xx <- xx + xxd
yy <- yy + yyd
xxv <- rep(xx, length(yy))
yyv <- rep(yy, each=length(xx))
cells <- unique(cellFromXY(r5, cbind(dat$PC1, dat$PC2)))
rect(xxv[cells] - xxd, yyv[cells] - yyd, xxv[cells] + xxd, yyv[cells] + yyd, border="grey")

dev.off()

png("figs/fig_sub10.png", width = 3.5, height = 4, units = 'in', res = 300)

plot(PC2 ~ PC1, dat, pch=20, cex=0.8, ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
xx <- seq(xmn, xmx, length=11)[1:10]
yy <- seq(ymx, ymn, length=11)[1:10]
xxd <- diff(xx)[1]/2
yyd <- diff(yy)[1]/2
xx <- xx + xxd
yy <- yy + yyd
xxv <- rep(xx, length(yy))
yyv <- rep(yy, each=length(xx))
cells <- unique(cellFromXY(r10, cbind(dat$PC1, dat$PC2)))
rect(xxv[cells] - xxd, yyv[cells] - yyd, xxv[cells] + xxd, yyv[cells] + yyd, border="grey")

dev.off()



