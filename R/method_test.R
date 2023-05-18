# This code looks at differences among the species selection methods (voronoi, nearest neighbor and hypervolume). SLOW

mat_hh <- matrix(NA, nrow(dat), nrow(dat)-4)
row.names(mat_hh) <- dat$species
mat_hh[,1] <- TRUE
mat_vv <- mat_nn <- mat_hh
flag <- TRUE

for (i in 1:(nrow(dat)-5)) {

  # # VORONOI
  if (flag) {
    dat_vv <- dat[mat_vv[,i],]
    v <- voronoi.mosaic(x=dat_vv[,'PC1'], y=dat_vv[,'PC2'], duplicate='error')
    info <- cells(v)
    vv <- unlist(lapply(info, function(x) x$area))
    vv[is.na(vv)] <- mean(vv, na.rm=TRUE)
    vv <- vv / max(vv)
    if (all(!is.na(vv))) {
      vvr <- which(vv == min(vv, na.rm=TRUE))
      if (length(vvr) > 1) { vvr <- sample(vvr, 1) }
      dat_vv <- dat_vv[-vvr,]
      mat_vv[,i+1] <- row.names(mat_vv) %in%  dat_vv$species
    } else { flag <- FALSE }
  }

  # NEAREST
  dat_nn <- dat[mat_nn[,i],]
  nn <- nndist(pp3(dat_nn$PC1, dat_nn$PC2, dat_nn$PC3, box3(c(0,1))), k=1)
  nn <- nn / max(nn)
  nnr <- which(nn == min(nn, na.rm=TRUE))
  if (length(nnr) > 1) { nnr <- sample(nnr, 1) }
  dat_nn <- dat_nn[-nnr,]
  mat_nn[,i+1] <- row.names(mat_nn) %in%  dat_nn$species

  # HYPERVOLUME
  dat_hh <- dat[mat_hh[,i],]
  com <- matrix(TRUE, 1, length(dat_hh$species))
  colnames(com) <- dat_hh$species
  traits <- dat_hh[c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter", "cat_skeletaldensity", "cat_colonyheight", "cat_SA_vol", "cat_spacesize")]
  rownames(traits) <- colnames(com)

  h <- kernel.build(com, traits, abund = FALSE)

ke <- kernel.evenness(h)
kec <- kernel.evenness.contribution(h)

  hh <- kernel.contribution(h, relative=TRUE)
  hh <- hh / max(hh)
  hhr <- which(hh == min(hh, na.rm=TRUE))
  if (length(hhr) > 1) { hhr <- sample(hhr, 1) }
  dat_hh <- dat_hh[-hhr,]
  mat_hh[,i+1] <- row.names(mat_hh) %in%  dat_hh$species

  write.csv(mat_vv, "output/mat_vv.csv", row.names=FALSE)
  write.csv(mat_nn, "output/mat_nn.csv", row.names=FALSE)
  write.csv(mat_hh, "output/mat_hh.csv", row.names=FALSE)

  print(i)
}


