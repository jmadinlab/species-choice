# Functions

cols <- hcl.colors(3, "Temps")

xmn <- -3.17
xmx <- 3.9
ymn <- -2.23
ymx <- 2.58

choice <- function(dat, s, vars, method="nearest", trait=TRUE, oppo=FALSE) {

  subset <- dat
  dropped <- vector()
  for (i in 1:(nrow(dat)-s)) {
    if (trait) {
      if (method=="voronoi") {
        v <- voronoi.mosaic(x=subset[,'PC1'], y=subset[,'PC2'], duplicate='error')
        info <- cells(v)
        areas <- unlist(lapply(info, function(x) x$area))
        areas[is.na(areas)] <- mean(areas, na.rm=TRUE)
      } 
      if (method=="nearest") {
        # areas <- apply(nndist(subset$xr, subset$yr, k=1:3), 1, mean)
        areas <- nndist(subset$PC1, subset$PC2, k=1)
        areas <- areas / max(areas)
      }
      if (method=="hypervolume") {
        com <- matrix(TRUE, 1, length(subset$species))
        colnames(com) <- subset$species
        
        traits <- subset[c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter", "cat_skeletaldensity", "cat_colonyheight", "cat_SA_vol", "cat_spacesize")]
        rownames(traits) <- colnames(com)
        hv <- kernel.build(com, traits, abund = FALSE)
        
        areas <- kernel.contribution(hv, relative=TRUE)
      }
    } else {
      areas <- rep(1, length(subset$PC2))
    }
    
    if ("abund" %in% vars) {
      areas <- areas * subset$abund
    }
    if ("bleach" %in% vars) {
      if (oppo) {
        areas <- areas * (1 - subset$bleach)
      } else {
        areas <- areas * (subset$bleach)
      }
    }
    if ("restore" %in% vars) {
      areas <- areas * (subset$restore)
    }
    if ("range" %in% vars) {
      areas <- areas * (subset$range)
    }
    
    smallest <- which(areas == min(areas, na.rm=TRUE))[1]
    dropped <- c(dropped, which(paste(dat[,'PC1'], dat[,'PC2'], sep='_') == paste(subset[smallest,'PC1'], subset[smallest,'PC2'], sep='_')))
    subset <- subset[-smallest,]
    areas <- areas[-smallest]
  }
  temp <- dat[-dropped,]
  temp$value <- areas
  return(temp)
}


