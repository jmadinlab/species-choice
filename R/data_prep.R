# Characteristics - Coral Traits Database  - www.coraltraits.org
ctd <- read.csv("data/data_20230510.csv", as.is=TRUE)

dat <- unique(ctd[c("specie_name", "family_molecules")])
names(dat) <- c("species", "family_molecules")
dat$genus <- sapply(strsplit(dat$species, " "), "[", 1)

# range size (hughes)
rs <- ctd[ctd$trait_name=="Range size" & ctd$resource_id==605, c("specie_name", "value")]
rs$value <- as.numeric(rs$value)
names(rs) <- c("species", "range")
rs <- aggregate(range ~ species, rs, mean)
dat <- merge(dat, rs, all.x=TRUE)
dat$range <- dat$range / max(dat$range, na.rm=TRUE) # normalize

# Local abundance
ab <- ctd[ctd$trait_name=="Abundance GBR", c("specie_name", "value")]
names(ab) <- c("species", "abund")
dat <- merge(dat, ab, all.y=TRUE)
dat$abund[dat$abund=="common"] <- 1  # normalize
dat$abund[dat$abund=="uncommon"] <- 0.5 # normalize
dat$abund[dat$abund=="rare"] <- 0.25 # normalize
dat$abund <- as.numeric(dat$abund)

# Bleaching susceptibility
bri <- read.csv("data/bleaching_response_index_swain.csv")[,c(1,2)]
names(bri) <- c("species", "BRI")
BRI_genus <- bri$BRI[match(dat$genus, bri$species)]
dat$bleach <- bri$BRI[match(dat$species, bri$species)]
dat$bleach[is.na(dat$bleach)] <- BRI_genus[is.na(dat$bleach)]
dat$bleach <- 1 - (dat$bleach / 100) # normalize

# Growth form typical
gf <- ctd[ctd$trait_name=="Growth form typical", c("specie_name", "value")]
names(gf) <- c("species", "growth_form")
dat <- merge(dat, gf, all.x=TRUE)

# Life history
lh <- ctd[ctd$trait_name=="Life history strategy", c("specie_name", "value")]
names(lh) <- c("species", "life_history")
lh$life_history[lh$life_history=="weedy"] <- "Weedy"
lh$life_history[lh$life_history=="generalist"] <- "Generalist"
lh$life_history[lh$life_history=="competitive"] <- "Competitive"
lh$life_history[lh$life_history=="stress-tolerant"] <- "Stress-tolerant"
dat <- merge(dat, lh, all.x=TRUE)

# Traits - McWilliam et al. 2018 PNAS
tra <- read.csv("data/mcwilliam2018.csv", as.is=TRUE)
dim(tra[tra$domain=="pacific",])
tra <- tra[,c("species", "cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")]
dat <- merge(dat, tra)

# Goreau. The trait loadings: growth rate (GR), corallite width (CW), rugosity/branch spacing (R), surface area per unit volume (SAV), colony height (CH), maximum colony size/diameter (MCS), and skeletal density (SD). 
# builders were classified as species with the highest values of size, height, and volume; 
# fillers with largest values of size and rugosity; and 
# cementers with the largest sizes and smallest height values 
dat$goreau <- NA
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_spacesize %in% c(4, 5)] <- "Fillers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(1, 2)] <- "Cementers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(3, 4, 5) & dat$cat_SA_vol  %in% c(1, 2)] <- "Builders"

# Restoration potential from PLoS paper
dat$restore <- 0
dat$restore[dat$growth_form %in% c("branching_open", "branching_closed")] <- 6
dat$restore[dat$growth_form %in% c("massive", "submassive")] <- 5
dat$restore[dat$growth_form %in% c("laminar")] <- 4
dat$restore[dat$growth_form %in% c("columnar", "encrusting", "encrusting_long_uprights")] <- 3
dat$restore[dat$growth_form %in% c("digitate", "corymbose", "tables_or_plates")] <- 2
dat$restore[dat$growth_form %in% c("hispidose")] <- 1
dat$restore <- (dat$restore) / 6 # normalize

# Make species anonymous
dat$species_n <- paste0("Species ", sample(1:nrow(dat)))

# Coral Nurture Program

cnp <- read.csv("data/cnp.csv", as.is=TRUE)
cnp <- as.character(unique(cnp$species))
cnp[cnp == "Acropora nobilis"] <- "Acropora intermedia"
cnp[cnp == "Galaxia fasicularis"] <- "Galaxea fascicularis"
cnp[cnp == "Isopora prolifera"] <- "Isopora palifera"
cnp[cnp == "Pocillopora verucossa"] <- "Pocillopora verrucosa"
cnp[cnp == "Seriatopora calliendrum"] <- "Seriatopora caliendrum"
cnp[!(cnp %in% dat$species)]
dat$cnp <- dat$species %in% cnp

# Lord Howe
lhi <- read.csv("data/lordhowe.csv", as.is=TRUE)
lhi$species[!(lhi$species %in% dat$species)]
dat$lhi <- dat$species %in% lhi$species
sum(dat$lhi)

# export
dat <- dat[!is.na(dat$bleach),]
write.csv(dat, "output/data.csv")


