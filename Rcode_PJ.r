#################################################################################
#################################################################################
#################################################################################
### Datasets:
###   1) species (x) traits/modalities 'matrix'; 108 species (x) 23 modalities (6 total traits)
###          Traits            Modalities
###          -------------------------------
###          Trophic level     2-2.5, 2.6-3.5, 3.6-4.5, 4.6-5
###          Feeding type      grazing, ambusher, active predator
###          Diet              herbivore, zoobenthivore, detritus, omnivore, zooplankton, piscivore
###          Spawning habitat  pelagic, benthic
###          Maximum length    <20cm, 20-50cm, 50-100cm, >100cm
###          Spawning period   spring, summer, autumn, winter
###   2) randomized MEDITS haul (x) species abundance 'data.frame'; 471 hauls (x) (108 species + 12 variables of haul info)
###   3) 'data.frame' storing intra-specific assymetric competition indices for all species
###   4) 'data.frame' storing all inter-specific assymetric competition indices for all species

data1 = readRDS("~/data1.rds")
data2 = readRDS("~/data2.rds")
data3 = readRDS("~/data3.rds")
data4 = readRDS("~/data4.rds")


### Calculate assymetric competition index ###
## column numbers containing haul info
index_of_non_species = 1:12

## compute competitions and store in 'list' of length 471 (one 108 x 108 matrix for each haul)
all_comps = vector("list", nrow(data2))
for(i in 1:nrow(data2)) {
  h = data2[i, ]
  index = match(colnames(h), rownames(data1))
  index = index[-index_of_non_species]
  m1 = data1[index, ]
  m2 = unlist(h[-index_of_non_species])
  m2[is.na(m2)] = 0
  mat = m1 * m2
  mat = t(mat)
  mat.q = mat / (colSums(mat)[col(mat)]); mat.q[!is.finite(mat.q)] = 0
  mat.p = mat / rowSums(mat); mat.p[!is.finite(mat.p)] = 0
  mat.comp = crossprod(mat.p, mat.q)
  
  # eliminate species not in haul
  no_catch_species = names(h[-index_of_non_species])[is.na(unlist(h[-index_of_non_species]))]
  fill_na = rownames(mat.comp) %in% no_catch_species
  mat.comp[fill_na, ] = NA
  mat.comp[, fill_na] = NA
  all_comps[[i]] = mat.comp
}
#################################################################################
#################################################################################
#################################################################################
###  Model intra-specific competitions
library(mgcv) 
library(ggplot2)

## model selected species 
csp = c("Arge sph", "Aris fol", "Cent gra", "Dasy pas", "Epin aen", "Etmo spi", "Gale mel", 
        "Hopl med", "Ille coi", "Loph bud", "Merl mer", "Mull bar", "Mura hel", "Poly ame",
        "Raja cla", "Scyl can", "Siga lur", "Spic sma", "Squa aca", "Squa bla", "Syno sau")
sub.data3 = data3[data3$species %in% csp, ]

mod = gam(intracomp ~ 
            s(Depth, k = 5) + 
            s(Depth, species, bs = "fs", k = 5) +
            s(species, bs = "re") +
            s(fyear, bs = "re") + 
            s(sh_lon, sh_lat, k = 10) + 
            s(Bottom.type.new, bs = "re") + 
            s(Bottom.type.new, species, bs = "re"),
          data = sub.data3,
          family = Gamma(link = "log"), 
          method = "REML",
          select = TRUE)


## predict/plot partial depth effect on intra-specific competition
newdata = expand.grid(Depth = seq(min(sub.data3$Depth, na.rm = TRUE), max(sub.data3$Depth, na.rm = TRUE), by = 10),
                      species = csp)
pred = predict(mod, newdata,
               type = "terms", terms = "s(Depth,species)", 
               newdata.guaranteed = TRUE, se = TRUE)
pred = data.frame(species = newdata$species, depth = newdata$Depth, 
                  fit = pred$fit[, 1], se.fit = pred$se.fit[, 1])
ggplot(data = pred) + 
  facet_wrap(~ species) + 
  geom_line(aes(x = depth, y = fit)) +
  geom_line(aes(x = depth, y = fit + se.fit), lty = 2) + 
  geom_line(aes(x = depth, y = fit - se.fit), lty = 2)

## predict/plot partial substrate effect on intra-specific competition
newdata = expand.grid(Bottom.type.new = levels(sub.data3$Bottom.type.new),
                      species = csp)
pred = predict(mod, newdata,
               type = "terms", terms = "s(Bottom.type.new,species)", 
               newdata.guaranteed = TRUE, se = TRUE)
pred = data.frame(species = newdata$species, bottom = newdata$Bottom.type.new, 
                  fit = pred$fit[, 1], se.fit = pred$se.fit[, 1])
ggplot(data = pred) + 
  facet_wrap(~ species) + 
  geom_point(aes(x = bottom, y = fit)) +
  geom_segment(aes(x = bottom, y = fit + se.fit, xend = bottom, yend = fit - se.fit))
#################################################################################
#################################################################################
#################################################################################
###  Model inter-specific competitions
library(DirichletReg)
library(splines)
library(mgcv)
library(ggplot2)
library(dplyr)

## choose species to model interspecific competition
sp = "Mull bar"
sub.data4 = data4[data4$species1 == sp, ]

fac.cols = c(1:15) ## variables holding haul info

## prepare dataset for dirichlet regression
jcomp = names(data4)[-fac.cols]
ikeep = rowSums(is.na(sub.data4[, jcomp])) != length(jcomp)
if(sum(ikeep) == 1) stop("cannot run model with 1 row")
sub.data4 = sub.data4[ikeep, ]
sub.data4zero = sub.data4
sub.data4zero[is.na(sub.data4zero)] = 0
pct = sub.data4zero[, jcomp] * 100
ddir = sub.data4zero[, fac.cols]
ddir$Y = DR_data(pct)

## only-intercept model
m0 = DirichReg(Y ~ 1, data = ddir, model = "common") 
## + depth
m1 = DirichReg(Y ~ bs(Depth, df = 4), data = ddir, model = "common")
## + substrate
m2 = DirichReg(Y ~ bs(Depth, df = 4) + Bottom.type.new, data = ddir, model = "common") 

## model comparison; select model with depth + bottom
sig = anova(m0, m1, m2)
mod = m2

## predict/plot depth*bottom effect on interspecific competition
newdata = expand.grid(Depth = seq(min(ddir$Depth, na.rm = TRUE), max(ddir$Depth, na.rm = TRUE), by = 10),
                      Bottom.type.new = levels(droplevels(ddir$Bottom.type.new)))
pred = predict(mod, newdata, mu = TRUE)
colnames(pred) = colnames(ddir$Y)
pred = stack(as.data.frame(pred))
pred = data.frame(newdata, pred); names(pred) = c("depth", "bottom", "COMP", "species")

## select, e.g., top 4 species by substrate for visualization
top_sp = as.character(unique(slice_max(group_by(aggregate(COMP ~ bottom + species, pred, min), bottom), order_by = COMP, n = 4)$species))

ggplot(data = subset(pred, species %in% top_sp),
       aes(x = depth, y = COMP, group = species, colour = species)) + 
  facet_wrap(~ bottom) + 
  geom_line()
