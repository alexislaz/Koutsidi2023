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
###   2) **randomized** (on catch abundance) MEDITS haul (x) species abundance 'data.frame'; 471 hauls (x) (108 species + 12 variables of haul info)
###   3) 'data.frame' storing intra-specific assymetric competition indices for all species
###   4) 'data.frame' storing all inter-specific assymetric competition indices for all species

#################################################################################
#################################################################################
#################################################################################
### Load datasets  (takes some time)
data1 = as.matrix(read.csv("https://raw.githubusercontent.com/alexislaz/Koutsidi2023/e6650ae4a8f3ae7f252808044fdc6ab682f60cd5/data1.csv", 
                           row.names = 1, check.names = FALSE))
data2 = read.csv("https://raw.githubusercontent.com/alexislaz/Koutsidi2023/e6650ae4a8f3ae7f252808044fdc6ab682f60cd5/data2.csv",
                 row.names = 1, check.names = FALSE)
data3 = read.csv("https://raw.githubusercontent.com/alexislaz/Koutsidi2023/e6650ae4a8f3ae7f252808044fdc6ab682f60cd5/data3.csv",
                row.names = 1, check.names = FALSE, 
                 colClasses = c(species = "factor", fyear = "factor", Bottom.type.new = "factor"))
data4 = read.csv("https://raw.githubusercontent.com/alexislaz/Koutsidi2023/e6650ae4a8f3ae7f252808044fdc6ab682f60cd5/data4.csv",
                row.names = 1, check.names = FALSE,
                colClasses = c(Bottom.type.new = "factor"))

#################################################################################
#################################################################################
#################################################################################
### Calculate assymetric competition index ###

## column numbers containing haul info variables
id_var = 1:12

## compute competitions and store in 'list' of length 471 (one 108 x 108 matrix for each haul)
all_comps = vector("list", nrow(data2))
names(all_comps) = sprintf("Haul%d", data2$aa)

for(i in 1:nrow(data2)) {
  h = data2[i, ]
  index = match(colnames(h), rownames(data1))
  index = index[-id_var]
  m1 = data1[index, ]
  m2 = unlist(h[-id_var])
  m2[is.na(m2)] = 0
  mat = m1 * m2
  mat = t(mat)
  mat.q = mat / (colSums(mat)[col(mat)]); mat.q[!is.finite(mat.q)] = 0
  mat.p = mat / rowSums(mat); mat.p[!is.finite(mat.p)] = 0
  mat.comp = crossprod(mat.p, mat.q)
  
  # eliminate species not in haul
  no_catch_species = names(h[-id_var])[is.na(unlist(h[-id_var]))]
  fill_na = rownames(mat.comp) %in% no_catch_species
  mat.comp[fill_na, ] = NA
  mat.comp[, fill_na] = NA
  all_comps[[i]] = mat.comp
}

str(all_comps)

#################################################################################
#################################################################################
#################################################################################
###  Model intra-specific competitions
library(mgcv) 
library(ggplot2)

## model selected species only 
sp.sel = c("Arge sph", "Aris fol", "Cent gra", "Dasy pas", "Epin aen", "Etmo spi", "Gale mel",
           "Hopl med", "Ille coi", "Loph bud", "Merl mer", "Mull bar", "Mura hel", "Poly ame",
           "Raja cla", "Scyl can", "Siga lur", "Spic sma", "Squa aca", "Squa bla", "Syno sau")
sub.data3 = data3[data3$species %in% sp.sel, ]

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
                      species = sp.sel)
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
                      species = sp.sel)
pred = predict(mod, newdata,
               type = "terms", terms = "s(Bottom.type.new,species)", 
               newdata.guaranteed = TRUE, se = TRUE)
pred = data.frame(species = newdata$species, substrate = newdata$Bottom.type.new, 
                  fit = pred$fit[, 1], se.fit = pred$se.fit[, 1])
ggplot(data = pred) + 
  facet_wrap(~ species) + 
  geom_point(aes(x = substrate, y = fit)) +
  geom_segment(aes(x = substrate, y = fit + se.fit, xend = substrate, yend = fit - se.fit))

#################################################################################
#################################################################################
#################################################################################
###  Model inter-specific competitions
library(DirichletReg)
library(splines)
library(ggplot2)

## choose species to model interspecific competition
sp = "Spic sma"
sub.data4 = data4[data4$species1 == sp, ]

id_var = 1:15 ## variables holding haul info

## prepare dataset for dirichlet regression
jcomp = names(data4)[-id_var]
ikeep = rowSums(is.na(sub.data4[, jcomp])) != length(jcomp)  # keep hauls with at least on competitor species
if(sum(ikeep) == 1) stop("cannot run model with 1 haul")
sub.data4 = sub.data4[ikeep, ]
sub.data4zero = sub.data4
sub.data4zero[is.na(sub.data4zero)] = 0
pct = sub.data4zero[, jcomp] * 100
ddir = sub.data4zero[, id_var]
ddir$Y = DR_data(pct)  # ignore warnings

## Run models -- may take some time for all
# only-intercept model
m0 = DirichReg(Y ~ 1, data = ddir, model = "common") 
# + depth
m1 = DirichReg(Y ~ bs(Depth, df = 4), data = ddir, model = "common")
# + substrate
m2 = DirichReg(Y ~ bs(Depth, df = 4) + Bottom.type.new, data = ddir, model = "common") 

## model comparison; select model with depth + bottom
sig = anova(m0, m1, m2)
sig
mod = m2

## predict/plot depth*bottom effect on interspecific competition
newdata = expand.grid(Depth = seq(min(ddir$Depth, na.rm = TRUE), 
                                  max(ddir$Depth, na.rm = TRUE), by = 10),
                      Bottom.type.new = levels(droplevels(ddir$Bottom.type.new)))
pred = predict(mod, newdata, mu = TRUE)
colnames(pred) = colnames(ddir$Y)
pred = stack(as.data.frame(pred))
pred = data.frame(newdata, pred); names(pred) = c("depth", "substrate", "COMP", "species")

## select, e.g., top 4 competitor species by substrate for easier visualization
ntop = 4
medCOMPs = aggregate(COMP ~ substrate + species, pred, median) # find median competition index by species across depths
ntop_sp = Reduce(union, lapply(split(medCOMPs, medCOMPs$substrate), function(d) d[order(d$COMP, decreasing = TRUE), "species"][1:ntop]))                                  

ggplot(data = subset(pred, species %in% ntop_sp),
       aes(x = depth, y = COMP, group = species, colour = species)) + 
  facet_wrap(~ substrate) + 
  geom_line()
