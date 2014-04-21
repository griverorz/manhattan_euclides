## Description: Estimation of the metric of utility functions
## Author: @griverorz
## Date: 16 Oct 2010

## READ AND FORMAT DATA

setwd('~/Documents/wip/manhattan_euclides/')

library(reshape)
library(lattice)
library(coda)
library(boa)
library(foreign)
library(mi)
library(mitools)
library(rjags)
library(R2jags)

# galicia
posGal <- read.spss('dta/Cis2796.sav', use.value.labels = TRUE, to.data.frame = T, reencode = "latin1")
names(posGal) <- tolower(names(posGal))

posGal <- data.frame(posGal$p26a,
                     posGal$p44,
                     posGal$p35,
                     posGal$p3601,
                     posGal$p3602,
                     posGal$p3603,
                     posGal$p4501,
                     posGal$p4502,
                     posGal$p4503,
                     posGal$p47,
                     posGal$p48,
                     posGal$p49a,
                     posGal$p3102,
                     posGal$p3101,
                     posGal$p3105,
                     posGal$p52)
posGal[posGal == "97" | posGal == "98" | posGal == "99"] <- NA
names(posGal) <- c("voto",
                   "ideol",
                   "nacl",
                   "nacl.pp",
                   "nacl.psdg",
                   "nacl.bng",
                   "ideol.pp",
                   "ideol.psdg",
                   "ideol.bng",
                   "genero",
                   "edad",
                   "educacion",
                   "valinc",
                   "valoposicion",
                   "valzapatero",
                   "class")
# candidato: tourino
# oposicion: feijoo
educacion <- c("menos5anos",
               "primaria",
               "eso",
               "fpmedio",
               "bachillerato",
               "fpsuperior",
               "diplomado",
               "diplomado",
               "licenciado",
               "licenciado",
               "doctor",
               NA,
               NA)
levels(posGal$educacion) <- educacion
levels(posGal$class) <- c("asalariado", "asalariado", "empresario", "autonomo", "otro", "otro", "otro", NA)
posGal$educacion <- as.numeric(posGal$educacion)
posGal$genero <- as.numeric(posGal$genero) - 1

# euskadi
posEus <- read.spss('dta/Cis2795.sav', use.value.labels = TRUE, to.data.frame = T, reencode = "latin1")
names(posEus) <- tolower(names(posEus))

posEus <- data.frame(posEus$p26a,
                     posEus$p48,
                     posEus$p38,
                     posEus$p3901,
                     posEus$p3902,
                     posEus$p3903,
                     posEus$p3904,
                     posEus$p4901,
                     posEus$p4902,
                     posEus$p4903,
                     posEus$p4904,
                     posEus$p51,
                     posEus$p52,
                     posEus$p53a,
                     posEus$p3403,
                     posEus$p3404,
                     posEus$p3409,
                     posEus$p56)
posEus[posEus == "97" | posEus == "98" | posEus == "99"] <- NA
names(posEus) <- c("voto",
                   "ideol",
                   "nacl",
                   "nacl.pnv",
                   "nacl.pse",
                   "nacl.pp",
                   "nacl.aralar",
                   "ideol.pnv",
                   "ideol.pse",
                   "ideol.pp",
                   "ideol.aralar",
                   "genero",
                   "edad",
                   "educacion",
                   "valinc",
                   "valoposicion",
                   "valzapatero",
                   "class")
# candidato: ibarretxe
# oposicion: patxi lopez
levels(posEus$educacion) <- educacion
levels(posEus$class) <- c("asalariado", "asalariado", "empresario", "autonomo", "otro", "otro", "otro", NA)
posEus$educacion <- as.numeric(posEus$educacion)
posEus$genero <- as.numeric(posEus$genero) - 1

# cataluna
posCat <- read.spss('dta/Cis2660.sav', use.value.labels = TRUE, to.data.frame = T, reencode = "latin1")
names(posCat) <- tolower(names(posCat))

posCat <- data.frame(posCat$p15a,
                     posCat$p31,
                     posCat$p29,
                     posCat$p3001,
                     posCat$p3002,
                     posCat$p3003,
                     posCat$p3004,
                     posCat$p3005,
                     posCat$p3201,
                     posCat$p3202,
                     posCat$p3203,
                     posCat$p3204,
                     posCat$p3205,
                     posCat$p35,
                     posCat$p36,
                     posCat$p37a,
                     posCat$p2203,
                     posCat$p2202,
                     posCat$p2301,
                     posCat$p42)
# candidato: montilla
# oposicion: mas
posCat[posCat == "98" | posCat == "99"] <- NA
names(posCat) <- c("voto",
                   "ideol",
                   "nacl",
                   "nacl.ciu",
                   "nacl.psc",
                   "nacl.erc",
                   "nacl.pp",
                   "nacl.icv",
                   "ideol.ciu",
                   "ideol.psc",
                   "ideol.erc",
                   "ideol.pp",
                   "ideol.icv",
                   "genero",
                   "edad",
                   "educacion",
                   "valinc",
                   "valoposicion",
                   "valzapatero",
                   "class")

levels(posCat$educacion) <- educacion
posCat$educacion <- as.numeric(posCat$educacion)
levels(posCat$class) <- c("asalariado", "asalariado", "empresario", "autonomo", "otro", "otro", "otro", NA)
posCat$genero <- as.numeric(posCat$genero) - 1

### final vote
# galicia
posGal$voto <- as.character(posGal$voto)
posGal$voto[posGal$voto != "PP" &
            posGal$voto != "PSdG" &
            posGal$voto != "BNG"] <- NA
posGal$voto <- as.factor(posGal$voto)

# euskadi
posEus$voto <- as.character(posEus$voto)
posEus$voto[posEus$voto == "PSE EE"] <- "PSE"
posEus$voto[posEus$voto != "PNV" &
            posEus$voto != "PSE" &
            posEus$voto != "PP" &
            posEus$voto != "Aralar"] <- NA
posEus$voto <- as.factor(posEus$voto)

# cataluna
posCat$voto <- as.character(posCat$voto)
posCat$voto[posCat$voto != "CiU" &
            posCat$voto != "PSC" &
            posCat$voto != "ERC" &
            posCat$voto != "PP" &
            posCat$voto != "ICV"] <- NA
posCat$voto <- as.factor(posCat$voto)


## grafico ubicacion
toGrid <- function(range, parties, dataset) {
  grid <- expand.grid(ideol = range, nacl = range, voto = toupper(parties))
  tab <- table(dataset[, "ideol"], dataset[, "nacl"], dataset[, "voto"])
  stab <- apply(tab, 3, sum)
  nparties <- dim(tab)[3]
  out <- vector("list", nparties)

  for (i in 1:nparties) {
    out[[i]] <- tab[, , i]/stab[i]
  }
  grid[,'freq'] <- as.vector(unlist(out))
  return(grid)
}

## Change labels and grid with region
grid <- toGrid(1:10, c("bng", "pp", "psdg"), posGal)
grid <- toGrid(1:10, c("ciu", "erc", "icv", "pp", "psc"), posCat)
grid <- toGrid(1:10, c("aralar", "pnv", "pp", "psoe"), posEus)

p <- ggplot(grid, aes(x = ideol, y = nacl, z = freq, group = voto))
pq <- p + geom_tile(aes(fill = freq)) + scale_fill_gradient("Freq.", low = "white", high = "red") +
  facet_wrap(~ voto) +
  theme_bw() + coord_equal(ratio = 1) +
  opts(title = "") + xlab("Ideological scale") + ylab("Nationalistic scale")
ggsave("paper/plots/autoub_cat.pdf", width = par("din")[1]*1.5, pq)

