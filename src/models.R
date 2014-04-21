## Description: Estimate models for Manhattan Euclides paper
## Author: @griverorz
## Date: 23 Oct 2010

library(doMC); registerDoMC(3)
library(rjags); load.module("glm"); load.module("lecuyer")
library(dclone)
library(RCurl)

source("~/Documents/wip/manhattan_euclides/src/organize_data.R")

## Test
runif(1)

## GALICIA
posGal <- na.omit(posGal)
ideol <- as.matrix(abs(posGal[, "ideol"] - posGal[, c("ideol.bng",
                                                       "ideol.pp",
                                                       "ideol.psdg")]))
nacl <- as.matrix(abs(posGal[, "nacl"] - posGal[, c("nacl.bng",
                                                    "nacl.pp",
                                                    "nacl.psdg")]))
voto <- factor(posGal$voto)
voto <- model.matrix( ~ voto - 1)

galicia.data <- list("voto" = voto,
                     "ideol" = ideol,
                     "nacl" = nacl,
                     "K" = ncol(ideol),
                     "N" = nrow(ideol))

## EUSKADI
posEus <- na.omit(posEus)
ideol <- as.matrix(abs(posEus[, "ideol"] - posEus[, c("ideol.aralar",
                                                       "ideol.pnv",
                                                       "ideol.pp",
                                                       "ideol.pse")]))
nacl <- as.matrix(abs(posEus[, "nacl"] - posEus[, c("nacl.aralar",
                                                    "nacl.pnv",
                                                    "nacl.pp",
                                                    "nacl.pse")]))
voto <- factor(posEus$voto)
voto <- model.matrix( ~ voto - 1)

euskadi.data <- list("voto" = voto,
                     "ideol" = ideol,
                     "nacl" = nacl,
                     "K" = ncol(ideol),
                     "N" = nrow(ideol))

## CATALONIA
posCat <- na.omit(posCat)
ideol <- as.matrix(abs(posCat[, "ideol"] -  posCat[, c("ideol.ciu",
                                                       "ideol.erc",
                                                       "ideol.icv",
                                                       "ideol.pp",
                                                       "ideol.psc")]))
nacl <- as.matrix(abs(posCat[, "nacl"] - posCat[, c("nacl.ciu",
                                                    "nacl.erc",
                                                    "nacl.icv",
                                                    "nacl.pp",
                                                    "nacl.psc")]))
voto <- factor(posCat$voto)
voto <- model.matrix( ~ voto - 1)

catalonia.data <- list("voto" = voto,
                       "ideol" = ideol,
                       "nacl" = nacl,
                       "K" = ncol(ideol),
                       "N" = nrow(ideol))

## MODELS

## Which model and which parameters are going to be estimated
mods <- list(#"manhattan",
             ## "penalized",
             ## "euclides",
             ## "fullnopenalized",
             "fullbern_nopenalized"
             ## "fullbern",
             ## "full"
             )

params <- list(
    ## c("beta", "omega"),
    ## c("beta", "omega"),
    ## c("alpha", "beta", "omega"),
    ## c("alpha", "beta", "omega", "rho"),
    c("alpha", "beta", "omega", "rho")
    ## c("alpha", "beta", "omega", "rho", "gamma"),
    ## c("alpha", "beta", "omega", "rho", "gamma")
    )

dtas <- list(galicia.data,
             euskadi.data,
             catalonia.data)

#################### MODELS WITHOUT RANDOM EFFECTS ####################

# PARAMETERS
niter = 50000
nburn = 5000
nthin = 5

for (m in 1:length(mods)) {
    for (d in 1:length(dtas)) {
        regname <- c("gal", "eus", "cat")[d]
        running <- paste(regname, "_", mods[[m]], ".RData", sep = "")
        cat("Now running: ", running, "\n")
        setwd("~/Documents/wip/manhattan_euclides/scripts/jags_scripts/no_re/")
        fitmodel <- jags(model.file = paste(mods[[m]], ".txt", sep = ""),
                         parameters.to.save = params[[m]],
                         data = dtas[[d]],
                         n.iter = niter,
                         n.burnin = nburn,
                         n.chains = 1,
                         n.thin = nthin,
                         DIC = TRUE)
        setwd("~/Documents/wip/manhattan_euclides/data/mcmc_data/")
        save(fitmodel, file = running)
        rm(fitmodel)
    }
}

#################### MODELS WITH RANDOM EFFECTS ####################
niter = 200000
nburn = 50000
nthin = 10

mods <- list(
    ## "manhattan",
    ## "manhattan_simplified",
    ## "penalized",
    "full",
    "full_rho_nogamma",
    "full_rho_gamma"
    )

params <- list(
    ## c("mu", "tau", "omega"),
    ## c("mu", "tau", "omega", "beta"),
    ## c("mu", "tau", "omega"),
    c("pbeta1", "pbeta2", "rho", "beta", "omega"),
    c("pbeta1", "pbeta2", "pi", "beta", "omega"),
    c("pbeta1", "pbeta2", "pi", "lambda", "beta", "omega")
    )

dtas <- list(galicia.data,
             euskadi.data,
             catalonia.data)

## USING FOREACH
for (m in 1:length(mods)) {
    for (d in 1:length(dtas)) {
        regname <- c("gal", "eus", "cat")[d]
        running <- paste(regname, "_re_", mods[[m]], ".RData", sep = "")
        cat("Now running: ", running, "\n")
        setwd("~/Documents/wip/manhattan_euclides/src/jags/re/")

        fitmodel <- jags(model = paste(mods[[m]], ".txt", sep = ""),
                         data = dtas[[d]],
                         n.burnin = nburn,
                         n.chains = 1,
                         n.iter = niter,
                         n.thin = nthin,
                         parameters.to.save = params[[m]],
                         DIC = TRUE)
        setwd("~/Documents/wip/manhattan_euclides/data/mcmc_data/")
        save(fitmodel, file = running)
     }
}
