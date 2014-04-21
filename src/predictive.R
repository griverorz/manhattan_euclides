## Description: Simulate predicted values of the model
## Author: @griverorz
## Date: 30 Nov 2012

setwd('~/Documents/wip/manhattan_euclides/')

library(reshape)
library(lattice)
library(coda)
library(boa)
library(foreign)
library(mlogit)
library(mi)
library(mitools)
library(msm)
library(rjags)
library(ggplot2)
library(R2jags)

theme_set(theme_bw())
# blue
myblue <- rgb(100, 170, 200, max = 255)
# red
myred <- rgb(240, 140, 100, max = 255)
# green
mygreen <- rgb(70, 120, 35, max = 255)

# location parties
source("~/Documents/wip/manhattan_euclides/src/organize_data.R")

extractloc <- function(data) {
    idloc <- data[, grep("\\<ideol.", names(data))]
    nacloc <- data[, grep("\\<nacl.", names(data))]
    idloc <- apply(idloc, 2, function(x) quantile(x, .5))
    nacloc <- apply(nacloc, 2, function(x) quantile(x, .5))
    partyloc <- cbind(idloc, nacloc)
    partynames <- do.call("rbind",
                          strsplit(rownames(partyloc), ".", fixed = TRUE))[,2]
    partyloc <- partyloc[order(partynames),]
    return(partyloc)
}

rowMedians <- function(x) apply(x, 1, function(x) quantile(x, .5))
rowQMins <- function(x) apply(x, 1, function(x) quantile(x, .025))
rowQMaxs <- function(x) apply(x, 1, function(x) quantile(x, .975))

colMedians <- function(x) apply(x, 2, function(x) quantile(x, .5))

## SIMULATION OF COOPTATION AREAS GIVEN LOCATION OF PARTIES
simcoefs <- function(model, prob = .5) {
    # change parameters according to the model
    betas <- model[, grepl("beta$", dimnames(model)[[2]])]
    pbeta1 <- model[, grepl("pbeta1", dimnames(model)[[2]])]
    pbeta2 <- model[, grepl("pbeta2", dimnames(model)[[2]])]
    omegas <- model[, grepl("omega", dimnames(model)[[2]])]
    rhos <- model[, grepl("rho", dimnames(model)[[2]])]

    alphas <- vector("numeric", length(pbeta1))
    ## rnorm(length, mus[1], taus[1])
    for (i in 1:length(pbeta1)) {
        alphas[i] <- qbeta(prob, pbeta1[i], pbeta2[i])
    }
    coefs <- cbind(betas, alphas, rhos, omegas)
    return(coefs)
}

simprob <- function(coefs, partyloc, gridData) {
    nparty <- nrow(partyloc)

    probs <- matrix(NA, nrow = nrow(coefs), ncol = nparty)
    prob_vote_1 <- matrix(NA, nrow = nrow(gridData), ncol = nparty)
    prob_vote_2 <- matrix(NA, nrow = nrow(gridData), ncol = nparty)
    prob_vote_3 <- matrix(NA, nrow = nrow(gridData), ncol = nparty)

    for (loc in 1:nrow(gridData)) {
        for (p in 1:nparty) {
            distance <- abs(partyloc[p,] - gridData[loc, ])
            intercept <- rep(0, nparty - 1)
            if (p > 1) intercept[p - 1] <- 1
            probs[, p] <- coefs[, "betas"] * (
                coefs[, "alphas"] %o% distance[1] +
                (1 - coefs[, "alphas"]) %o% distance[2])^coefs[, "rhos"] +
                    coefs[,4:ncol(coefs)] %*% intercept
            ## probs[, p] <- coefs[, 1] %o% distance[1] +
            ##               coefs[, 2] %o% distance[2] +
            ##               coefs[,3:ncol(coefs)] %*% intercept
        }
        prob_vote_1[loc,] <- rowQMins(apply(probs, 1,
                                            function(x) exp(x)/sum(exp(x))))
        prob_vote_2[loc,] <- rowMedians(apply(probs, 1,
                                              function(x) exp(x)/sum(exp(x))))
        prob_vote_3[loc,] <- rowQMaxs(apply(probs, 1,
                                            function(x) exp(x)/sum(exp(x))))

    }
    return(list("qmin" = prob_vote_1,
                "median" = prob_vote_2,
                "qmax" = prob_vote_3))
}

############ APPLY TO DIFFERENT DATASETS (manually change names) ############

## FAKE DATASET OF ALL POSITIONS IN THE IDxNAC SPACE
gridData <- expand.grid(1:10, 1:10)
names(gridData) <- c("id", "nac")

load("data/mcmc_data/eus_re_full.RData")
e_estimation <- as.mcmc(fitmodel)

posEus <- na.omit(posEus)
e_partyloc <- extractloc(posEus)

e_coefs <- simcoefs(e_estimation)
prob_vote <- simprob(e_coefs, e_partyloc, gridData)

plotdata <- cbind(gridData, prob_vote$median)
## names(plotdata) <- c("id", "nac", "BNG", "PP", "PSOE")
names(plotdata) <- c("id", "nac", "Aralar", "PNV", "PP", "PSOE")
## names(plotdata) <- c("id", "nac", "CIU", "ERC", "ICV", "PP", "PSC")
plotdata <- melt(plotdata, c("id", "nac"))

pq <- ggplot(plotdata, aes(x = as.integer(id),
                           y = as.integer(nac),
                           z = value,
                           group = variable)) +
    geom_tile(aes(fill = value)) +
    facet_wrap(~ variable, nrow = 2) +
    scale_fill_gradient("Prob.", low = "white", high = "red") +
    ## labs(title = "Predicted vote \n Catalonian sample") +
    xlab("Ideological scale") +
    ylab("Nationalistic scale")
ggsave("img/pred_eus.pdf", width = par("din")[1],  pq)

##### SIMULATE EFFECT OF DISTANCE WITH RESPECT TO A GIVEN PARTY (PSOE)
simprob_party <- function(coefs, partyloc, data, party) {
    prob_vote <- simprob(coefs, partyloc, data)
    locparty <- grep(party, rownames(partyloc))
    posparty <- partyloc[grep(party, rownames(partyloc)), ]

    probdist <- sapply(prob_vote, function(x) x[, locparty])

    prob_nac_axis <- data[, "id"] == posparty["idloc"]
    prob_id_axis <- data[, "nac"] == posparty["nacloc"]

    prob_nac_axis <- probdist[prob_nac_axis, ]
    prob_id_axis <- probdist[prob_id_axis, ]

    prob_axis <- as.data.frame(rbind(prob_id_axis, prob_nac_axis))
    prob_axis$location <- rep(1:10, 2)
    prob_axis$dimension <- rep(c("id", "nac"), each = 10)
    return(prob_axis)
}


it <- 1
g_probs <- g_coefs <- vector("list", 3)
e_probs <- e_coefs <- vector("list", 3)
c_probs <- c_coefs <- vector("list", 3)
for (ci in c(.025, .5, .975)) {
    cat("Galicia: ", it, " \n")
    g_coefs[[it]] <- simcoefs(g_estimation, ci)
    g_probs[[it]] <- simprob_party(g_coefs[[it]], g_partyloc, gridData, "psdg")

    cat("Euskadi: ", it, " \n")
    e_coefs[[it]] <- simcoefs(e_estimation, ci)
    e_probs[[it]] <- simprob_party(e_coefs[[it]], e_partyloc, gridData, "pse")

    cat("Catalonia: ", it, " \n")
    c_coefs[[it]] <- simcoefs(c_estimation, ci)
    c_probs[[it]] <- simprob_party(c_coefs[[it]], c_partyloc, gridData, "psc")
    it <- it + 1
}

g_probs <- do.call(rbind, g_probs)
e_probs <- do.call(rbind, e_probs)
c_probs <- do.call(rbind, c_probs)

probs <- as.data.frame(rbind(g_probs, e_probs, c_probs))
probs$region <- rep(c("Galicia", "Basque Country", "Catalonia"), each = 20*3)
probs$voter <- rep(c("2.5% percentile", "50% percentile", "97.5% percentile"),
                   each = 20)

p <- ggplot(probs,
            aes(x = as.integer(location),
                y = median, ymin = qmin, ymax = qmax,
                colour = dimension, group = region))
pq <- p +
    geom_line(aes(group = dimension), size = .25) +
    geom_point(size = 2.5, aes(shape = dimension)) +
    facet_grid(voter ~ region) +
    scale_colour_manual(name = "Dimension", breaks = c("id", "nac"),
                        labels = c("Ideology", "Nationalism"),
                        values = c("id" = myblue, "nac" = myred)) +
    scale_shape_manual(name = "Dimension", breaks = c("id", "nac"),
                       labels = c("Ideology", "Nationalism"),
                       values = c("id" = 19, "nac" = 17)) +
    geom_linerange(size = 1.15) +
    xlab("Location") +
    ylab("Probability (voting for PSOE)")
ggsave("img/simfull.pdf", pq)

### DISTRIBUTION OF THE 0.5 OF THE REGION
dist_quant <- function(model, quantile = .5) {
    pbeta1 <- model[, grepl("pbeta1", dimnames(model)[[2]])]
    pbeta2 <- model[, grepl("pbeta2", dimnames(model)[[2]])]

    dist <- vector("numeric", nrow(model))
    for (i in 1:nrow(model)) {
        dist[i] <- pbeta(quantile, pbeta1[i], pbeta2[i])
    }
    return(dist)
}

dist_equalw <- lapply(list(g_estimation, e_estimation, c_estimation), dist_quant)
lapply(dist_equalw, function(x) quantile(x, c(.025, .5, .975)))
