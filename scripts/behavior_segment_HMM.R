## Segment Crested Tern data into behaviors - HMM ## --------------------------

pacman::p_load(momentuHMM, track2KBA, dplyr, sf, sp, lubridate, mapview)


trips <- readRDS("data/analysis/B_tripsplit_tracks.rds")

## custom projection ## 
# trips_proj <- projectTracks(trips, projType = "azim", custom = T)
## Projection for Guinea-Bissau (EPSG 2095)
trips_proj <- spTransform(trips, 
                  CRSobj = CRS("+proj=utm +zone=28 +ellps=intl +towgs84=-173,253,27,0,0,0,0 +units=m +no_defs"))
n_ID <- n_distinct(trips_proj$ID)

prj <- crs(trips_proj)

ncores <- 5

# use default starting values (theta) and explore likelihood surface using retryFits
tracks4fit <- trips_proj
# tracks4fit <- subset(tracks4fit, tracks4fit$ID %in% unique(tracks4fit$ID)[c(1)])

crwOut_all <- crawlWrap(obsData=tracks4fit,
                    mov.model= ~1, 
                    timeStep = "10 min",
                    ncores=ncores,
                    retryFits=50,
                    # err.model=list(x=~rep(10, nrow(track)),y=~rep(10, nrow(track)), rho=~rep(0, nrow(track))),
                    # fixPar=c(1,1,NA,NA),
                    Time.name = "DateTime",
                    attempts=50
)

plot(crwOut)

logliks <- data.frame( # inspect fits for each ind. (no outlier log-likelis. (positives are bad))
  ID     = unique(tracks4fit$ID),
  loglik = sprintf("%.0f", do.call(rbind, 
                                   lapply(crwOut$crwFits, function(x) x$loglik)
  ))
)
logliks

badfits <- logliks$ID[which(logliks$loglik > 0)]
badfits

## Fit HMM for single imputation of track (i.e. easy and fast version) ~~~~ ####
# set up data for fitting
hmmData <- prepData(data=crwOut)

hmmData <- subset(hmmData, (!hmmData$ID %in% badfits) )

## Fit 2-state model (for starters) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# nbStates <- 2
# stateNames <- c("moving", "resting") # naming states
# 
# # Set assumed distributions of 'data streams' # 
# # dist <- list(step = "gamma", angle = "wrpcauchy") 
# # dist <- list(step = "gamma", angle = "vm")
# 
# # visualize temporal auto-correlation structure
# # acf(hmmData$step[!is.na(hmmData$step)], lag.max=72)
# 
# # check out data (useful for setting intial parameters)
# plot(hmmData)
# 
# # vignette for initial paramter setting https://cran.r-project.org/web/packages/moveHMM/vignettes/moveHMM-starting-values.pdf
# # sd should be same order of mag. as mean
# # zero-mass parameter gives proportion of step lengths == 0 in data
# stepPar0 <- c(1000,5000, 4000,10000) # (mu_1, mu_2, sd_1, sd_2, zeromass_1, zeromass_2)
# anglePar0 <- c(1.5,0.3) # (mean1, mean2) ## wrp
# 
# ## model no covariates, no DM
# m0 <- fitHMM(data = hmmData, nbStates = 2,
#              dist = dist,
#              stateNames = stateNames,
#              Par0 = list(step = stepPar0, angle = anglePar0),
#              formula = ~ 1)
# plot(m0)


## Try random set of initial parameter values (2-state) ~~~~~~~~~~~~~~~~~~~~~~~
dist <- list(step = "gamma", angle = "vm")
# dist <- list(step = "gamma", angle = "wrpcauchy") 

# For reproducibility
set.seed(12345)
# Number of tries with different starting values
niter <- 10
# Save list of fitted models
allm <- list()
par0s <- list()

for(i in 1:niter) {
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(100, 3000),
                     max = c(3000, 15000))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(500, 5000),
                   max = c(5000, 20000))
  # Turning angle concentration
  # anglePar0 <- runif(2,
  #                    min = c(0.1, 0.1),
  #                    max = c(1, 5))
  anglePar0 <- runif(2, ## von mises
                     min = c(0.1, 1),
                     max = c(1, 5))
  # anglePar0 <- runif(2, ## wrapped cauchy
  #                    min = c(0.1, 0.1),
  #                    max = c(1, 5))
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0)
  allm[[i]] <- fitHMM(data = hmmData, nbStates = 2, 
                      Par0 = list(step = stepPar0, angle = anglePar0),
                      dist = dist,
                      stateNames = stateNames)
  par0s[[i]] <- list(step = stepPar0, angle = anglePar0)
  
}

allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
allnllk

# Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)
startPars <- par0s[whichbest][[1]]
# Best fitting model
mbest <- allm[[whichbest]]
mbest
plot(mbest)

# compute the pseudo-residuals (visual check of model-fit)
pr <- pseudoRes(mbest)
shapiro.test(pr$stepRes)  # both are non-normal --> indicating poor fit...
shapiro.test(pr$angleRes)
# time series, qq-plots, and ACF of the pseudo-residuals (visual check of model-fit)
plotPR(mbest)

## add state info back into tracks and map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hmmData$state <- viterbi(mbest) 

tracks_seg <- hmmData %>% mutate(
  state = ifelse( state == 2, "moving", "resting")
) %>% st_as_sf(coords = c("x", "y"), crs = prj, agr = "constant")

# mapview(tracks_seg, zcol = "state", col.regions=scales::brewer_pal("div"))
tracks_seg %>% filter(state == "resting") %>% mapview(col.regions = "red")


### Run for multiple individuals ## -------------------------------------------
