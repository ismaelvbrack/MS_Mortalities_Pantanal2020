
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##* Estimating mortality for Medium/Large size mammals
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##* Multi-spp Zero-inflated Poisson N-mixture model
##* 423 plots of 1ha
##* Bayesian approach with JAGS

library(jagsUI)

# Import data -------------------------------------------------------------
mammInfo <- read.table(file.path("data","mammal_info.txt"),h=T)
mammals <- read.table(file.path("data","mammal_records.txt"),h=T)
plotsCovs <- read.table(file.path("data","plots_covariates.txt"),h=T)


# Organize data -----------------------------------------------------
taxa <- sort(unique(mammals$Genero))

# count matrix
Y <- as.array(table(mammals$plotID,mammals$Dupla,mammals$Genero))


dimnames(Y)[[3]] <- mammInfo[which(do.call(rbind,strsplit(mammInfo$species,"_"))[,1] %in% taxa),"sp.id"]

# Augmented array
Y.DA <- abind::abind(Y,
                     array(data=0,dim=c(nrow(Y),2,nrow(mammInfo)-length(taxa)),
                           dimnames=list(rownames(Y),colnames(Y),mammInfo[-which(do.call(rbind,strsplit(mammInfo$species,"_"))[,1] %in% taxa),"sp.id"])),
                     along=3)

Y.DA <- Y.DA[,,order(dimnames(Y.DA)[[3]])]

# number of detectec carcasses
ncap <- apply(Y.DA,c(1,3),sum)

# number of species
nspp <- nrow(mammInfo)

# Scale covariates
plotsCovs.s <- scale(plotsCovs[,c(3,6:16)])


# Arrange data for JAGS ---------------------------------------------------
#* Bundle data
my.data <- list(Y=Y.DA,ncap=ncap, # counts
                forest=plotsCovs.s[,c("forest0","forest200","forest400","forest600","forest800","forest1000")],
                dist2wiw=plotsCovs.s[,"dist2wiw"],
                deltaBurn=plotsCovs.s[,"deltaBurn"],
                tanque=plotsCovs$tanque,
                greenVeg=plotsCovs.s[,"greenVeg"],
                S=nrow(Y),n.spp=nspp) # sampling design

#* Define initial values
inits <- function(){list(w=rep(1,nspp),a=matrix(1,nrow=nrow(Y.DA),ncol=nspp),N=ncap+3,
                         mu.delta0=rnorm(1),mu.alpha0=rnorm(1),mu.beta0=rnorm(1), # intercepts
                         #omega0=rnorm(1),omega1=rnorm(1),gamma0=rnorm(1),gamma1=rnorm(1),
                         mu.delta1=rnorm(1),mu.delta2=rnorm(1),
                         alpha1=rnorm(1),alpha2=rnorm(1),
                         beta1=rnorm(1),#beta2=rnorm(1), # slopes
                         sd.delta0=runif(1,0,1),sd.alpha0=runif(1,0,1), 
                         sd.beta0=runif(1,0,1),sd.delta1=runif(1,0,1),sd.delta2=runif(1,0,1)
)}

# Define parameters of interest to be monitored (these parameters will appear on the JAGS output)
params <- c("scale_for",
            "mu.delta0","mu.delta1",#"mu.delta2",
            "mu.alpha0","alpha1","alpha2","mean.p",                
            #"mu.beta0","beta0","beta1",#"beta2",
            "alpha0","delta0","delta1",#"delta2",
            #"omega0","omega1","gamma0","gamma1",  
            "sd.delta0","sd.delta1","sd.alpha0",#"sd.beta0","sd.delta2",
            "w")

# MCMC settings
# na=100; ni=100; nt=1; nb=50; nc=3
ni <- 200000  # Iterations: number of steps for each chain
nt <- 20      # Thinning: number of steps to keep (save) every nt steps
nb <- 20000  # Burning: number of initial steps to throw out
nc <- 3       # Chains: number of chains
na <- 10000  # Adaptive: number of steps for the adaptive (pre) phase

#* Running!

out <- jags(my.data, inits, params, model=here::here("R","JAGS_models","multisppDA_fatalNmix_covars9_6scales.txt"),
            n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=T)


saveRDS(out, file.path("outputs","resu_modelCovars3Pool.rds"))

