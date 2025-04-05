library(tidyverse)

# filename of parameters for simulation
parameters.fname <- "parameters.csv"

# change 'sim' to correspond to the simulation number
sim <- 1

out.fname <- paste0("data_",sim,".RData")
set.seed(sim)

# read in corresponding row of the parameters file
parameters.df <- read.csv(parameters.fname)
parameters.row <- (sim %% nrow(parameters.df)) + 1
parameters <- parameters.df[parameters.row,]

# establish hyperparameters
n.individuals.total <- parameters$n.ind
n.symptoms <- parameters$n.symp
nuisance.pct <- parameters$npct
K <- parameters$K
beta.high <- parameters$beta.high
bin.high <- parameters$bin.high
majority.ratio <- parameters$maj.ratio
tw <- parameters$tw
data.type <- parameters$data.type #cts, bin

# establish distribution parameters
high <- c(beta.high,3)
medium <- c(10,10)
low <- c(3,8)
cts.dist.params <- list(H = high, M = medium, L = low)
bin.dist.params <- list(H = bin.high, M = 0.5, L = 0.1)

# establish class memberships
etas <- rep(1,K)
etas[1] <- majority.ratio
etas <- etas/(sum(etas))
class.memberships <- sample(1:K,n.individuals.total,replace = TRUE,prob = etas)

# establish thresholds
thresholds <- runif(n.individuals.total,0.5-tw,0.5+tw)

# establish class profiles
make.class.profile <- function(i,n.symptoms,nuisance.pct){
  n.relevant <- round(n.symptoms*(1-nuisance.pct))
  v <- sample(c("H","L"),n.relevant,replace = TRUE)
  return(c(v,rep("M",n.symptoms-n.relevant)))
}

# make sure no relevant variable is the same across all conditions
ensure.relevant <- function(i,chosen.profiles){
  v <- chosen.profiles[,i]
  uv <- unique(v)
  if(length(uv) > 1){
    return(v)
  }
  if(uv == "H"){
    v[sample(1:length(v), size = 1)] <- "L"
    return(v)
  }
  if(uv == "L"){
    v[sample(1:length(v), size = 1)] <- "H"
    return(v)
  }
}


n.relevant <- round(n.symptoms*(1-nuisance.pct))
available.profiles <- expand.grid(replicate(n.relevant, c("H","L"), simplify=FALSE))
chosen.profiles <- available.profiles[sample(1:nrow(available.profiles),K),]
mod.chosen.profiles <- sapply(1:ncol(chosen.profiles),ensure.relevant,chosen.profiles)
class.profiles <- t(cbind(mod.chosen.profiles,matrix("M",nrow = K,ncol = n.symptoms-n.relevant)))


# generate data
get.cts.probs <- function(i,HL,dist.params){
  params.key <- HL[i]
  params <- dist.params[[params.key]]
  alpha <- params[1]
  beta <- params[2]
  return(rbeta(1,alpha,beta))
}

get.bin.probs <- function(i,HL,dist.params){
  params.key <- HL[i]
  p <- dist.params[[params.key]]
  return(rbinom(1,1,p))
}

generate.data <- function(i,class.memberships, thresholds, class.profiles, cts.dist.params, bin.dist.params, data.type){
  class <- class.memberships[i]
  HL <- class.profiles[,class]
  if(data.type == "cts"){
    p <- thresholds[i]
    cts <- sapply(1:length(HL),get.cts.probs,HL,cts.dist.params)
    v <- cts
    v[v < p] <- 0
    v[v > 0] <- 1
    return(v)
  }
  if(data.type == "bin"){
    v <- sapply(1:length(HL),get.bin.probs,HL,bin.dist.params)
    return(v)
  }
}
df <- data.frame(t(sapply(1:n.individuals.total, generate.data, class.memberships, thresholds, class.profiles, cts.dist.params, bin.dist.params, data.type)))


# save data
df$disease.vec <- as.character(class.memberships)
df$thresholds <- thresholds
save(df, file=out.fname)



