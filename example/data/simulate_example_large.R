## Simulates a realistic labor market with many error terms
# Pre-requisites
rm(list = ls())

library(Matrix)
# library(pcg)
library(igraph)
library(tictoc)
library(doParallel)
library(data.table)

# Set the seed to reproduce
set.seed(5)

# Change directory to where the current script is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #works when running or sourcing the whole script 

# Set path
sourcepath <- paste0(dirname(getwd()),'\\')

##### FUNCTIONS #####
unique.index <- function(A) {
  match(A,unique(A))
}

pmax.sparse <- function(..., na.rm = FALSE) {
  
  # check that all matrices have conforming sizes
  num.rows <- unique(sapply(list(...), nrow))
  num.cols <- unique(sapply(list(...), ncol))
  stopifnot(length(num.rows) == 1)
  stopifnot(length(num.cols) == 1)
  
  cat.summary <- rbindlist(lapply(list(...), summary))
  out.summary <- cat.summary[, list(x = max(x)), by = c("i","j")]
  
  sparseMatrix(i = out.summary$i,
               j = out.summary$j,
               x = out.summary$x,
               dims = c(num.rows, num.cols))
}

#### Main function simulating data
gendata <- function(N1,N2, years=6:8, movers=5,cp=0.3) { #default values
  
  # probability of not moving is 1-p. Probability of not moving in k years
  # is (1-p)^k.  The share of non-movers is (1-p)^k. If a firm of size
  # F should have m movers, it must have F-m non-movers, i.e. (F-m)/F = (1-p)^k
  # i.e. m = F*(1-(1-p)^k), so that p = 1-(1-m/F)^(1/k)
  
  if(N2*movers >= N1) stop('more movers than firm-size')
  hazard <- 1-(1-movers*N2/N1)^(1/mean(years))
  message('moving hazard/year=',hazard)
  message('movers per firm: ',m_per_f=N1/N2*(1- (1-hazard)^mean(years)))
  
  # Simulate fixed effects for individuals and firms and sort
  eff1 <- sort(rnorm(N1), decreasing=cp < 0) # if cp is negative, there will be negative correlation between worker and firm effects
  eff2 <- sort(rnorm(N2))
  
  # Simulate number of observations for each individual
  obs <- sample(years, N1, replace=TRUE)
  id <- allid <- 1:N1 # those are movers and stayers. whole population
  
  # random draw for age -1 and normalize as in Card, Heining, Kline (2013)
  age0 <- sample(20:60,N1,replace = T)-1-40

  # Function to pick probabilities of firm for each individual,
  # i.e. make correlation between high wage workers and high wage firms
  # Note: qnorm: quantile function from prob to random variable, qnorm(0.5) = 0 ;
  #       rnorm: simulates normal random variable; 
  #       pnorm: ditribution function, pnorm(0) = 0.5 and pnorm(qnorm(x)) = x
  drawcor <- function(k,v) pnorm(rnorm(length(k), sd=qnorm(v))+qnorm(k)) # density of normal is maximized at 0, pnorm(0) = 0.5
  # for each 0<k<1, draw an x so that k and x are correlated.
  # With v=0 we have k=x, with v=1 they are uncorrelated.
  draw <- function(k,v) {
    (1-v)*drawcor(k,0.5 + 0.5*v) + v*runif(length(k)) # if v=1, random allocation
  }

  # draw firm-sizes, so that the sum equals the number of individuals
  # there are N2 firms
  firmsize <- 2+rchisq(N2,N1/N2-2)
  
  # draw initial firms for each individual
  firm <- as.integer(id/N1*N2)
  
  # initialize variables
  mover <- rep(TRUE, N1)
  DATA <- NULL
  f1 <- NULL
  f2 <- NULL
  year <- NULL
  allmovers <- NULL
  
  ### simulate mobility for different years
  for(yr in 1:max(obs)) {
    message('yr ',yr, ' move ', sum(mover))
    # to where
    # make an array of length N1 with the firm ids
    fsize <- as.integer(round(firmsize/sum(firmsize)*length(id)))
    if(sum(fsize) < length(id)) {
      fi <- sample(length(fsize),length(id)-sum(fsize))
      fsize[fi] <- fsize[fi] + 1
    }
    firmids <- unlist(sapply(seq_along(fsize), function(i) rep(i,fsize[i]))) # firmids from which to draw for movers
    firm[mover] <- firmids[as.integer(draw((id[mover]-1)/(N1-1),1-abs(cp))*length(id) + 1)]
    f1 <- c(f1,id)
    f2 <- c(f2,firm)
    
    year <- rep(yr,length(id))
    tmp <- data.frame(id=id, t=year, age=age0+yr)
    
    if (yr==1) {id_year <- tmp
    } else {id_year <- rbind(id_year,tmp)
    }
    
    # Indicator to keep observation for next year
    keep <- obs[id] > yr
    id <- id[keep]
    age0 <- age0[keep]
    
    firm <- firm[keep]
    # Indicator to define movers
    mover <- runif(length(id)) < hazard
    allmovers <- c(allmovers,id[mover])
  }

  # wages of all workers
  f1 <- factor(f1)
  f2 <- factor(f2)
  
  NT <- length(f1)
  
  eff1 <- eff1[as.numeric(levels(f1))]
  eff2 <- eff2[as.numeric(levels(f2))]
  
  levels(f1) <- 1:nlevels(f1)
  levels(f2) <- 1:nlevels(f2)
  # normalize variances
  eff1 <- sqrt(2)*eff1/sd(eff1[f1])
  eff2 <- sqrt(2)*eff2/sd(eff2[f2])

  # stack the simulation into a data.table
  dt <- data.table(id=id_year$id, firmid=f2, pe_t=eff1[f1], fe_t=eff2[f2], year=id_year$t, age=id_year$age)
  names(dt)[1] <-"id"
  setkey(dt,id,year)
  dt[ ,lagfirmid := dt[J(id,year-1),firmid]]
  
  N <- nrow(dt)
  
  # simulate errors
  resvar <- 1
  u <- matrix(rnorm(NT*1,0,sqrt(runif(NT,resvar-resvar/2,resvar+resvar/2))),NT,1)
  
  beta_age <- 0.008
  beta_age_2 <- -0.002
  
  # log wages
  dt[ ,y:= pe_t + fe_t + age*beta_age + age^2*beta_age_2 + u]
  
  structure(list(dt,hazard=hazard))
}


##### MAIN CODE #####
#### a) Parameters ####
N_id <- 5000*1000/8     # Number of unique ids when starting the simulation
N_firm <- 400*1000/8    # Number of unique firmids
N_mover_per_firm <- 3   # Number of movers per firm
N_years <- 2:5          # Number of years to simulate
v_cp <- 0.215           # Parameter governing the assortativeness of worker and firm effects. If >0 there is some degree of PAM. 
                        # Should be in [-1,1]. If 0, the allocation is random, if 1 perfect PAM and if -1 perfect NAM

# Probability to move in a given year
hazard <- 1-(1-N_mover_per_firm*N_firm/N_id)^(1/mean(N_years))

# path to save simulation output from parent directory
outpath <- paste0(sourcepath,"data/")



#### b) Simulate data ####
# data with whole population. output 2 data.tables: 1. Whole simulated sample, 2. Hazard rate
dt_list <- gendata(N_id, N_firm, years=N_years, movers=N_mover_per_firm, cp=v_cp)

# Whole sample
dt_all <- dt_list[[1]]

# Hazard rate
hazard <- matrix(0,1,2)
hazard[1] <- dt_list[[2]]
hazard[2] <- N_mover_per_firm
rm(dt_list)


#### c) Export and clear memory ####
# print tables
write.csv(dt_all,paste0(outpath,"/example_large.csv"), row.names = F)


