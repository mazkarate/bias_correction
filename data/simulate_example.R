## Simulates a realistic labor market without the need of a 'main' file
# Pre-requisites

rm(list = ls())
# rm(list = setdiff(ls(),c('params0','params','sourcepath')))


library(Matrix)
# library(pcg)
library(igraph)
library(tictoc)
library(doParallel)
library(data.table)

# Set the seed to reproduce
set.seed(5)

#Change directory to where the current script is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #works when running or sourcing the whole script

sourcepath <- getwd()

# FUNCTIONS

unique.index <- function(A) {
  match(A,unique(A))
}

pmax.sparse <- function(..., na.rm = FALSE) {
  
  # check that all matrices have conforming sizes
  num.rows <- unique(sapply(list(...), nrow))
  num.cols <- unique(sapply(list(...), ncol))
  stopifnot(length(num.rows) == 1)
  stopifnot(length(num.cols) == 1)
  
  # cat.summary <- do.call(rbind, lapply(list(...), summary))
  # out.summary <- aggregate(x ~ i + j, data = cat.summary, max, na.rm)
  cat.summary <- rbindlist(lapply(list(...), summary))
  out.summary <- cat.summary[, list(x = max(x)), by = c("i","j")]
  
  sparseMatrix(i = out.summary$i,
               j = out.summary$j,
               x = out.summary$x,
               dims = c(num.rows, num.cols))
}



gendata <- function(N1,N2, years=6:8, movers=5,cp=0.3) { #default values
  
  # probability of not moving is 1-p. Probability of not moving in k years
  # is (1-p)^k.  The share of non-movers is (1-p)^k. If a firm of size
  # F should have m movers, it must have F-m non-movers, i.e. (F-m)/F = (1-p)^k
  # i.e. m = F*(1-(1-p)^k), so that p = 1-(1-m/F)^(1/k)
  
  if(N2*movers >= N1) stop('more movers than firm-size')
  hazard <- 1-(1-movers*N2/N1)^(1/mean(years))
  message('moving hazard/year=',hazard)
  message('movers per firm: ',m_per_f=N1/N2*(1- (1-hazard)^mean(years)))
  # pick fixed effects for individuals and firms
  eff1 <- sort(rnorm(N1), decreasing=cp < 0)
  eff2 <- sort(rnorm(N2))
  
  # pick number of observations for each individual
  # should we make the observation length dependent on the fixed effect?
  obs <- sample(years, N1, replace=TRUE)
  #  obs <- obs - as.integer(sel*eff1)
  # pick probabilities of firm for each individual,
  # i.e. make correlation between high wage workers and high wage firms
  id <- allid <- 1:N1 # those are movers and stayers. whole population

  # random draw for age -1 and normalize as in Card, Heining, Kline (2013)
  age0 <- sample(20:60,N1,replace = T)-1-40
  

  
  
  drawcor <- function(k,v) pnorm(rnorm(length(k), sd=qnorm(v))+qnorm(k))
  # for each 0<k<1, draw an x so that k and x are correlated.
  # With v=0 we have k=x, with v=1 they are uncorrelated.
  draw <- function(k,v) {
    (1-v)*drawcor(k,0.5 + 0.5*v) + v*runif(length(k))
  }
  normal <- function(x) if(length(x)==0) numeric(0) else (x-min(x)+1e-7)/(max(x) - min(x) + 1e-7)
  # draw firm-sizes, so that the sum equals the number of individuals
  # there are N2 firms
  firmsize <- 2+rchisq(N2,N1/N2-2)
  # draw initial firms for each individual
  
  firm <- as.integer(id/N1*N2)
  mover <- rep(TRUE, N1)
  mover_occup <- rep(TRUE, N1)
  DATA <- NULL
  f1 <- NULL
  f2 <- NULL
  f3 <- NULL
  year <- NULL
  N_occup <- 20 # different occupation identifiers
  occup <- NULL
  
  allmovers <- NULL
  for(yr in 1:max(obs)) {
    message('yr ',yr, ' move ', sum(mover))
    # to where
    # make an array of length N1 with the firm ids
    fsize <- as.integer(round(firmsize/sum(firmsize)*length(id)))
    if(sum(fsize) < length(id)) {
      fi <- sample(length(fsize),length(id)-sum(fsize))
      fsize[fi] <- fsize[fi] + 1
    }
    firmids <- unlist(sapply(seq_along(fsize), function(i) rep(i,fsize[i])))
    firm[mover] <- firmids[as.integer(draw((id[mover]-1)/(N1-1),1-abs(cp))*length(id) + 1)]
    f1 <- c(f1,id)
    f2 <- c(f2,firm)
    
    year <- rep(yr,length(id))
    occup[mover_occup] <- sample(1:N_occup,sum(mover_occup),replace = T) # update occupation of those who can change
    f3 <- c(f3,occup)
    tmp <- data.frame(id=id,t=year,age=age0+yr,occup=occup)

    
    
    if (yr==1) {id_year <- tmp
    
    } else {id_year <- rbind(id_year,tmp)
    }
    
    # which to keep next iter
    keep <- obs[id] > yr
    id <- id[keep]
    age0 <- age0[keep]
    occup <- occup[keep]
    
    firm <- firm[keep]
    # which should change job?
    # should high eff1 people be more prone to change job?
    # or the product?
    #    mover <- runif(length(id)) < p*(1+2*sel*(pnorm((eff1[id]*eff2[firm]))-0.5))
    #    mover <- runif(length(id)) < p*(1+ 2*pnorm((normal(eff1[id]+eff2[firm])-0.5), sd=sel))
    mover <- runif(length(id)) < hazard
    mover_occup <- runif(length(occup)) < hazard
    allmovers <- c(allmovers,id[mover])
    
  }
  #  cat("5num eff1 of movers ", fivenum(eff1[unique(allmovers)]), '\n')
  
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
  
  # occupation fixed effects
  beta_occ <- rnorm(N_occup,0,1)

  
  dt <- data.table(id=id_year$id,firmid=f2,occ_id=f3,pe_t=eff1[f1],fe_t=eff2[f2],occ_t=beta_occ[f3],year=id_year$t,age=id_year$age)
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
  dt[ ,y:=pe_t + fe_t + occ_t + age*beta_age + age^2*beta_age_2 + u]
  
  
  

  structure(list(dt,hazard=hazard))

}




##### MAIN CODE #####
#### a) Parameters ####
n_id <- 5000*1000/8 # to give about 2.5 million observations
n_firm <- 400*1000/8

# n_id <- 5000 # to give about 22K observations
# n_firm <- 400

v_move <- 3

v_cp <- 0.215
resvar <- 1
hetero <- 1


outpath <- sourcepath


#### b) Generate data ####

cat("\n")
# data with whole population. output data.table
dt_list <- gendata(n_id,n_firm,years=2:7, movers=v_move, cp=v_cp)
hazard <- matrix(0,1,2)
# merge and filter data.table
# dt_con <- merge(x=dt_list[[1]],y=dt_list[[2]],by.x=c("id","firmid"),by.y=c("id_old","firmid_old"))
dt_all <- dt_list[[1]]
hazard[1] <- dt_list[[2]]
hazard[2] <- v_move
rm(dt_list)


# dt_all_out <- data.table(id=dt_all$id,year=dt_all$year,firmid=dt_all$firmid,pe_t=dt_all$pe_t,fe_t=dt_all$fe_t,idmatch=dt_all$idmatch)
# cat('nrow all:', nrow(dt_all_out))
# cat("\n")



# print tables
write.csv(dt_all,paste0(outpath,"/example.csv"), row.names = F)
# write.csv(dt_all_out,paste0(outpath,i,"graph.csv"), row.names = F)
write.csv(hazard,paste0(outpath,"/example_hazard.csv"), row.names = F)


 



