mcp_ci <- function(success, trials, alpha){
  ## Copyright (C) 2001 Frank E Harrell Jr
  ## Modified by Matthew G Johnston 2020 from binconf in the Hmisc package
  ## Distributed under GNU General Public License v2 (or later) without any warranty. See the GNU General Public License for more details.
  zcrit <-  - qnorm(alpha/2)
  z2 <- zcrit * zcrit
  mc_p <- success/trials
  cl <- (mc_p + z2/2/trials + c(-1, 1) * zcrit *
           sqrt((mc_p * (1 - mc_p) + z2/4/trials)/trials))/(1 + z2/trials)
  if(success == 1)
    cl[1] <-  - log(1 - alpha)/trials
  if(success == (trials - 1))
    cl[2] <- 1 + log(1 - alpha)/trials
  return(cbind(mc_p, lower_ci=cl[1], upper_ci=cl[2]))
}

medianBootstrap<- function(data1, data2, N=5000, alpha=0.05){
  ## Calculate observed test statistic
  mediandiff<-median(data1)-median(data2)
  ## Generate the null distribution
  boots<-replicate(N, median(sample(data1,length(data1), replace=T))-median(sample(data2,length(data2),  replace=T))-mediandiff)
  ## Count the number of at resampled observations which are at least as extreme
  above <- sum(abs(boots)>=abs(mediandiff))
  ## Calculate p value and confidence intervals
  mcp<-mcp_ci(above+1,N+1, alpha)
  return(mcp)
}