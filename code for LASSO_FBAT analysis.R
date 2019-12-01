## 3 supporting functions for LASSO_FBAT analysis
## 2 supporting package are needed: MASS && glmnet
## If you haven't install the needed packages, please install first.

library(glmnet)
library(MASS)

##calculate estimated covariance matrix
cal_va=function(sample,genenum)
{
  va=(sample[1:genenum,,3]-sample[1:genenum,,4])%*%t(sample[1:genenum,,3]-sample[1:genenum,,4])
  return(va)
}


##calculate statistic FBAT-MM
Stat_chi=function(sample,va,genenum)
{
  u=rep(0,genenum)
  for(i in 1:genenum)
    u[i]=sum(sample[i,,3]-sample[i,,4])
  Kai=u%*%ginv(va)%*%u
  return(Kai)
}


##calculate number ovserved
num=function(sample,genenum,sampnum)
{
  n=matrix(0,nrow=genenum,ncol=7)
  for(i in 1:genenum)
  {
    for(j in 1:sampnum)
    {
      if((sample[i,j,1]==2 & sample[i,j,2]==1 & sample[i,j,3]==2)||(sample[i,j,1]==1 & sample[i,j,2]==2 & sample[i,j,3]==2)) n[i,1]=n[i,1]+1
      if((sample[i,j,1]==2 & sample[i,j,2]==1 & sample[i,j,3]==1)||(sample[i,j,1]==1 & sample[i,j,2]==2 & sample[i,j,3]==1)) n[i,2]=n[i,2]+1
      if((sample[i,j,1]==1 & sample[i,j,2]==1 & sample[i,j,3]==2)) n[i,3]=n[i,3]+1
      if((sample[i,j,1]==1 & sample[i,j,2]==1 & sample[i,j,3]==1)) n[i,4]=n[i,4]+1
      if((sample[i,j,1]==1 & sample[i,j,2]==1 & sample[i,j,3]==0)) n[i,5]=n[i,5]+1
      if((sample[i,j,1]+sample[i,j,2])==1 & sample[i,j,3]==1) n[i,6]=n[i,6]+1
      if((sample[i,j,1]+sample[i,j,2])==1 & sample[i,j,3]==0) n[i,7]=n[i,7]+1
    }
  }
  return(n)
}



#main function for LASSO_FBAT analysis (supporting function are necessary)
LASSO_FBAT <- function(sample,diserisk)
{
  genenum <- dim(sample)[1]
  sampnum <- dim(sample)[2]
  va <- cal_va(sample,genenum)
  
  
  #calculate the mariginal statistical vector
  n <- num(sample,genenum,sampnum)
  Za <- rep(0,genenum)
  for(i in 1:genenum)
  {
    Za[i] <- (n[i,1]-n[i,2])+2*(n[i,3]-n[i,5])+(n[i,6]-n[i,7])
  }
  
  
  #equal weight
  weight_equal=rep(1,genenum)
  Za_equal <- weight_equal%*%Za/sqrt(weight_equal%*%(va*4)%*%weight_equal)
  
  
  #frequency weight
  weight_freq=rep(0,genenum)
  for(i in 1:genenum)
  {
    weight_freq[i]=sum(sample[i,,1])/(4*sampnum)+sum(sample[i,,2])/(4*sampnum)
    weight_freq[i]=1/sqrt(weight_freq[i]*(1-weight_freq[i]))
  }
  Za_freq <- weight_freq%*%Za/sqrt(weight_freq%*%(va*4)%*%weight_freq)
  
  
  #logistic weight
  data.logi <- data.frame(Y=c(diserisk[,1],diserisk[,2],rep(1,sampnum)))
  for(i in 1:genenum)
  {
    new <- matrix(c(sample[i,,1],sample[i,,2],sample[i,,4]),ncol=1)
    newname <- paste0("X",i)
    colnames(new) <- c(newname)
    data.logi <- cbind(data.logi,new)
  }
  
  markname <- paste0("X", 1:genenum)
  fmla <- vector(mode = "character", length = 0)
  fmla <- paste0(markname, collapse = "+")
  formula.X <- as.formula(paste0("Y~",fmla))
  
  logi <- glm(formula.X,family=binomial,data=data.logi)
  weight_logi <- coef(logi)[-1]
  Za_logi<- weight_logi%*%Za/sqrt(weight_logi%*%(va*4)%*%weight_logi)
  if(is.na(Za_logi))
  {
    weight_logi <- weight_equal
    Za_logi <- Za_equal
  }
  
  
  #lasso logistic weight
  X <- matrix(c(sample[1,,1],sample[1,,2],sample[1,,4]),ncol=1)
  for(i in 2:genenum)
  {
    X <- cbind(X,c(sample[i,,1],sample[i,,2],sample[i,,4]))
  }
  Y <- matrix(c(diserisk[,1],diserisk[,2],rep(1,sampnum)),ncol=1)
  cv.fit <- cv.glmnet(X,Y,family='binomial',alpha=1)
  weight_lasso <- as.matrix(coef(cv.fit,s=cv.fit$lambda.min))[-1]
  Za_lasso <- weight_lasso%*%Za/sqrt(weight_lasso%*%(va*4)%*%weight_lasso)
  if(is.na(Za_lasso))
  {
    weight_lasso <- weight_logi
    Za_lasso <- Za_logi
  }
  
  #chi-square test statistic
  S_chi <- Stat_chi(sample,va,genenum)
  
  return(list(T_equal=Za_equal,T_frequency=Za_freq,T_logistic=Za_logi,T_lasso=Za_lasso))
}


LASSO.FBAT.result <- LASSO_FBAT(sample,diserisk)

LASSO.FBAT.result