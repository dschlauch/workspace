## simulations for one condition
#prop: abundance proportions; 
#numspl: the number of samples;
#lss: lower sample scale; 
#uss: upper sample scale; 
#byss: increment unit of sample scales
#lvlm.coef: the estimated coefficients of lm for logvariance-logmean relationship
#lvlm.sd: the estimated sigma of lm for logvariance-logmean relationship
#lvlm.range: the range of the ratio of logvariance and logmean
#seed1: set a value to have the same expected counts for simulated data
#seed2: set a value to have the same NB distribution generated count if seed1 has been set a value already 
#cond: for the colnames of samples

SimuOneCond<-function(prop,
                      numspl,
                      lss,uss,byss,
                      lvlm.coef=c(2,1.5), lvlm.sd=1.3, lvlm.range=c(0.5,3.5),
                      seed1=F,seed2=F,cond=NULL)
{
  if (seed1!=FALSE) set.seed(seed1)  
  spls<-sample(seq(lss,uss,by=byss), numspl, replace=FALSE)
  SimuExp<-round(matrix(prop,nrow=length(prop),ncol=1)%*%spls)
  rownames(SimuExp)<-names(prop)
  SimuExp<-subset(SimuExp,apply(SimuExp,1,function(x) !all(x==0))) #the rows with all zeros are deleted
  
  # NB distribution generated counts
  if (seed2!=FALSE) set.seed(seed2)
  SimuCount<-matrix(0,nrow=nrow(SimuExp),ncol=ncol(SimuExp))
  rownames(SimuCount)<-rownames(SimuExp)
  for (i in 1:nrow(SimuCount)) 
    for (j in 1:ncol(SimuCount))
    {
      V=exp(lvlm.coef[1]+lvlm.coef[2]*log(SimuExp[i,j])+rnorm(1,0,lvlm.sd))
      V=min(V,SimuExp[i,j]^lvlm.range[2])
      V=max(V,SimuExp[i,j]^lvlm.range[1])
      if (V>SimuExp[i,j]) 
        r=SimuExp[i,j]^2/(V-SimuExp[i,j]) else r=1e6
      SimuCount[i,j]<-rnegbin(1,SimuExp[i,j],r)
    }
  
  ids2rem<-apply(SimuCount,1,function(x) !all(x==0)) ##the rows with all zeros are deleted
  SimuCount<-SimuCount[ids2rem,]
  SimuExp<-SimuExp[ids2rem,]
  
  if (is.null(cond))
  {
    colnames(SimuCount)<-paste("s",1:numspl,sep="")
    colnames(SimuExp)<-paste("s",1:numspl,sep="")
  } else
  {
    colnames(SimuCount)<-paste("s",1:numspl,"_",cond, sep="")
    colnames(SimuExp)<-paste("s",1:numspl,"_",cond, sep="")
  }
  
  return(list(SimuCount=SimuCount,SimuExp=SimuExp))
}
