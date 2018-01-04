library(mclust)
library(RcppHMM)
MclustAnalysis<-function( Observation, g, name){
  starttime=proc.time()
  
  Dens <- densityMclust ( Observation , G = g, modelName = "V")
  ng=10
  Qt   <- quantileMclust( Dens, seq(0,1,by=1/ng)[-1])
#  Qt <- cdfMclust(Dens,Observation,ngrid = ng)$x
  a    <- sort ( c( Qt, Observation), index.return=TRUE)
  o    <- match( c(1:ng), a$ix)
  ob   <- o - c(0,o[-ng]) - 1
  GoodnessFit <- chisq.test(ob,p=rep(1/ng,ng))
  print(paste(name,": Goodness of Fit Chi-square test"))
  print(GoodnessFit)
  #  print(summary(GoodnessFit))
  
  png(filename=paste("/Users/liuxinning/Documents/R/Run20y/",name,"qqplot.png"))
  
  qqplot( Observation, quantileMclust(Dens ,ppoints(length(Observation)) ), 
          main = paste( "QQ plot for", name, "and Gaussian Mixture"),
          ylab="Fit distribution"
  )
  lines(seq(-20,20,by=1),seq(-20,20,by=1),col=rgb(1,0,0))
  dev.off()
  
  fn<-ecdf(Observation)
  cdf<-cdfMclust(Dens)
  
  png(filename=paste("/Users/liuxinning/Documents/R/Run20y/",name,"cdfplot.png"))
  plot(fn,col=rgb(1,0,0), main=paste("CDF of empirical and Gaussian Mixture for",name))
  lines(cdf,col=rgb(0,0,1))
  legend("bottomright", 
         legend = c("Estimate", "Gaussian Mixture"), 
         col = c(rgb(1,0,0),rgb(0,0,1)), 
         lty=c(1,1),
         bty='n')
  dev.off()
  
  #  list(Dens$parameters,GoodnessFit)
  #  Dens$parameters
  print(paste(name,": Fitting Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  
  
  dens<-c()
  dens$pro = Dens$parameters$pro
  dens$mean = as.vector(Dens$parameters$mean)
  dens$variance = Dens$parameters$variance$sigmasq
  dens
}

plotModelSimulation.GM <- function(Observation, parameters, g, name){
  plot(Observation,type='l', main = paste(name,"observation"),ylim=c(-15,15))
  
  starttime=proc.time()
  
  n = length(Observation)
  a=c()
  for(i in 1:g){
    mean=parameters$mean[i]
    sd=sqrt(parameters$variance[i])
    a[(n*(i-1)+1):(n*i)]=rnorm(n,mean,sd)
  }
  a=matrix(a,ncol = g)
  s=sample(g,size=n,replace = TRUE, prob = parameters$pro)
  ob=c()
  for(j in 1:n){
    ob[j]=a[j,s[j]]
  }
  
  png(filename=paste("/Users/liuxinning/Documents/R/Run20y/",name,"GMSimulation.png"))
  plot(ob,type='l',main=paste(name,"Gaussian Mixture Simulation"),xlab="time",ylab="Simulation",ylim=c(-15,15))
  dev.off()
  
  print(paste(name,": Simulating Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  
}
HMMGMAnalysis<-function(Observation, parameters, g){
  starttime=proc.time()
  
  N=c()
  for (i in 1:g){
    N[i]=paste("State",i)
  }
  A <- matrix(rep(parameters$pro,each=g),ncol=g)
  Mu<- matrix(parameters$mean,ncol=g)
  Sigma <- array(parameters$variance, dim = c(1,1,length(N)))
#  Pi <- rep(1/g, g)
  Pi<-parameters$pro
  HMM.cont.univariate <- verifyModel(list( "Model"="GHMM",
                                           "StateNames" = N,
                                           "A" = A,
                                           "Mu" = Mu,
                                           "Sigma" = Sigma,
                                           "Pi" = Pi))
  observationSequences = array(Observation,c(1,length(Observation),1))
  HMM.cont.univariate <- learnEM(HMM.cont.univariate,
                                 observationSequences,
                                 iter= 5000,
                                 delta = 1E-10,
                                 print = FALSE)
  #  print(HMM.cont.univariate)
  print(paste("Fitting HMM Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  
  # hmm<-c()
  # hmm$StateNames=HMM.cont.univariate$StateNames
  # hmm$A=HMM.cont.univariate$A
  # hmm$Mu=as.vector(HMM.cont.univariate$Mu)
  # hmm$Sigma=as.vector(HMM.cont.univariate$Sigma)
  # hmm$Pi=round(as.vector(HMM.cont.univariate$Pi),digits = 5)
  HMM.cont.univariate
}

plotModelSimulation.HMM.GM <- function(Observation, hmm, g, name){
  #  plot(Observation,type='l', main = paste(name,"observation"))
  
  starttime=proc.time()
  
  n = length(Observation)
  # s=c(sample(g,size=1,replace = TRUE, prob = hmm$Pi))
  # for(i in 2:n){
  #   s[i] = sample(g,1,TRUE, prob=hmm$A[s[i-1],])
  # }
  # a<-c()
  # for(i in 1:g){
  #   mean=hmm$Mu[i]
  #   sd=sqrt(hmm$Sigma[i])
  #   a[(n*(i-1)+1):(n*i)]=rnorm(n,mean,sd)
  # }
  # a=matrix(a,ncol = g)
  # ob=c()
  # for(j in 1:n){
  #   ob[j]=a[j,s[j]]
  # }
  # plot(ob,type='l',main=paste(name,"HMM Gaussian Mixture Simulation"),xlab="time",ylab="Simulation")
  
  
  observationSequence <- generateObservations(hmm, n)
  png(filename=paste("/Users/liuxinning/Documents/R/Run20y/",name,"HMMSimulation.png"))
  
  plot(observationSequence$Y[1,],type='l',main=paste(name,"HMM Gaussian Mixture Simulation"),xlab="time",ylab="Simulation",ylim=c(-15,15))
  dev.off()
  
  print(paste("Simulating HMM Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  
}

plotModelSimulation.HMM.useparameterGM <- function(Observation, hmm, parameters, g, name){
  
  starttime=proc.time()
  
  n = length(Observation)
  s=c(sample(g,size=1,replace = TRUE, prob = hmm$Pi))
  for(i in 2:n){
    s[i] = sample(g,1,TRUE, prob=hmm$A[s[i-1],])
  }
  a<-c()
  for(i in 1:g){
    mean=parameters$mean[i]
    sd=sqrt(parameters$variance[i])
    a[(n*(i-1)+1):(n*i)]=rnorm(n,mean,sd)
  }
  a=matrix(a,ncol = g)
  ob=c()
  for(j in 1:n){
    ob[j]=a[j,s[j]]
  }
  png(filename=paste("/Users/liuxinning/Documents/R/Run20y/",name,"HMMwithGMSimulation.png"))
  
  plot(ob,type='l',main=paste(name,"HMM using Gaussian Mixture parameter"),xlab="time",ylab="Simulation",ylim=c(-15,15))
  dev.off()
  print(paste("Simulating HMM using Gaussian Mixture parameter use",(proc.time()-starttime)[3],"sec."))
  
}

plotStateSimulation<-function(Observation,hmm,g,name){
  starttime=proc.time()
  hiddenStates <- viterbi(hmm, matrix(Observation,ncol=length(Observation)) )
  mean=hmm$Mu
  sd=sqrt(hmm$Sigma)
  n=length(Observation)
  a<-c()
  for(i in 1:g){
    a[(n*(i-1)+1):(n*i)]=rnorm(n,mean[i],sd[i])
  }
  a=matrix(a,ncol = g)
  ob=c()
  for(j in 1:n){
    s=as.numeric(strsplit(hiddenStates[j],' ')[[1]][2]) %% 10
    ob[j]=a[j,s]
  }
  png(filename=paste("/Users/liuxinning/Documents/R/Run20y/",name,"StateSimulation.png"))
  
  plot(ob,type='l',main=paste(name,"State Simulation"),xlab="time",ylab="Simulation")
  dev.off()
  print(paste("State Simulation use",(proc.time()-starttime)[3],"sec."))
}

library(readr)
#library(doParallel)
#library(foreach)
DailyReturn <- read_csv("~/Documents/R/DailyReturnPercent.csv")[2:7]
#source("/Users/liuxinning/Documents/R/HMMGaussianMixture.R")
parameter<-c()
name<-names(DailyReturn)
hmm<-c()
g=3

#no_cores <- detectCores()
#cl<-makeCluster(no_cores)
#registerDoParallel(cl)

sink("/Users/liuxinning/Documents/R/Run20y/Run6stocks20y.output.txt")

for(i in 1:6){  
#foreach (i = 1:6, .packages = c("mclust","RcppHMM")) %dopar% {
  print(paste("Analysing ",name[i]))
  
  t=MclustAnalysis(DailyReturn[[i]], g, name[i])
  
  parameter[[i]]<-t
  
  print("Gaussian Mixture parameter:")
  print(parameter[[i]])
  
  plotModelSimulation.GM(DailyReturn[[i]],parameter[[i]],g,name[i])
  
  hmm[[i]]<-HMMGMAnalysis(DailyReturn[[i]],parameter[[i]],g)
  
  hmmp<-c()
  hmmp$StateNames=hmm[[i]]$StateNames
  hmmp$A=round(hmm[[i]]$A,digits=7)
  hmmp$Mu=as.vector(hmm[[i]]$Mu)
  hmmp$Sigma=as.vector(hmm[[i]]$Sigma)
  hmmp$Pi=round(as.vector(hmm[[i]]$Pi),digits = 7)
  
  
  print("HMM parameter:")
  print(hmmp)
  
  plotModelSimulation.HMM.GM(DailyReturn[[i]], hmm[[i]], g, name[i])
  
  plotModelSimulation.HMM.useparameterGM(DailyReturn[[i]], hmm[[i]], parameter[[i]], g, name[i])
  
  plotStateSimulation(DailyReturn[[i]], hmm[[i]], g, name[i])
}

#stopCluster(cl)
sink()