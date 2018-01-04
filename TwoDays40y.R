library(mclust)
library(RcppHMM)
library(e1071)
library(randtests)
MclustAnalysis<-function( Observation, g, name){
  starttime=proc.time()
  
  Dens <- densityMclust ( Observation , G = c(2:4))
  
  ng=10
  Qt   <- quantileMclust( Dens, seq(0,1,by=1/ng)[-1])

  a    <- sort ( c( Qt, Observation), index.return=TRUE)
  o    <- match( c(1:10), a$ix)
  ob   <- o - c(0,o[-10]) - 1
  GoodnessFit <- chisq.test(ob,p=rep(0.1,10))
  print(GoodnessFit)

  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"qqplot.png"))

  ng=length(Observation)%/%10
  Qt   <- quantileMclust( Dens, seq(0,1,by=1/ng)[-1])

  qqplot( Observation, Qt,
          main = paste( "QQ plot for", name, "and Gaussian Mixture"),
          ylab="Fit distribution",
          xlim=c(-10,10),
          ylim=c(-10,10)
  )
  dev.off()

  fn<-ecdf(Observation)
  cdf<-cdfMclust(Dens)

  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"cdfplot.png"))

  plot(fn,col=rgb(1,0,0), main=paste("CDF of empirical and Gaussian Mixture for",name))
  lines(cdf,col=rgb(0,0,1))
  legend("bottomright",
         legend = c("Estimate", "Gaussian Mixture"),
         col = c(rgb(1,0,0),rgb(0,0,1)),
         lty=c(1,1),
         bty='n')
  dev.off()
  print(paste("Fitting Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  
  dens<-c()
  dens$pro = Dens$parameters$pro
  dens$mean = as.vector(Dens$parameters$mean)
  if (Dens$parameters$variance$modelName=="E"){
    dens$variance = rep(Dens$parameters$variance$sigmasq,3)
  }else{
    dens$variance = Dens$parameters$variance$sigmasq
  }
  dens$G=Dens$G
  dens
}

plotModelSimulation.GM <- function(Observation, parameters, g, name){
  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"Observation.png"))
  plot(Observation,type='l', main = paste(name,"observation"),ylim=c(-15,15))
  dev.off()
  
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
  
  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"GMSimulation.png"))
  
  plot(ob,type='l',main=paste(name,"Gaussian Mixture Simulation"),xlab="time",ylab="Simulation",ylim=c(-15,15))
  dev.off()
  
  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"GMHistogram.png"))
  
  hist(ob,breaks=20,freq=FALSE,main=paste(name,"GM Simulation Histogram"),xlab=name,xlim=c(-15,15))
  dev.off()
  
  print(paste("GM Simulation, Kurtosis=",kurtosis(ob),"; Skewness=", skewness(ob),sep=""))
  
  rt<-runs.test(ob, plot=TRUE)
  print(paste(name,"GM Simulation Wald-Wolfowitz Runs Test:"))
  print(rt)
  
  
  print(paste("Simulating Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  
}

HMMGMAnalysis<-function(Observation, parameters, g){
  starttime=proc.time()
  
  N=c()
  for (i in 1:(g*g)){
    N[i]=paste("State ",(i+g-1)%/%g,(i-1)%%g+1,sep="")
  }
  A <- matrix(rep(0,g^4),ncol=g*g)
  for (i in 1:g){
    for(j in 0:(g-1)){
      A[i+j*g,(i*g-g+1):(i*g)]=parameters$pro
    }
  }
  Mu<- matrix(rep(parameters$mean,g),ncol=g*g)
  Sigma <- array(rep(parameters$variance,g), dim = c(1,1,length(N)))
  Pi <- rep(parameters$pro, g)/g
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
  print(paste("Fitting Two days HMM Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  

  HMM.cont.univariate
}

plotModelSimulation.HMM.GM <- function(Observation, hmm, g, name){
  #  plot(Observation,type='l', main = paste(name,"observation"))
  
  starttime=proc.time()
  
  n = length(Observation)

  observationSequence <- generateObservations(hmm, n)
  
  ob=observationSequence$Y[1,]
  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"HMMSimulation.png"))
  plot(ob,type='l',main=paste(name,"HMM Gaussian Mixture Simulation"),xlab="time",ylab="Simulation",ylim=c(-15,15))
  dev.off()
  
  
  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"HMMHistogram.png"))
  
  hist(ob,breaks=20,freq=FALSE,main=paste(name,"HMM Simulation Histogram"),xlab=name,xlim=c(-15,15))
  dev.off()
  
  print(paste("HMM Simulation, Kurtosis=",kurtosis(ob),"; Skewness=", skewness(ob),sep=""))
  rt<-runs.test(ob, plot=TRUE)
  print(paste(name,"HMM Simulation Wald-Wolfowitz Runs Test:"))
  print(rt)
  
  print(paste("Simulating HMM Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  
}

plotStateSimulation<-function(Observation,hmm,g,name){
  starttime=proc.time()
  hiddenStates <- viterbi(hmm, matrix(Observation,ncol=length(Observation)) )
  print(paste(name,"last state is",hiddenStates[length(hiddenStates)]))
  mean=hmm$Mu
  sd=sqrt(hmm$Sigma)
  n=length(Observation)
  a<-c()
  for(i in 1:(g*g)){
    a[(n*(i-1)+1):(n*i)]=rnorm(n,mean[i],sd[i])
  }
  a=matrix(a,ncol = (g*g))
  ob=c()
  for(j in 1:n){
    s=as.numeric(strsplit(hiddenStates[j],' ')[[1]][2])
    s1=s  %% 10
    s2=s %/% 10
    ob[j]=a[j,s1+(s2-1)*g]
  }
  
  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"StateSimulation.png"))
  
  plot(ob,type='l',main=paste(name,"State Simulation"),xlab="time",ylab="Simulation",ylim=c(-15,15))
  dev.off()
  
  png(filename=paste("/Users/liuxinning/Documents/R/TwoDayDependent/",name,"StateSimulationHistogram.png"))
  
  hist(ob,breaks=20,freq=FALSE,main=paste(name,"State Simulation Histogram"),xlab=name,xlim=c(-15,15))
  dev.off()
  
  print(paste("State Simulation, Kurtosis=",kurtosis(ob),"; Skewness=", skewness(ob),sep=""))
  rt<-runs.test(ob, plot=TRUE)
  print(paste(name,"State Simulation Wald-Wolfowitz Runs Test:"))
  print(rt)
  
  
  
  print(paste("State Simulation, Kurtosis=",kurtosis(ob),"; Skewness=", skewness(ob),sep=""))
  print(paste("State Simulation use",(proc.time()-starttime)[3],"sec."))
}


library(readr)
DailyReturn <- read_csv("~/Documents/R/DailyReturnPercent40y.csv")[2:7]
#source("/Users/liuxinning/Documents/R/HMMGaussianMixturegg.R")
parameter<-c()
name<-names(DailyReturn)
hmm<-c()
g=3

sink("/Users/liuxinning/Documents/R/TwoDayDependent/TwoDayDependent.output.txt")
for (i in 1:6) {
 
  print(paste("Analysing ",name[i]))
  
  t=MclustAnalysis(DailyReturn[[i]], g, name[i])
  if(g!=t$G){
    g=t$G
  }
  
  parameter[[i]]<-t
  
  print("Gaussian Mixture parameter:")
  print(parameter[[i]])
  
#  plotModelSimulation.GM(DailyReturn[[i]],parameter[[i]],g,name[i])
  
  hmm[[i]]<-HMMGMAnalysis(DailyReturn[[i]],parameter[[i]],g)
  
  hmmp<-c()
  hmmp$StateNames=hmm[[i]]$StateNames
  hmmp$A=round(hmm[[i]]$A,digits = 7)
  hmmp$Mu=as.vector(hmm[[i]]$Mu)
  hmmp$Sigma=as.vector(hmm[[i]]$Sigma)
  hmmp$Pi=round(as.vector(hmm[[i]]$Pi),digits = 7)
  
  print("HMM parameter:")
  print(hmmp)
  
  plotModelSimulation.HMM.GM(DailyReturn[[i]], hmm[[i]], g, name[i])
  
  plotStateSimulation(DailyReturn[[i]], hmm[[i]], g, name[i])
    
}

sink()


