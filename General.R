library(mclust)
library(RcppHMM)
library(randtests)
library(e1071)
MclustAnalysis<-function( Observation, g, name, year){
  starttime=proc.time()
  
  
  
#  Dens <- densityMclust ( Observation , G = g, modelName = "V")
  Dens <- densityMclust ( Observation , G=c(2:4))
  #  Qt   <- quantileMclust( Dens, seq(0.1,1,by=0.1))
  
  ng=10
  Qt   <- quantileMclust( Dens, seq(0,1,by=1/ng)[-1])
  a    <- sort ( c( Qt, Observation), index.return=TRUE)
  o    <- match( c(1:ng), a$ix)
  ob   <- o - c(0,o[-ng]) - 1
  GoodnessFit <- chisq.test(ob,p=rep(1/ng,ng))
  print(paste(name,": Goodness of Fit Chi-square test"))
  print(GoodnessFit)
  #  print(summary(GoodnessFit))
  
  png(filename=paste("/Users/liuxinning/Documents/R/General/",name,year,"qqplot.png"))

  ng=length(Observation)%/%10
  Qt   <- quantileMclust( Dens, seq(0,1,by=1/ng)[-1])
  
  qqplot( Observation, Qt, 
          main = paste( "QQ plot for", name, "and Gaussian Mixture"),
          ylab="Fit distribution",
          xlim=c(-10,10),
          ylim=c(-10,10)
  )
  lines(seq(-20,20,by=1),seq(-20,20,by=1),col=rgb(1,0,0))
  dev.off()
  
  fn<-ecdf(Observation)
  cdf<-cdfMclust(Dens)
  png(filename=paste("/Users/liuxinning/Documents/R/General/",name,year,"cdfplot.png"))
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
  
  if (Dens$parameters$variance$modelName=="E"){
    dens$variance = rep(Dens$parameters$variance$sigmasq,3)
  }else{
    dens$variance = Dens$parameters$variance$sigmasq
  }  
  dens$G=Dens$G
#  dens$variance = Dens$parameters$variance$sigmasq
  dens
}

plotModelSimulation.GM <- function(Observation, parameters, g, name,year){

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
  png(filename=paste("/Users/liuxinning/Documents/R/General/",name,year,"GMSimulation.png"))
  plot(ob,type='l',main=paste(name,"Gaussian Mixture Simulation"),xlab="time",ylab="Simulation",ylim=c(-15,15))
  dev.off()
  print(paste(name,": Simulating Gaussian Mixture use",(proc.time()-starttime)[3],"sec."))
  
}
library(readr)
DailyReturnTable <- read_csv("~/Documents/R/DailyReturnPercent.csv")
DailyReturn<-DailyReturnTable[2:7]
library(e1071)
DailyReturnTable40 <- read_csv("~/Documents/R/DailyReturnPercent40y.csv")
DailyReturn40y=DailyReturnTable40[2:7]
name<-names(DailyReturn)
parameter20y<-c()
parameter40y<-c()
g=3

sink("/Users/liuxinning/Documents/R/General/General.output.txt")

for (i in 1:6){
  print("------------------------")
  print(paste("Analysing ",name[i]))
  
  par(mfrow=c(2,2))
#  png(filename=paste("/Users/liuxinning/Documents/R/General/",name[i],"Histogram.png"))
  
  plot(DailyReturn[[i]],type='l', main = paste("20-year daily return of:",name[i]),xlab=paste("1997.11 - 2017.11 "),ylab="return in percent")
  hist(DailyReturn[[i]],breaks=20,freq=FALSE,main="Histogram",xlab=name[i],xlim=c(-15,15))
  
  plot(DailyReturn40y[[i]],type='l', main = paste("40-year daily return of:",name[i]),xlab=paste("1980.03 - 2017.11 "),ylab="return in percent")
  hist(DailyReturn40y[[i]],breaks=20,freq=FALSE,main="Histogram",xlab=name[i],xlim=c(-15,15))

#  dev.off()
  print(paste("20 years: Kurtosis and Skewness:",kurtosis(DailyReturn[[i]]),skewness(DailyReturn[[i]])))
  print(paste("40 years: Kurtosis and Skewness:",kurtosis(DailyReturn40y[[i]]),skewness(DailyReturn40y[[i]])))

  rt<-runs.test(DailyReturn[[i]], plot=TRUE)
  print(paste(name[i],"20 years Wald-Wolfowitz Runs Test:"))
  print(rt)  
    
  rt<-runs.test(DailyReturn40y[[i]], plot=TRUE)
  print(paste(name[i],"40 years Wald-Wolfowitz Runs Test:"))
  print(rt)
  
  par(mfrow=c(1,1))
  t=MclustAnalysis(DailyReturn[[i]], g, name[i],"20y")
  if(t$G!=g){
    g=t$G
  }
  parameter20y[[i]]<-t
  plotModelSimulation.GM(DailyReturn[[i]],parameter20y[[i]],g,name[i],"20y")
  
  
  t=MclustAnalysis(DailyReturn40y[[i]], g, name[i],"40y")
  if(t$G!=g){
    g=t$G
  }  
  parameter40y[[i]]<-t
  plotModelSimulation.GM(DailyReturn40y[[i]],parameter40y[[i]],g,name[i],"40y")
  
}
sink()